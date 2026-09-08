package nicloj.affine;

import java.util.Arrays;
import nicloj.header.NiftiError;

/**
 * Axis orientation bookkeeping.
 *
 * <p>An <em>orientation</em> is a {@code (p, 2)} table with one row per input
 * (voxel) axis: column 0 is the closest output (world) axis, column 1 is
 * {@code +1} or {@code -1} for direction. A row of {@code NaN} marks an input
 * axis with no corresponding output axis.
 *
 * <p>Ported from nibabel's {@code nibabel.orientations} (MIT).
 */
public final class Orientations {
    private Orientations() {}

    /** Default axis labels: negative then positive end of each output axis. */
    public static final String[][] RAS_LABELS = {{"L", "R"}, {"P", "A"}, {"I", "S"}};

    /**
     * Orientation of the input axes of {@code affine} in terms of its output
     * axes. Pass a non-positive {@code tol} to use numpy's default rank cutoff.
     */
    public static double[][] ioOrientation(double[][] affine, double tol) {
        int q = affine.length - 1;
        int p = affine[0].length - 1;
        double[][] rzs = Mat.block(affine, q, p);
        double[] zooms = Mat.columnNorms(rzs);
        double[][] rs = new double[q][p];
        for (int i = 0; i < q; i++) {
            for (int j = 0; j < p; j++) rs[i][j] = rzs[i][j] / (zooms[j] == 0 ? 1 : zooms[j]);
        }
        Svd svd = Svd.of(rs);
        double[][] r = svd.polar(tol > 0 ? tol : svd.defaultTol(q, p));

        double[][] ornt = new double[p][2];
        for (double[] row : ornt) Arrays.fill(row, Double.NaN);
        // Resolve the strongest input axes first so ties break consistently.
        for (int in : byDescendingWeight(r, p)) {
            double[] col = Mat.column(r, in);
            int out = -1;
            double best = 0;
            for (int i = 0; i < col.length; i++) {
                if (Math.abs(col[i]) > best) { best = Math.abs(col[i]); out = i; }
            }
            if (out < 0 || best == 0.0) continue;
            ornt[in][0] = out;
            ornt[in][1] = col[out] < 0 ? -1 : 1;
            // Retire the claimed output axis.
            Arrays.fill(r[out], 0.0);
        }
        return ornt;
    }

    private static Integer[] byDescendingWeight(double[][] r, int p) {
        double[] weight = new double[p];
        for (int j = 0; j < p; j++) {
            for (double[] row : r) weight[j] = Math.max(weight[j], row[j] * row[j]);
        }
        Integer[] order = new Integer[p];
        for (int j = 0; j < p; j++) order[j] = j;
        Arrays.sort(order, (a, b) -> Double.compare(weight[b], weight[a]));
        return order;
    }

    /** Axis labels for an orientation; {@code null} entries mark dropped axes. */
    public static String[] toAxcodes(double[][] ornt, String[][] labels) {
        String[] out = new String[ornt.length];
        for (int i = 0; i < ornt.length; i++) {
            double axis = ornt[i][0];
            double dir = ornt[i][1];
            if (Double.isNaN(axis)) continue;
            int a = (int) Math.round(axis);
            if (a != axis) throw new NiftiError("non-integer axis number " + axis);
            if (dir == 1) out[i] = labels[a][1];
            else if (dir == -1) out[i] = labels[a][0];
            else throw new NiftiError("direction should be -1 or 1, got " + dir);
        }
        return out;
    }

    /** Orientation table for a list of axis labels such as {@code {"R","A","S"}}. */
    public static double[][] fromAxcodes(String[] axcodes, String[][] labels) {
        double[][] ornt = new double[axcodes.length][2];
        for (double[] row : ornt) Arrays.fill(row, Double.NaN);
        for (int i = 0; i < axcodes.length; i++) {
            String code = axcodes[i];
            if (code == null) continue;
            boolean found = false;
            for (int a = 0; a < labels.length && !found; a++) {
                for (int d = 0; d < 2; d++) {
                    if (labels[a][d].equalsIgnoreCase(code)) {
                        ornt[i][0] = a;
                        ornt[i][1] = d == 0 ? -1 : 1;
                        found = true;
                        break;
                    }
                }
            }
            if (!found) throw new NiftiError("unknown axis code '" + code + "'");
        }
        return ornt;
    }

    /** Axis labels of an affine's input axes, i.e. {@code ioOrientation} then {@code toAxcodes}. */
    public static String[] axcodes(double[][] affine) {
        return toAxcodes(ioOrientation(affine, 0), RAS_LABELS);
    }

    /** The orientation taking data from {@code start} to {@code end}. */
    public static double[][] transform(double[][] start, double[][] end) {
        if (start.length != end.length) {
            throw new NiftiError("orientations must have the same length");
        }
        double[][] result = new double[start.length][2];
        for (int endIn = 0; endIn < end.length; endIn++) {
            boolean matched = false;
            for (int startIn = 0; startIn < start.length; startIn++) {
                if (start[startIn][0] == end[endIn][0]) {
                    result[startIn][0] = endIn;
                    result[startIn][1] = start[startIn][1] == end[endIn][1] ? 1 : -1;
                    matched = true;
                    break;
                }
            }
            if (!matched) {
                throw new NiftiError("output axis " + end[endIn][0] + " missing from start orientation");
            }
        }
        return result;
    }

    /**
     * Affine that undoes the flips and transpose described by {@code ornt} on an
     * array of {@code shape}, so that
     * {@code new_affine = old_affine * invOrntAff(ornt, shape)}.
     */
    public static double[][] invOrntAff(double[][] ornt, int[] shape) {
        int p = ornt.length;
        for (double[] row : ornt) {
            if (Double.isNaN(row[0]) || Double.isNaN(row[1])) {
                throw new NiftiError("cannot invert an orientation with dropped axes");
            }
        }
        double[][] undoReorder = new double[p + 1][p + 1];
        for (int i = 0; i < p; i++) undoReorder[i][(int) ornt[i][0]] = 1.0;
        undoReorder[p][p] = 1.0;

        double[][] undoFlip = Mat.identity(p + 1);
        for (int i = 0; i < p; i++) {
            double flip = ornt[i][1];
            double centre = -(shape[i] - 1) / 2.0;
            undoFlip[i][i] = flip;
            undoFlip[i][p] = flip * centre - centre;
        }
        return Mat.mul(undoFlip, undoReorder);
    }

    /** Inverse permutation, i.e. the transpose that {@code applyOrientation} performs. */
    public static int[] axisPermutation(double[][] ornt) {
        int p = ornt.length;
        int[] perm = new int[p];
        for (int i = 0; i < p; i++) {
            double axis = ornt[i][0];
            if (Double.isNaN(axis)) {
                throw new NiftiError("cannot reorder data with dropped axes");
            }
            perm[(int) axis] = i;
        }
        return perm;
    }
}
