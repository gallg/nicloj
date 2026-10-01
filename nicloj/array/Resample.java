package nicloj.array;

import nicloj.header.DataType;

/**
 * Resampling an array through an affine map of voxel indices, the way
 * {@code scipy.ndimage.affine_transform} does it with {@code mode='constant'}:
 * a point outside {@code [0, dim - 1]} on any axis takes {@code cval}, and
 * nothing is interpolated past the edge.
 *
 * <p>Unlike scipy, a point within {@link #EDGE} of a voxel past the edge counts
 * as on it. Composing affines leaves rounding residue of that order, and
 * without the slack resampling an image onto its own grid can lose a slice.
 */
public final class Resample {
    private Resample() {}

    static final double EDGE = 1e-9;

    /**
     * Sample {@code in} at {@code m · (i, 1)} for every index {@code i} of an
     * array of {@code outShape}. {@code m} has at least {@code n} rows of
     * {@code n + 1} columns and maps output indices to input ones. Order 0 is
     * nearest neighbour, rounding halves up; order 1 is linear.
     *
     * <p>The result keeps the input's voxel type. Integer types are rounded
     * half away from zero and clamped to their range, as scipy does.
     */
    public static NdArray affine(NdArray in, double[][] m, int[] outShape, int order, double cval) {
        int n = in.ndim();
        if (outShape.length != n) {
            throw new IllegalArgumentException(
                    "output shape needs " + n + " axes, got " + outShape.length);
        }
        if (m.length < n || m[0].length != n + 1) {
            throw new IllegalArgumentException(
                    "need an affine with " + (n + 1) + " columns for " + n + " axes");
        }
        if (order != 0 && order != 1) {
            throw new IllegalArgumentException("order must be 0 or 1, got " + order);
        }
        Store src = in.store();
        DataType type = src.type();
        Store out = src.alloc(NdArray.checkedSize(outShape));
        int[] dims = in.shape();
        int[] strides = NdArray.columnMajorStrides(dims);
        int[] idx = new int[n];
        double[] x = new double[n];
        int[] lo = new int[n];
        double[] frac = new double[n];
        for (int o = 0, size = out.size(); o < size; o++) {
            boolean inside = true;
            for (int j = 0; j < n && inside; j++) {
                double v = m[j][n];
                for (int k = 0; k < n; k++) v += m[j][k] * idx[k];
                inside = v >= -EDGE && v <= dims[j] - 1 + EDGE;   // false for NaN too
                x[j] = Math.min(dims[j] - 1, Math.max(0, v));
            }
            double v;
            if (!inside) {
                v = cval;
            } else if (order == 0) {
                int off = 0;
                for (int j = 0; j < n; j++) off += (int) Math.floor(x[j] + 0.5) * strides[j];
                v = src.get(off);
            } else {
                for (int j = 0; j < n; j++) {
                    lo[j] = (int) Math.floor(x[j]);
                    frac[j] = x[j] - lo[j];
                }
                // Sum over the corners of the enclosing cell, skipping axes the
                // point sits exactly on, which also keeps the last index in range.
                v = 0;
                for (int corner = 0; corner < (1 << n); corner++) {
                    double w = 1;
                    int off = 0;
                    for (int j = 0; j < n && w != 0; j++) {
                        boolean up = (corner & (1 << j)) != 0;
                        if (up && frac[j] == 0) w = 0;
                        w *= up ? frac[j] : 1 - frac[j];
                        off += (lo[j] + (up ? 1 : 0)) * strides[j];
                    }
                    if (w != 0) v += w * src.get(off);
                }
            }
            out.set(o, type.isInteger() ? fitInteger(v, type) : v);
            for (int a = 0; a < n; a++) {
                if (++idx[a] < outShape[a]) break;
                idx[a] = 0;
            }
        }
        return new NdArray(outShape, out);
    }

    private static double fitInteger(double v, DataType type) {
        double r = Math.copySign(Math.floor(Math.abs(v) + 0.5), v);
        return Math.min(type.maxValue(), Math.max(type.minValue(), r));
    }
}
