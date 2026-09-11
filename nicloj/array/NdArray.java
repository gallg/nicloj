package nicloj.array;

import java.util.Arrays;
import nicloj.header.DataType;
import java.util.Locale;

/**
 * A dense n-dimensional array of doubles in column-major order, matching the
 * NIfTI on-disk layout where the first axis varies fastest.
 *
 * <p>Instances are treated as immutable by the rest of nicloj; the mutating
 * setters exist for building arrays before handing them off.
 */
public final class NdArray {
    private final int[] shape;
    private final int[] strides;
    private final Store data;

    public NdArray(int[] shape) {
        this(shape, new Store.Doubles(checkedSize(shape)));
    }

    public NdArray(int[] shape, double[] data) {
        this(shape, new Store.Doubles(data));
    }

    public NdArray(int[] shape, Store data) {
        int n = checkedSize(shape);
        if (data.size() != n) {
            throw new IllegalArgumentException(
                    "data length " + data.size() + " does not match shape " + Arrays.toString(shape));
        }
        this.shape = shape.clone();
        this.data = data;
        this.strides = columnMajorStrides(this.shape);
    }

    /** Wrap {@code data} without copying; caller must not mutate it afterwards. */
    public static NdArray wrap(int[] shape, double[] data) { return new NdArray(shape, data); }

    /** Wrap a store without copying; caller must not mutate it afterwards. */
    public static NdArray wrap(int[] shape, Store data) { return new NdArray(shape, data); }


    public int[] shape() { return shape.clone(); }
    public int ndim() { return shape.length; }
    public int size() { return data.size(); }
    public int dim(int axis) { return shape[axis]; }

    /** The backing store, in column-major order. Not copied. */
    public Store store() { return data; }

    /** The voxel type the elements are held in. */
    public DataType dtype() { return data.type(); }

    /** A fresh {@code double[]} of every element, in column-major order. */
    public double[] toDoubleArray() {
        double[] out = new double[data.size()];
        for (int i = 0; i < out.length; i++) out[i] = data.get(i);
        return out;
    }

    public double flat(int i) { return data.get(i); }
    public void setFlat(int i, double v) { data.set(i, v); }

    public double get(int... idx) { return data.get(offset(idx)); }
    public void set(double v, int... idx) { data.set(offset(idx), v); }

    public int offset(int[] idx) {
        if (idx.length != shape.length) {
            throw new IllegalArgumentException(
                    "expected " + shape.length + " indices, got " + idx.length);
        }
        int off = 0;
        for (int i = 0; i < idx.length; i++) {
            int k = idx[i];
            if (k < 0 || k >= shape[i]) {
                throw new IndexOutOfBoundsException(
                        "index " + k + " out of bounds for axis " + i + " of size " + shape[i]);
            }
            off += k * strides[i];
        }
        return off;
    }


    // -------------------------------------------------------------- reshaping

    /** Reinterpret the same buffer under a new shape of equal total size. */
    public NdArray reshape(int[] newShape) {
        if (checkedSize(newShape) != data.size()) {
            throw new IllegalArgumentException(
                    "cannot reshape size " + data.size() + " into " + Arrays.toString(newShape));
        }
        return new NdArray(newShape, data);
    }

    /** Drop every axis of length 1. */
    public NdArray squeeze() {
        int n = 0;
        for (int s : shape) if (s != 1) n++;
        int[] out = new int[n];
        int j = 0;
        for (int s : shape) if (s != 1) out[j++] = s;
        return reshape(out);
    }

    /** Permute axes; {@code perm[i]} is the source axis for output axis {@code i}. */
    public NdArray transpose(int[] perm) {
        if (perm.length != shape.length) {
            throw new IllegalArgumentException("permutation length must equal ndim");
        }
        int[] outShape = new int[perm.length];
        for (int i = 0; i < perm.length; i++) outShape[i] = shape[perm[i]];
        Store out = data.alloc(data.size());
        int[] outStrides = columnMajorStrides(outShape);
        int[] idx = new int[perm.length];
        for (int o = 0; o < out.size(); o++) {
            int src = 0;
            for (int i = 0; i < perm.length; i++) {
                idx[i] = o / outStrides[i] % outShape[i];
                src += idx[i] * strides[perm[i]];
            }
            copy(data, src, out, o);
        }
        return new NdArray(outShape, out);
    }

    /** Reverse the order of elements along {@code axis}. */
    public NdArray flip(int axis) {
        Store out = data.alloc(data.size());
        int n = shape[axis];
        int stride = strides[axis];
        for (int o = 0; o < out.size(); o++) {
            int k = o / stride % n;
            copy(data, o + (n - 1 - 2 * k) * stride, out, o);
        }
        return new NdArray(shape, out);
    }

    /**
     * Extract a strided sub-block. {@code start} is inclusive, {@code stop}
     * exclusive, and {@code step} must be positive.
     */
    public NdArray slice(int[] start, int[] stop, int[] step) {
        int nd = shape.length;
        if (start.length != nd || stop.length != nd || step.length != nd) {
            throw new IllegalArgumentException("start, stop and step must have length " + nd);
        }
        int[] outShape = new int[nd];
        for (int i = 0; i < nd; i++) {
            if (step[i] <= 0) throw new IllegalArgumentException("step must be positive");
            if (start[i] < 0 || stop[i] > shape[i] || start[i] > stop[i]) {
                throw new IndexOutOfBoundsException("bad slice on axis " + i);
            }
            outShape[i] = (stop[i] - start[i] + step[i] - 1) / step[i];
        }
        Store out = data.alloc(checkedSize(outShape));
        int[] outStrides = columnMajorStrides(outShape);
        for (int o = 0; o < out.size(); o++) {
            int src = 0;
            for (int i = 0; i < nd; i++) {
                src += (start[i] + step[i] * (o / outStrides[i] % outShape[i])) * strides[i];
            }
            copy(data, src, out, o);
        }
        return new NdArray(outShape, out);
    }

    /** Join arrays of matching shape (except along {@code axis}) end to end. */
    public static NdArray concat(NdArray[] parts, int axis) {
        if (parts.length == 0) throw new IllegalArgumentException("nothing to concatenate");
        int nd = parts[0].ndim();
        int total = 0;
        for (NdArray p : parts) {
            if (p.ndim() != nd) throw new IllegalArgumentException("mismatched ndim");
            for (int i = 0; i < nd; i++) {
                if (i != axis && p.shape[i] != parts[0].shape[i]) {
                    throw new IllegalArgumentException(
                            "mismatched shape on axis " + i + ": " + Arrays.toString(p.shape));
                }
            }
            total += p.shape[axis];
        }
        int[] outShape = parts[0].shape();
        outShape[axis] = total;
        // Keep the voxel type when every part agrees on it, widen otherwise.
        Store store = parts[0].data.alloc(checkedSize(outShape));
        for (NdArray p : parts) {
            if (p.dtype() != parts[0].dtype()) {
                store = new Store.Doubles(checkedSize(outShape));
                break;
            }
        }
        NdArray out = new NdArray(outShape, store);
        int[] outStrides = columnMajorStrides(outShape);
        int base = 0;
        for (NdArray p : parts) {
            for (int i = 0; i < p.size(); i++) {
                int off = 0;
                for (int a = 0; a < nd; a++) {
                    int k = i / p.strides[a] % p.shape[a];
                    off += (a == axis ? k + base : k) * outStrides[a];
                }
                copy(p.data, i, out.data, off);
            }
            base += p.shape[axis];
        }
        return out;
    }

    // ------------------------------------------------------------- statistics

    /** Minimum and maximum over finite entries, or {@code {+inf, -inf}} if none. */
    public double[] finiteRange() {
        double mn = Double.POSITIVE_INFINITY;
        double mx = Double.NEGATIVE_INFINITY;
        boolean check = data.canBeNonFinite();
        for (int i = 0; i < data.size(); i++) {
            double v = data.get(i);
            if (check && !Double.isFinite(v)) continue;
            if (v < mn) mn = v;
            if (v > mx) mx = v;
        }
        return new double[] {mn, mx};
    }

    public boolean hasNaN() {
        if (!data.canBeNonFinite()) return false;
        for (int i = 0; i < data.size(); i++) if (Double.isNaN(data.get(i))) return true;
        return false;
    }

    /** True when every finite entry is a whole number. */
    public boolean isIntegral() {
        if (data.exact()) return true;
        for (int i = 0; i < data.size(); i++) {
            double v = data.get(i);
            if (Double.isFinite(v) && v != Math.rint(v)) return false;
        }
        return true;
    }

    /**
     * Elementwise {@code v * slope + inter}. Identity scaling returns this array
     * untouched, which is what keeps {@code fdata} free on an unscaled image.
     */
    public NdArray scaled(double slope, double inter) {
        if (slope == 1.0 && inter == 0.0) return this;
        double[] out = new double[data.size()];
        for (int i = 0; i < out.length; i++) out[i] = data.get(i) * slope + inter;
        return new NdArray(shape, out);
    }

    public boolean closeTo(NdArray other, double atol) {
        if (!Arrays.equals(shape, other.shape)) return false;
        for (int i = 0; i < data.size(); i++) {
            double a = data.get(i);
            double b = other.data.get(i);
            if (Double.isNaN(a) && Double.isNaN(b)) continue;
            if (!(Math.abs(a - b) <= atol)) return false;
        }
        return true;
    }

    /**
     * Move one element between stores, exactly when both hold integers, so a
     * reshape or a transpose cannot round an {@code int64} past 2^53.
     */
    private static void copy(Store from, int i, Store to, int j) {
        if (from.exact() && to.exact()) to.setLong(j, from.getLong(i));
        else to.set(j, from.get(i));
    }

    /** Elements shown at each end of an abbreviated axis, as in numpy. */
    private static final int EDGE = 3;

    /** Arrays larger than this are abbreviated rather than printed whole. */
    private static final int THRESHOLD = 1000;

    /**
     * A numpy-style repr: nested brackets one level per axis, and for anything
     * over {@link #THRESHOLD} elements only the first and last {@link #EDGE}
     * along each long axis, with {@code ...} standing for the rest.
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder("#nicloj/ndarray[");
        for (int i = 0; i < shape.length; i++) {
            if (i > 0) sb.append(' ');
            sb.append(shape[i]);
        }
        return sb.append("]\n").append(layout()).toString();
    }

    /**
     * The bracketed element layout alone, without the shape header, for callers
     * that print a view onto this array under their own name. Only the elements
     * actually shown are visited, so this stays cheap on a large array.
     */
    public String layout() {
        boolean brief = data.size() > THRESHOLD;
        double[] shown = new double[shownCount(brief)];
        collect(new int[shape.length], 0, brief, shown, new int[1]);
        String fmt = numberFormat(shown);
        int width = 0;
        for (double v : shown) width = Math.max(width, String.format(Locale.ROOT, fmt, v).length());

        StringBuilder sb = new StringBuilder();
        render(sb, new int[shape.length], 0, brief, fmt, width);
        return sb.toString();
    }

    /** Indices to print along {@code axis}; -1 marks the elided middle. */
    private int[] axisIndices(int axis, boolean brief) {
        int n = shape[axis];
        if (!brief || n <= 2 * EDGE) {
            int[] all = new int[n];
            for (int i = 0; i < n; i++) all[i] = i;
            return all;
        }
        int[] ends = new int[2 * EDGE + 1];
        for (int i = 0; i < EDGE; i++) {
            ends[i] = i;
            ends[EDGE + 1 + i] = n - EDGE + i;
        }
        ends[EDGE] = -1;
        return ends;
    }

    private int shownCount(boolean brief) {
        int n = 1;
        for (int dim : shape) n *= (brief && dim > 2 * EDGE) ? 2 * EDGE : dim;
        return n;
    }

    private void collect(int[] idx, int axis, boolean brief, double[] out, int[] at) {
        if (axis == shape.length) {
            out[at[0]++] = data.get(offset(idx));
            return;
        }
        for (int i : axisIndices(axis, brief)) {
            if (i < 0) continue;
            idx[axis] = i;
            collect(idx, axis + 1, brief, out, at);
        }
    }

    /** Integers print as {@code 6195.}, anything else to six significant digits. */
    private static String numberFormat(double[] shown) {
        for (double v : shown) {
            if (!Double.isFinite(v) || v != Math.rint(v) || Math.abs(v) >= 1e15) return "%.6g";
        }
        return "%.0f.";
    }

    private void render(StringBuilder sb, int[] idx, int axis, boolean brief, String fmt, int width) {
        if (axis == shape.length) {
            String s = String.format(Locale.ROOT, fmt, data.get(offset(idx)));
            for (int p = s.length(); p < width; p++) sb.append(' ');
            sb.append(s);
            return;
        }
        sb.append('[');
        int[] ix = axisIndices(axis, brief);
        for (int k = 0; k < ix.length; k++) {
            if (k > 0) {
                if (axis == shape.length - 1) {
                    sb.append(' ');
                } else {
                    for (int n = shape.length - 1 - axis; n > 0; n--) sb.append('\n');
                    for (int p = 0; p <= axis; p++) sb.append(' ');
                }
            }
            if (ix[k] < 0) {
                sb.append("...");
            } else {
                idx[axis] = ix[k];
                render(sb, idx, axis + 1, brief, fmt, width);
            }
        }
        sb.append(']');
    }

    // ----------------------------------------------------------------- helpers

    /** Column-major strides: {@code stride[0] == 1}. */
    static int[] columnMajorStrides(int[] shape) {
        int[] s = new int[shape.length];
        int acc = 1;
        for (int i = 0; i < shape.length; i++) {
            s[i] = acc;
            acc *= shape[i];
        }
        return s;
    }

    public static int checkedSize(int[] shape) {
        long n = 1;
        for (int s : shape) {
            if (s < 0) throw new IllegalArgumentException("negative dimension " + s);
            n *= s;
            if (n > Integer.MAX_VALUE) {
                throw new IllegalArgumentException("array too large: " + Arrays.toString(shape));
            }
        }
        return (int) n;
    }
}
