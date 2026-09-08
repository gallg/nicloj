package nicloj.array;

import java.util.Arrays;

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
    private final double[] data;

    public NdArray(int[] shape) {
        this(shape, new double[checkedSize(shape)]);
    }

    public NdArray(int[] shape, double[] data) {
        int n = checkedSize(shape);
        if (data.length != n) {
            throw new IllegalArgumentException(
                    "data length " + data.length + " does not match shape " + Arrays.toString(shape));
        }
        this.shape = shape.clone();
        this.data = data;
        this.strides = columnMajorStrides(this.shape);
    }

    /** Wrap {@code data} without copying; caller must not mutate it afterwards. */
    public static NdArray wrap(int[] shape, double[] data) { return new NdArray(shape, data); }

    public static NdArray scalar(double v) { return new NdArray(new int[0], new double[] {v}); }

    public int[] shape() { return shape.clone(); }
    public int ndim() { return shape.length; }
    public int size() { return data.length; }
    public int dim(int axis) { return shape[axis]; }

    /** The backing buffer, in column-major order. Not copied. */
    public double[] data() { return data; }

    public double flat(int i) { return data[i]; }
    public void setFlat(int i, double v) { data[i] = v; }

    public double get(int... idx) { return data[offset(idx)]; }
    public void set(double v, int... idx) { data[offset(idx)] = v; }

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

    /** Decompose a flat offset back into per-axis indices. */
    public int[] indices(int offset) {
        int[] idx = new int[shape.length];
        for (int i = 0; i < shape.length; i++) {
            idx[i] = offset / strides[i] % shape[i];
        }
        return idx;
    }

    // -------------------------------------------------------------- reshaping

    /** Reinterpret the same buffer under a new shape of equal total size. */
    public NdArray reshape(int[] newShape) {
        if (checkedSize(newShape) != data.length) {
            throw new IllegalArgumentException(
                    "cannot reshape size " + data.length + " into " + Arrays.toString(newShape));
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
        double[] out = new double[data.length];
        int[] outStrides = columnMajorStrides(outShape);
        int[] idx = new int[perm.length];
        for (int o = 0; o < out.length; o++) {
            int src = 0;
            for (int i = 0; i < perm.length; i++) {
                idx[i] = o / outStrides[i] % outShape[i];
                src += idx[i] * strides[perm[i]];
            }
            out[o] = data[src];
        }
        return new NdArray(outShape, out);
    }

    /** Reverse the order of elements along {@code axis}. */
    public NdArray flip(int axis) {
        double[] out = new double[data.length];
        int n = shape[axis];
        int stride = strides[axis];
        for (int o = 0; o < out.length; o++) {
            int k = o / stride % n;
            out[o] = data[o + (n - 1 - 2 * k) * stride];
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
        double[] out = new double[checkedSize(outShape)];
        int[] outStrides = columnMajorStrides(outShape);
        for (int o = 0; o < out.length; o++) {
            int src = 0;
            for (int i = 0; i < nd; i++) {
                src += (start[i] + step[i] * (o / outStrides[i] % outShape[i])) * strides[i];
            }
            out[o] = data[src];
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
        NdArray out = new NdArray(outShape);
        int[] outStrides = columnMajorStrides(outShape);
        int base = 0;
        for (NdArray p : parts) {
            for (int i = 0; i < p.size(); i++) {
                int off = 0;
                for (int a = 0; a < nd; a++) {
                    int k = i / p.strides[a] % p.shape[a];
                    off += (a == axis ? k + base : k) * outStrides[a];
                }
                out.data[off] = p.data[i];
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
        for (double v : data) {
            if (!Double.isFinite(v)) continue;
            if (v < mn) mn = v;
            if (v > mx) mx = v;
        }
        return new double[] {mn, mx};
    }

    public boolean hasNaN() {
        for (double v : data) if (Double.isNaN(v)) return true;
        return false;
    }

    /** True when every finite entry is a whole number. */
    public boolean isIntegral() {
        for (double v : data) {
            if (Double.isFinite(v) && v != Math.rint(v)) return false;
        }
        return true;
    }

    /** Elementwise {@code v * slope + inter}, returning a new array. */
    public NdArray scaled(double slope, double inter) {
        if (slope == 1.0 && inter == 0.0) return this;
        double[] out = new double[data.length];
        for (int i = 0; i < out.length; i++) out[i] = data[i] * slope + inter;
        return new NdArray(shape, out);
    }

    public boolean closeTo(NdArray other, double atol) {
        if (!Arrays.equals(shape, other.shape)) return false;
        for (int i = 0; i < data.length; i++) {
            double a = data[i];
            double b = other.data[i];
            if (Double.isNaN(a) && Double.isNaN(b)) continue;
            if (!(Math.abs(a - b) <= atol)) return false;
        }
        return true;
    }

    @Override
    public String toString() {
        return "NdArray" + Arrays.toString(shape) + "(" + data.length + " values)";
    }

    // ----------------------------------------------------------------- helpers

    /** Column-major strides: {@code stride[0] == 1}. */
    public static int[] columnMajorStrides(int[] shape) {
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
