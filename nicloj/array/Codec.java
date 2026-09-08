package nicloj.array;

import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import nicloj.header.DataType;
import nicloj.header.NiftiError;

/** Conversion between raw NIfTI voxel bytes and {@link NdArray}. */
public final class Codec {
    private Codec() {}

    private static final double TWO_POW_64 = 0x1p64;
    private static final double TWO_POW_63 = 0x1p63;

    /**
     * Decode {@code raw} into an array of {@code shape}, applying
     * {@code v * slope + inter}. Pass slope 1 / inter 0 for unscaled values.
     */
    public static NdArray decode(byte[] raw, int[] shape, DataType type, ByteOrder order,
                                 double slope, double inter) {
        int n = NdArray.checkedSize(shape);
        long need = (long) n * type.itemSize();
        if (raw.length < need) {
            throw new NiftiError("voxel data truncated: need " + need + " bytes, have " + raw.length);
        }
        if (!type.isReal()) {
            throw new NiftiError(type.label() + " voxel data is not supported yet");
        }
        double[] out = new double[n];
        ByteBuffer b = ByteBuffer.wrap(raw, 0, (int) need).order(order);
        switch (type) {
            case UINT8:
                for (int i = 0; i < n; i++) out[i] = b.get() & 0xFF;
                break;
            case INT8:
                for (int i = 0; i < n; i++) out[i] = b.get();
                break;
            case INT16: {
                var v = b.asShortBuffer();
                for (int i = 0; i < n; i++) out[i] = v.get();
                break;
            }
            case UINT16: {
                var v = b.asShortBuffer();
                for (int i = 0; i < n; i++) out[i] = v.get() & 0xFFFF;
                break;
            }
            case INT32: {
                var v = b.asIntBuffer();
                for (int i = 0; i < n; i++) out[i] = v.get();
                break;
            }
            case UINT32: {
                var v = b.asIntBuffer();
                for (int i = 0; i < n; i++) out[i] = v.get() & 0xFFFFFFFFL;
                break;
            }
            case INT64: {
                var v = b.asLongBuffer();
                for (int i = 0; i < n; i++) out[i] = v.get();
                break;
            }
            case UINT64: {
                var v = b.asLongBuffer();
                for (int i = 0; i < n; i++) out[i] = unsigned(v.get());
                break;
            }
            case FLOAT32: {
                var v = b.asFloatBuffer();
                for (int i = 0; i < n; i++) out[i] = v.get();
                break;
            }
            case FLOAT64: {
                var v = b.asDoubleBuffer();
                for (int i = 0; i < n; i++) out[i] = v.get();
                break;
            }
            default:
                throw new NiftiError("no decoder for " + type.label());
        }
        if (slope != 1.0 || inter != 0.0) {
            for (int i = 0; i < n; i++) out[i] = out[i] * slope + inter;
        }
        return NdArray.wrap(shape, out);
    }

    /**
     * Encode {@code arr} as raw voxel bytes, inverting {@code v * slope + inter}.
     * Integer targets are rounded half-to-even and clamped to the type range;
     * NaN becomes zero.
     */
    public static byte[] encode(NdArray arr, DataType type, ByteOrder order,
                                double slope, double inter) {
        if (!type.isReal()) {
            throw new NiftiError(type.label() + " voxel data is not supported yet");
        }
        int n = arr.size();
        long size = (long) n * type.itemSize();
        if (size > Integer.MAX_VALUE) {
            throw new NiftiError("voxel data would be " + size
                    + " bytes, over the 2 GiB nicloj can hold in one array");
        }
        ByteBuffer b = ByteBuffer.allocate((int) size).order(order);
        double[] src = arr.data();
        boolean unscale = slope != 1.0 || inter != 0.0;
        if (type.isInteger()) {
            double lo = type.minValue();
            double hi = type.maxValue();
            for (int i = 0; i < n; i++) {
                double v = unscale ? (src[i] - inter) / slope : src[i];
                v = Double.isNaN(v) ? 0.0 : Math.min(hi, Math.max(lo, Math.rint(v)));
                putInt(b, type, v);
            }
        } else {
            for (int i = 0; i < n; i++) {
                double v = unscale ? (src[i] - inter) / slope : src[i];
                if (type == DataType.FLOAT32) b.putFloat((float) v); else b.putDouble(v);
            }
        }
        return b.array();
    }

    /**
     * Pick {@code {slope, inter}} so that {@code arr} survives storage as
     * {@code type}. Returns {@code {1, 0}} when the values already fit, and
     * otherwise linearly maps the data range onto the type range.
     *
     * <p>Mirrors nibabel's {@code SlopeInterArrayWriter} for the common cases;
     * unlike nibabel it does not chase exact-representability of the scalers,
     * relying on clamping in {@link #encode} instead.
     */
    public static double[] autoScale(NdArray arr, DataType type, boolean float32Scalers) {
        if (!type.isInteger()) return new double[] {1.0, 0.0};
        double[] range = arr.finiteRange();
        double mn = range[0];
        double mx = range[1];
        if (mn > mx) return new double[] {1.0, 0.0};  // no finite data
        if (arr.hasNaN()) {
            mn = Math.min(mn, 0.0);
            mx = Math.max(mx, 0.0);
        }
        double lo = type.minValue();
        double hi = type.maxValue();
        if (arr.isIntegral() && mn >= lo && mx <= hi) return new double[] {1.0, 0.0};
        if (mn == mx) return new double[] {1.0, mn};

        double slope = (mx - mn) / (hi - lo);
        double inter;
        if (lo == 0 && Math.abs(mx) < Math.abs(mn)) {
            inter = mx;
            slope = -slope;
        } else {
            inter = mn - lo * slope;
        }
        if (float32Scalers) {
            slope = (float) slope;
            inter = (float) inter;
            // Widen slightly so float32 rounding cannot push values past the range.
            slope = (float) (slope * (1 + Math.ulp(1.0f)));
        }
        if (!Double.isFinite(slope) || !Double.isFinite(inter) || slope == 0) {
            throw new NiftiError("cannot scale data into " + type.label());
        }
        return new double[] {slope, inter};
    }

    private static void putInt(ByteBuffer b, DataType type, double v) {
        switch (type) {
            case UINT8: case INT8: b.put((byte) (long) v); break;
            case INT16: case UINT16: b.putShort((short) (long) v); break;
            case INT32: case UINT32: b.putInt((int) (long) v); break;
            case INT64: b.putLong((long) v); break;
            case UINT64: b.putLong(v >= TWO_POW_63 ? (long) (v - TWO_POW_64) : (long) v); break;
            default: throw new NiftiError("no integer encoder for " + type.label());
        }
    }

    private static double unsigned(long v) { return v >= 0 ? v : v + TWO_POW_64; }
}
