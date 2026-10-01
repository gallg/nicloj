package nicloj.array;

import nicloj.header.DataType;

/**
 * Backing storage for an {@link NdArray}, one implementation per NIfTI voxel
 * type.
 *
 * <p>Everything above this layer speaks {@code double}; the store keeps the
 * values at their on-disk width, so an {@code int16} volume costs two bytes a
 * voxel rather than eight. {@link #getLong} gives the writer an exact path for
 * integer types, so {@code int64} values past 2^53 survive a read-write round
 * trip.
 */
public abstract class Store {
    public abstract int size();

    public abstract double get(int i);

    public abstract void set(int i, double v);

    /** A new, zeroed store of the same type, so array ops preserve the type. */
    public abstract Store alloc(int n);

    public abstract DataType type();

    /** Exact integer read; meaningful only when {@link #exact} is true. */
    public long getLong(int i) { return (long) get(i); }

    /** Exact integer write; meaningful only when {@link #exact} is true. */
    public void setLong(int i, long v) { set(i, v); }

    /** True when every element is a whole number held exactly, at any magnitude. */
    public boolean exact() { return false; }

    /** True when the type can hold NaN or an infinity. */
    public boolean canBeNonFinite() { return !exact(); }

    /** A zeroed store for {@code type}. */
    public static Store of(DataType type, int n) {
        switch (type) {
            case UINT8: return new UBytes(n);
            case INT8: return new Bytes(n);
            case INT16: return new Shorts(n);
            case UINT16: return new UShorts(n);
            case INT32: return new Ints(n);
            case UINT32: return new UInts(n);
            case INT64: return new Longs(n);
            case UINT64: return new ULongs(n);
            case FLOAT32: return new Floats(n);
            case FLOAT64: return new Doubles(n);
            default: throw new IllegalArgumentException("no store for " + type.label());
        }
    }

    private static final double TWO_POW_64 = 0x1p64;
    private static final double TWO_POW_63 = 0x1p63;

    public static final class Doubles extends Store {
        final double[] a;
        public Doubles(int n) { this.a = new double[n]; }
        public Doubles(double[] a) { this.a = a; }
        @Override public int size() { return a.length; }
        @Override public double get(int i) { return a[i]; }
        @Override public void set(int i, double v) { a[i] = v; }
        @Override public Store alloc(int n) { return new Doubles(n); }
        @Override public DataType type() { return DataType.FLOAT64; }
    }

    public static final class Floats extends Store {
        final float[] a;
        public Floats(int n) { this.a = new float[n]; }
        @Override public int size() { return a.length; }
        @Override public double get(int i) { return a[i]; }
        @Override public void set(int i, double v) { a[i] = (float) v; }
        @Override public Store alloc(int n) { return new Floats(n); }
        @Override public DataType type() { return DataType.FLOAT32; }
    }

    public static final class Bytes extends Store {
        final byte[] a;
        public Bytes(int n) { this.a = new byte[n]; }
        @Override public int size() { return a.length; }
        @Override public double get(int i) { return a[i]; }
        @Override public void set(int i, double v) { a[i] = (byte) (long) v; }
        @Override public long getLong(int i) { return a[i]; }
        @Override public void setLong(int i, long v) { a[i] = (byte) v; }
        @Override public Store alloc(int n) { return new Bytes(n); }
        @Override public DataType type() { return DataType.INT8; }
        @Override public boolean exact() { return true; }
    }

    public static final class UBytes extends Store {
        final byte[] a;
        public UBytes(int n) { this.a = new byte[n]; }
        @Override public int size() { return a.length; }
        @Override public double get(int i) { return a[i] & 0xFF; }
        @Override public void set(int i, double v) { a[i] = (byte) (long) v; }
        @Override public long getLong(int i) { return a[i] & 0xFF; }
        @Override public void setLong(int i, long v) { a[i] = (byte) v; }
        @Override public Store alloc(int n) { return new UBytes(n); }
        @Override public DataType type() { return DataType.UINT8; }
        @Override public boolean exact() { return true; }
    }

    public static final class Shorts extends Store {
        final short[] a;
        public Shorts(int n) { this.a = new short[n]; }
        @Override public int size() { return a.length; }
        @Override public double get(int i) { return a[i]; }
        @Override public void set(int i, double v) { a[i] = (short) (long) v; }
        @Override public long getLong(int i) { return a[i]; }
        @Override public void setLong(int i, long v) { a[i] = (short) v; }
        @Override public Store alloc(int n) { return new Shorts(n); }
        @Override public DataType type() { return DataType.INT16; }
        @Override public boolean exact() { return true; }
    }

    public static final class UShorts extends Store {
        final short[] a;
        public UShorts(int n) { this.a = new short[n]; }
        @Override public int size() { return a.length; }
        @Override public double get(int i) { return a[i] & 0xFFFF; }
        @Override public void set(int i, double v) { a[i] = (short) (long) v; }
        @Override public long getLong(int i) { return a[i] & 0xFFFF; }
        @Override public void setLong(int i, long v) { a[i] = (short) v; }
        @Override public Store alloc(int n) { return new UShorts(n); }
        @Override public DataType type() { return DataType.UINT16; }
        @Override public boolean exact() { return true; }
    }

    public static final class Ints extends Store {
        final int[] a;
        public Ints(int n) { this.a = new int[n]; }
        @Override public int size() { return a.length; }
        @Override public double get(int i) { return a[i]; }
        @Override public void set(int i, double v) { a[i] = (int) (long) v; }
        @Override public long getLong(int i) { return a[i]; }
        @Override public void setLong(int i, long v) { a[i] = (int) v; }
        @Override public Store alloc(int n) { return new Ints(n); }
        @Override public DataType type() { return DataType.INT32; }
        @Override public boolean exact() { return true; }
    }

    public static final class UInts extends Store {
        final int[] a;
        public UInts(int n) { this.a = new int[n]; }
        @Override public int size() { return a.length; }
        @Override public double get(int i) { return a[i] & 0xFFFFFFFFL; }
        @Override public void set(int i, double v) { a[i] = (int) (long) v; }
        @Override public long getLong(int i) { return a[i] & 0xFFFFFFFFL; }
        @Override public void setLong(int i, long v) { a[i] = (int) v; }
        @Override public Store alloc(int n) { return new UInts(n); }
        @Override public DataType type() { return DataType.UINT32; }
        @Override public boolean exact() { return true; }
    }

    public static final class Longs extends Store {
        final long[] a;
        public Longs(int n) { this.a = new long[n]; }
        @Override public int size() { return a.length; }
        @Override public double get(int i) { return a[i]; }
        @Override public void set(int i, double v) { a[i] = (long) v; }
        @Override public long getLong(int i) { return a[i]; }
        @Override public void setLong(int i, long v) { a[i] = v; }
        @Override public Store alloc(int n) { return new Longs(n); }
        @Override public DataType type() { return DataType.INT64; }
        @Override public boolean exact() { return true; }
    }

    /** {@code uint64} kept as raw bits; {@link #get} widens, {@link #getLong} does not. */
    public static final class ULongs extends Store {
        final long[] a;
        public ULongs(int n) { this.a = new long[n]; }
        @Override public int size() { return a.length; }
        @Override public double get(int i) { return a[i] >= 0 ? a[i] : a[i] + TWO_POW_64; }
        @Override public void set(int i, double v) {
            a[i] = v >= TWO_POW_63 ? (long) (v - TWO_POW_64) : (long) v;
        }
        @Override public long getLong(int i) { return a[i]; }
        @Override public void setLong(int i, long v) { a[i] = v; }
        @Override public Store alloc(int n) { return new ULongs(n); }
        @Override public DataType type() { return DataType.UINT64; }
        @Override public boolean exact() { return true; }
    }
}
