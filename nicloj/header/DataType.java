package nicloj.header;

/**
 * NIfTI voxel datatype codes.
 *
 * <p>Codes, bit widths and names follow the public-domain {@code nifti1.h}
 * reference header (https://nifti.nimh.nih.gov/nifti-1).
 */
public enum DataType {
    UINT8(2, 8, "uint8", Kind.UINT),
    INT16(4, 16, "int16", Kind.INT),
    INT32(8, 32, "int32", Kind.INT),
    FLOAT32(16, 32, "float32", Kind.FLOAT),
    COMPLEX64(32, 64, "complex64", Kind.COMPLEX),
    FLOAT64(64, 64, "float64", Kind.FLOAT),
    RGB24(128, 24, "rgb24", Kind.RGB),
    INT8(256, 8, "int8", Kind.INT),
    UINT16(512, 16, "uint16", Kind.UINT),
    UINT32(768, 32, "uint32", Kind.UINT),
    INT64(1024, 64, "int64", Kind.INT),
    UINT64(1280, 64, "uint64", Kind.UINT),
    COMPLEX128(1792, 128, "complex128", Kind.COMPLEX),
    RGBA32(2304, 32, "rgba32", Kind.RGB);

    public enum Kind { INT, UINT, FLOAT, COMPLEX, RGB }

    private final int code;
    private final int bitpix;
    private final String label;
    private final Kind kind;

    DataType(int code, int bitpix, String label, Kind kind) {
        this.code = code;
        this.bitpix = bitpix;
        this.label = label;
        this.kind = kind;
    }

    public int code() { return code; }
    public int bitpix() { return bitpix; }
    public int itemSize() { return bitpix / 8; }
    public String label() { return label; }
    public Kind kind() { return kind; }

    /** True for types {@link nicloj.array.Codec} can read and write. */
    public boolean isReal() { return kind == Kind.INT || kind == Kind.UINT || kind == Kind.FLOAT; }

    public boolean isInteger() { return kind == Kind.INT || kind == Kind.UINT; }

    /** Lowest representable value, for integer types only. */
    public double minValue() {
        switch (this) {
            case INT8: return -128;
            case INT16: return -32768;
            case INT32: return Integer.MIN_VALUE;
            case INT64: return Long.MIN_VALUE;
            case UINT8: case UINT16: case UINT32: case UINT64: return 0;
            default: throw new NiftiError(label + " is not an integer type");
        }
    }

    /** Highest representable value, for integer types only. */
    public double maxValue() {
        switch (this) {
            case INT8: return 127;
            case UINT8: return 255;
            case INT16: return 32767;
            case UINT16: return 65535;
            case INT32: return Integer.MAX_VALUE;
            case UINT32: return 4294967295.0;
            case INT64: return Long.MAX_VALUE;
            // 2^64 - 1 is not a double; use the largest double below 2^64 so
            // that clamping cannot wrap round to zero on encode.
            case UINT64: return Math.nextDown(0x1p64);
            default: throw new NiftiError(label + " is not an integer type");
        }
    }

    public static DataType fromCode(int code) {
        for (DataType t : values()) if (t.code == code) return t;
        throw new NiftiError("unsupported NIfTI datatype code " + code);
    }

    /** Look up by canonical label, numpy-style short code, or NIfTI macro name. */
    public static DataType fromLabel(String name) {
        String s = name.trim().toLowerCase();
        for (DataType t : values()) if (t.label.equals(s)) return t;
        switch (s) {
            case "u1": case "b": case "ubyte": case "mask": return UINT8;
            case "i1": case "byte": return INT8;
            case "i2": case "short": return INT16;
            case "u2": case "ushort": return UINT16;
            case "i4": case "int": return INT32;
            case "u4": case "uint": return UINT32;
            case "i8": case "long": return INT64;
            case "u8": case "ulong": return UINT64;
            case "f4": case "float": case "single": return FLOAT32;
            case "f8": case "double": return FLOAT64;
            case "c8": return COMPLEX64;
            case "c16": return COMPLEX128;
            case "rgb": return RGB24;
            case "rgba": return RGBA32;
            default: throw new NiftiError("unknown NIfTI datatype '" + name + "'");
        }
    }
}
