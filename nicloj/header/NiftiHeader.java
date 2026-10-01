package nicloj.header;

import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * A mutable NIfTI header held in memory, covering both the NIfTI-1 (348 byte)
 * and NIfTI-2 (540 byte) on-disk layouts and the legacy Analyze 7.5 fields
 * that NIfTI-1 inherits.
 *
 * <p>Fields use the wider of the two representations ({@code long} /
 * {@code double}); narrowing happens on write. Field names, offsets and
 * defaults follow the public-domain {@code nifti1.h} and {@code nifti2.h}
 * reference headers (https://nifti.nimh.nih.gov/).
 *
 * <p>This class is a struct plus binary codec only: geometric interpretation
 * of qform/sform lives in {@code nicloj.affine} and {@code nicloj.core.header}.
 */
public final class NiftiHeader {
    public static final int NIFTI1_SIZE = 348;
    public static final int NIFTI2_SIZE = 540;

    /** Byte-swapped 348, i.e. a NIfTI-1 header written big-endian. */
    private static final int NIFTI1_SWAPPED = 0x5C010000;
    private static final int NIFTI2_SWAPPED = 0x1C020000;

    /** Trailing bytes of a NIfTI-2 magic: a CR/LF/EOF pair guarding against text mangling. */
    private static final byte[] NIFTI2_EOL = {0x0D, 0x0A, 0x1A, 0x0A};

    /** 1 for NIfTI-1, 2 for NIfTI-2. */
    public int version = 1;
    public ByteOrder order = ByteOrder.LITTLE_ENDIAN;
    /** True for n+1/n+2 magic, where header and voxels share one file. */
    public boolean singleFile = true;

    // --- Analyze 7.5 legacy fields, preserved for faithful round trips ---
    public String legacyDataType = "";
    public String dbName = "";
    public int extents;
    public int sessionError;
    public int regular;
    public int glmax;
    public int glmin;

    public int dimInfo;
    public long[] dim = {0, 1, 1, 1, 1, 1, 1, 1};
    public double intentP1;
    public double intentP2;
    public double intentP3;
    public int intentCode;
    public int datatype = DataType.FLOAT32.code();
    public int bitpix = DataType.FLOAT32.bitpix();
    public long sliceStart;
    public long sliceEnd;
    public double[] pixdim = {1, 1, 1, 1, 1, 1, 1, 1};
    public long voxOffset;
    public double sclSlope = 1.0;
    public double sclInter = 0.0;
    public int sliceCode;
    public int xyztUnits;
    public double calMax;
    public double calMin;
    public double sliceDuration;
    public double toffset;
    public String descrip = "";
    public String auxFile = "";
    public int qformCode;
    public int sformCode;
    public double quaternB;
    public double quaternC;
    public double quaternD;
    public double qoffsetX;
    public double qoffsetY;
    public double qoffsetZ;
    public double[][] srow = new double[3][4];
    public String intentName = "";
    public List<Extension> extensions = new ArrayList<>();

    public static NiftiHeader of(int version) {
        NiftiHeader h = new NiftiHeader();
        h.version = version;
        return h;
    }

    // ------------------------------------------------------------------ shape

    /** Voxel grid shape, i.e. dim[1..dim[0]]. */
    public int[] shape() {
        int n = (int) dim[0];
        if (n < 0 || n > 7) throw new NiftiError("dim[0] out of range: " + dim[0]);
        int[] s = new int[n];
        for (int i = 0; i < n; i++) {
            if (dim[i + 1] < 0 || dim[i + 1] > Integer.MAX_VALUE) {
                throw new NiftiError("dim[" + (i + 1) + "] out of range: " + dim[i + 1]);
            }
            s[i] = (int) dim[i + 1];
        }
        return s;
    }

    public void setShape(int[] shape) {
        if (shape.length > 7) throw new NiftiError("NIfTI supports at most 7 dimensions");
        Arrays.fill(dim, 1);
        dim[0] = shape.length;
        for (int i = 0; i < shape.length; i++) dim[i + 1] = shape[i];
    }

    /** Voxel sizes along each used dimension, i.e. pixdim[1..dim[0]]. */
    public double[] zooms() {
        int n = (int) dim[0];
        if (n == 0) return new double[] {1.0};
        return Arrays.copyOfRange(pixdim, 1, n + 1);
    }

    public void setZooms(double[] zooms) {
        int n = (int) dim[0];
        if (zooms.length != n) {
            throw new NiftiError("expected " + n + " zooms for a " + n + "-d image");
        }
        for (double z : zooms) if (z < 0) throw new NiftiError("zooms must be positive");
        System.arraycopy(zooms, 0, pixdim, 1, n);
    }

    public DataType dataType() { return DataType.fromCode(datatype); }

    public void setDataType(DataType t) {
        datatype = t.code();
        bitpix = t.bitpix();
    }

    public long voxelCount() {
        long n = 1;
        for (int s : shape()) n *= s;
        return n;
    }

    /** Bytes of voxel data implied by shape and datatype. */
    public long dataBytes() { return voxelCount() * dataType().itemSize(); }

    /** pixdim[0], the qform handedness flag; 0 is treated as 1. */
    public double qfac() { return pixdim[0] < 0 ? -1.0 : 1.0; }

    /** Size of the fixed part of the header for this version. */
    public int structSize() { return version == 2 ? NIFTI2_SIZE : NIFTI1_SIZE; }

    /** Lowest legal vox_offset for a single-file image. */
    public long minVoxOffset() { return structSize() + Extension.totalSize(extensions); }

    /** The 3- or 4-character magic tag, without the NIfTI-2 end-of-line guard. */
    public String magic() { return (singleFile ? "n+" : "ni") + version; }

    // ------------------------------------------------------------------- read

    /**
     * Peek at sizeof_hdr to tell which layout and byte order a buffer holds.
     * Returns {version, swapped} where swapped is 1 for big-endian.
     */
    public static int[] probe(byte[] buf) {
        if (buf.length < 4) throw new NiftiError("too short to be a NIfTI header");
        int le = ByteBuffer.wrap(buf).order(ByteOrder.LITTLE_ENDIAN).getInt(0);
        if (le == NIFTI1_SIZE) return new int[] {1, 0};
        if (le == NIFTI2_SIZE) return new int[] {2, 0};
        if (le == NIFTI1_SWAPPED) return new int[] {1, 1};
        if (le == NIFTI2_SWAPPED) return new int[] {2, 1};
        throw new NiftiError("not a NIfTI file: sizeof_hdr is " + le + ", expected 348 or 540");
    }

    /**
     * Parse a header, and any extensions present, from the start of {@code buf}.
     * The buffer may extend past the header.
     */
    public static NiftiHeader read(byte[] buf) {
        int[] p = probe(buf);
        NiftiHeader h = new NiftiHeader();
        h.version = p[0];
        h.order = p[1] == 0 ? ByteOrder.LITTLE_ENDIAN : ByteOrder.BIG_ENDIAN;
        if (buf.length < h.structSize()) {
            throw new NiftiError("truncated NIfTI-" + h.version + " header: " + buf.length + " bytes");
        }
        ByteBuffer b = ByteBuffer.wrap(buf).order(h.order);
        String magic = h.version == 2 ? h.read2(b) : h.read1(b);
        if (!magic.equals("n+" + h.version) && !magic.equals("ni" + h.version)) {
            throw new NiftiError("bad NIfTI magic '" + magic + "'");
        }
        h.singleFile = magic.charAt(1) == '+';
        // nibabel's load-time fix-ups: qfac other than +-1 becomes 1, and
        // spatial pixdims become positive, with 0 read as 1.
        if (h.pixdim[0] != 1 && h.pixdim[0] != -1) h.pixdim[0] = 1;
        for (int i = 1; i <= 3; i++) h.pixdim[i] = h.pixdim[i] == 0 ? 1 : Math.abs(h.pixdim[i]);
        long limit = h.singleFile && h.voxOffset > 0
                ? Math.min(h.voxOffset, buf.length) : buf.length;
        b.position(h.structSize());
        h.extensions = Extension.readAll(b, limit);
        return h;
    }

    private String read1(ByteBuffer b) {
        b.position(4);
        legacyDataType = str(b, 10);
        dbName = str(b, 18);
        extents = b.getInt();
        sessionError = b.getShort();
        regular = b.get() & 0xFF;
        dimInfo = b.get() & 0xFF;
        for (int i = 0; i < 8; i++) dim[i] = b.getShort();
        intentP1 = b.getFloat();
        intentP2 = b.getFloat();
        intentP3 = b.getFloat();
        intentCode = b.getShort();
        datatype = b.getShort();
        bitpix = b.getShort();
        sliceStart = b.getShort();
        for (int i = 0; i < 8; i++) pixdim[i] = b.getFloat();
        voxOffset = (long) b.getFloat();
        sclSlope = b.getFloat();
        sclInter = b.getFloat();
        sliceEnd = b.getShort();
        sliceCode = b.get() & 0xFF;
        xyztUnits = b.get() & 0xFF;
        calMax = b.getFloat();
        calMin = b.getFloat();
        sliceDuration = b.getFloat();
        toffset = b.getFloat();
        glmax = b.getInt();
        glmin = b.getInt();
        descrip = str(b, 80);
        auxFile = str(b, 24);
        qformCode = b.getShort();
        sformCode = b.getShort();
        quaternB = b.getFloat();
        quaternC = b.getFloat();
        quaternD = b.getFloat();
        qoffsetX = b.getFloat();
        qoffsetY = b.getFloat();
        qoffsetZ = b.getFloat();
        for (int r = 0; r < 3; r++) for (int c = 0; c < 4; c++) srow[r][c] = b.getFloat();
        intentName = str(b, 16);
        return str(b, 4);
    }

    private String read2(ByteBuffer b) {
        b.position(4);
        String magic = str(b, 8);
        datatype = b.getShort();
        bitpix = b.getShort();
        for (int i = 0; i < 8; i++) dim[i] = b.getLong();
        intentP1 = b.getDouble();
        intentP2 = b.getDouble();
        intentP3 = b.getDouble();
        for (int i = 0; i < 8; i++) pixdim[i] = b.getDouble();
        voxOffset = b.getLong();
        sclSlope = b.getDouble();
        sclInter = b.getDouble();
        calMax = b.getDouble();
        calMin = b.getDouble();
        sliceDuration = b.getDouble();
        toffset = b.getDouble();
        sliceStart = b.getLong();
        sliceEnd = b.getLong();
        descrip = str(b, 80);
        auxFile = str(b, 24);
        qformCode = b.getInt();
        sformCode = b.getInt();
        quaternB = b.getDouble();
        quaternC = b.getDouble();
        quaternD = b.getDouble();
        qoffsetX = b.getDouble();
        qoffsetY = b.getDouble();
        qoffsetZ = b.getDouble();
        for (int r = 0; r < 3; r++) for (int c = 0; c < 4; c++) srow[r][c] = b.getDouble();
        sliceCode = b.getInt();
        xyztUnits = b.getInt();
        intentCode = b.getInt();
        intentName = str(b, 16);
        dimInfo = b.get() & 0xFF;
        return magic;
    }

    // ------------------------------------------------------------------ write

    /** Serialise the fixed header only. */
    public byte[] toBytes() {
        ByteBuffer b = ByteBuffer.allocate(structSize()).order(order);
        if (version == 2) write2(b); else write1(b);
        return b.array();
    }

    /** Serialise the fixed header followed by the extender and extension blocks. */
    public byte[] toBytesWithExtensions() {
        byte[] head = toBytes();
        byte[] ext = Extension.writeAll(extensions, order);
        byte[] out = new byte[head.length + ext.length];
        System.arraycopy(head, 0, out, 0, head.length);
        System.arraycopy(ext, 0, out, head.length, ext.length);
        return out;
    }

    private void write1(ByteBuffer b) {
        b.putInt(NIFTI1_SIZE);
        putStr(b, legacyDataType, 10);
        putStr(b, dbName, 18);
        b.putInt(extents);
        b.putShort((short) sessionError);
        b.put((byte) regular);
        b.put((byte) dimInfo);
        for (int i = 0; i < 8; i++) b.putShort(toShort(dim[i], "dim[" + i + "]"));
        b.putFloat((float) intentP1);
        b.putFloat((float) intentP2);
        b.putFloat((float) intentP3);
        b.putShort((short) intentCode);
        b.putShort((short) datatype);
        b.putShort((short) bitpix);
        b.putShort(toShort(sliceStart, "slice_start"));
        for (int i = 0; i < 8; i++) b.putFloat((float) pixdim[i]);
        b.putFloat((float) voxOffset);
        b.putFloat((float) sclSlope);
        b.putFloat((float) sclInter);
        b.putShort(toShort(sliceEnd, "slice_end"));
        b.put((byte) sliceCode);
        b.put((byte) xyztUnits);
        b.putFloat((float) calMax);
        b.putFloat((float) calMin);
        b.putFloat((float) sliceDuration);
        b.putFloat((float) toffset);
        b.putInt(glmax);
        b.putInt(glmin);
        putStr(b, descrip, 80);
        putStr(b, auxFile, 24);
        b.putShort((short) qformCode);
        b.putShort((short) sformCode);
        b.putFloat((float) quaternB);
        b.putFloat((float) quaternC);
        b.putFloat((float) quaternD);
        b.putFloat((float) qoffsetX);
        b.putFloat((float) qoffsetY);
        b.putFloat((float) qoffsetZ);
        for (int r = 0; r < 3; r++) for (int c = 0; c < 4; c++) b.putFloat((float) srow[r][c]);
        putStr(b, intentName, 16);
        putStr(b, magic(), 4);
    }

    private void write2(ByteBuffer b) {
        b.putInt(NIFTI2_SIZE);
        putStr(b, magic(), 4);
        b.put(NIFTI2_EOL);
        b.putShort((short) datatype);
        b.putShort((short) bitpix);
        for (int i = 0; i < 8; i++) b.putLong(dim[i]);
        b.putDouble(intentP1);
        b.putDouble(intentP2);
        b.putDouble(intentP3);
        for (int i = 0; i < 8; i++) b.putDouble(pixdim[i]);
        b.putLong(voxOffset);
        b.putDouble(sclSlope);
        b.putDouble(sclInter);
        b.putDouble(calMax);
        b.putDouble(calMin);
        b.putDouble(sliceDuration);
        b.putDouble(toffset);
        b.putLong(sliceStart);
        b.putLong(sliceEnd);
        putStr(b, descrip, 80);
        putStr(b, auxFile, 24);
        b.putInt(qformCode);
        b.putInt(sformCode);
        b.putDouble(quaternB);
        b.putDouble(quaternC);
        b.putDouble(quaternD);
        b.putDouble(qoffsetX);
        b.putDouble(qoffsetY);
        b.putDouble(qoffsetZ);
        for (int r = 0; r < 3; r++) for (int c = 0; c < 4; c++) b.putDouble(srow[r][c]);
        b.putInt(sliceCode);
        b.putInt(xyztUnits);
        b.putInt(intentCode);
        putStr(b, intentName, 16);
        b.put((byte) dimInfo);
        putStr(b, "", 15);
    }

    // ------------------------------------------------------------------- misc

    public NiftiHeader copy() {
        NiftiHeader h = new NiftiHeader();
        h.version = version;
        h.order = order;
        h.singleFile = singleFile;
        h.legacyDataType = legacyDataType;
        h.dbName = dbName;
        h.extents = extents;
        h.sessionError = sessionError;
        h.regular = regular;
        h.glmax = glmax;
        h.glmin = glmin;
        h.dimInfo = dimInfo;
        h.dim = dim.clone();
        h.intentP1 = intentP1;
        h.intentP2 = intentP2;
        h.intentP3 = intentP3;
        h.intentCode = intentCode;
        h.datatype = datatype;
        h.bitpix = bitpix;
        h.sliceStart = sliceStart;
        h.sliceEnd = sliceEnd;
        h.pixdim = pixdim.clone();
        h.voxOffset = voxOffset;
        h.sclSlope = sclSlope;
        h.sclInter = sclInter;
        h.sliceCode = sliceCode;
        h.xyztUnits = xyztUnits;
        h.calMax = calMax;
        h.calMin = calMin;
        h.sliceDuration = sliceDuration;
        h.toffset = toffset;
        h.descrip = descrip;
        h.auxFile = auxFile;
        h.qformCode = qformCode;
        h.sformCode = sformCode;
        h.quaternB = quaternB;
        h.quaternC = quaternC;
        h.quaternD = quaternD;
        h.qoffsetX = qoffsetX;
        h.qoffsetY = qoffsetY;
        h.qoffsetZ = qoffsetZ;
        h.srow = new double[][] {srow[0].clone(), srow[1].clone(), srow[2].clone()};
        h.intentName = intentName;
        h.extensions = new ArrayList<>(extensions);
        return h;
    }

    @Override
    public String toString() {
        return "NiftiHeader[v" + version + " " + Arrays.toString(shape()) + " "
                + dataType().label() + " " + (order == ByteOrder.LITTLE_ENDIAN ? "LE" : "BE") + "]";
    }

    private static short toShort(long v, String field) {
        if (v < Short.MIN_VALUE || v > Short.MAX_VALUE) {
            throw new NiftiError(field + " = " + v + " does not fit NIfTI-1; use NIfTI-2");
        }
        return (short) v;
    }

    /** Read {@code n} bytes as a NUL-terminated Latin-1 string. */
    private static String str(ByteBuffer b, int n) {
        byte[] raw = new byte[n];
        b.get(raw);
        int len = 0;
        while (len < n && raw[len] != 0) len++;
        return new String(raw, 0, len, StandardCharsets.ISO_8859_1);
    }

    /** Write {@code s} into an {@code n}-byte NUL-padded field, truncating if needed. */
    private static void putStr(ByteBuffer b, String s, int n) {
        byte[] raw = s.getBytes(StandardCharsets.ISO_8859_1);
        int len = Math.min(raw.length, n);
        b.put(raw, 0, len);
        for (int i = len; i < n; i++) b.put((byte) 0);
    }
}
