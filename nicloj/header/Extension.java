package nicloj.header;

import java.io.ByteArrayOutputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.ArrayList;
import java.util.List;

/**
 * A NIfTI header extension: an {@code ecode} tag plus an opaque payload.
 *
 * <p>On disk the fixed header is followed by a 4-byte extender. A non-zero
 * first extender byte means extensions follow, each one a 16-byte-aligned
 * block of {@code esize:int32, ecode:int32, payload}.
 */
public final class Extension {
    /** Extension payload alignment mandated by the NIfTI spec. */
    public static final int ALIGN = 16;

    public final int code;
    public final byte[] content;

    public Extension(int code, byte[] content) {
        this.code = code;
        this.content = content;
    }

    /** Total on-disk size, including the 8-byte block header and padding. */
    public int diskSize() {
        int need = 8 + content.length;
        return (need + ALIGN - 1) / ALIGN * ALIGN;
    }

    /** Read the extender byte and any extension blocks that follow it. */
    public static List<Extension> readAll(ByteBuffer buf, ByteOrder order, long limit) {
        List<Extension> out = new ArrayList<>();
        if (buf.remaining() < 4 || buf.get(buf.position()) == 0) return out;
        buf.position(buf.position() + 4);
        buf.order(order);
        while (buf.position() + 8 <= limit) {
            int esize = buf.getInt();
            int ecode = buf.getInt();
            if (esize < 8 || buf.position() + esize - 8 > limit) break;
            byte[] payload = new byte[esize - 8];
            buf.get(payload);
            out.add(new Extension(ecode, payload));
        }
        return out;
    }

    /** Serialise the extender byte plus every block, zero-padded to alignment. */
    public static byte[] writeAll(List<Extension> exts, ByteOrder order) {
        ByteArrayOutputStream out = new ByteArrayOutputStream();
        out.write(exts.isEmpty() ? 0 : 1);
        out.write(0); out.write(0); out.write(0);
        for (Extension e : exts) {
            int size = e.diskSize();
            ByteBuffer block = ByteBuffer.allocate(size).order(order);
            block.putInt(size).putInt(e.code).put(e.content);
            out.writeBytes(block.array());
        }
        return out.toByteArray();
    }

    /** Bytes needed for the extender plus every block. */
    public static int totalSize(List<Extension> exts) {
        int n = 4;
        for (Extension e : exts) n += e.diskSize();
        return n;
    }

    @Override
    public String toString() { return "Extension[code=" + code + ", " + content.length + " bytes]"; }
}
