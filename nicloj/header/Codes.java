package nicloj.header;

import java.util.LinkedHashMap;
import java.util.Map;

/**
 * Symbolic names for the small enumerated fields of a NIfTI header
 * (transform spaces, measurement units, slice orders, statistical intents).
 *
 * <p>Values are from the public-domain {@code nifti1.h} reference header.
 */
public final class Codes {
    private Codes() {}

    public static final Map<Integer, String> XFORM = table(
            0, "unknown", 1, "scanner", 2, "aligned", 3, "talairach", 4, "mni", 5, "template");

    public static final Map<Integer, String> UNITS = table(
            0, "unknown", 1, "meter", 2, "mm", 3, "micron",
            8, "sec", 16, "msec", 24, "usec", 32, "hz", 40, "ppm", 48, "rads");

    public static final Map<Integer, String> SLICE = table(
            0, "unknown", 1, "sequential increasing", 2, "sequential decreasing",
            3, "alternating increasing", 4, "alternating decreasing",
            5, "alternating increasing 2", 6, "alternating decreasing 2");

    public static final Map<Integer, String> INTENT = table(
            0, "none", 2, "correlation", 3, "t test", 4, "f test", 5, "z score",
            6, "chi2", 7, "beta", 8, "binomial", 9, "gamma", 10, "poisson",
            11, "normal", 12, "f test noncentral", 13, "chi2 noncentral",
            14, "logistic", 15, "laplace", 16, "uniform", 17, "t test noncentral",
            18, "weibull", 19, "chi", 20, "invgauss", 21, "extval", 22, "p value",
            23, "log p value", 24, "log10 p value", 1001, "estimate", 1002, "label",
            1003, "neuroname", 1004, "matrix", 1005, "symmetric matrix",
            1006, "displacement vector", 1007, "vector", 1008, "pointset",
            1009, "triangle", 1010, "quaternion", 1011, "dimensionless",
            2001, "time series", 2002, "node index", 2003, "rgb vector",
            2004, "rgba vector", 2005, "shape");

    /** Human-readable name for a code, or {@code "<n>"} when unrecognised. */
    public static String label(Map<Integer, String> t, int code) {
        String s = t.get(code);
        return s != null ? s : "<" + code + ">";
    }

    /** Reverse lookup; accepts either a known name or a decimal string. */
    public static int code(Map<Integer, String> t, String label) {
        for (Map.Entry<Integer, String> e : t.entrySet()) {
            if (e.getValue().equalsIgnoreCase(label)) return e.getKey();
        }
        try {
            return Integer.parseInt(label.trim());
        } catch (NumberFormatException ex) {
            throw new NiftiError("unknown code name '" + label + "'");
        }
    }

    /** Split {@code xyzt_units} into its spatial and temporal halves. */
    public static int[] splitUnits(int xyztUnits) {
        int xyz = xyztUnits & 0x07;
        return new int[] {xyz, xyztUnits & 0x38};
    }

    private static Map<Integer, String> table(Object... kvs) {
        Map<Integer, String> m = new LinkedHashMap<>();
        for (int i = 0; i < kvs.length; i += 2) m.put((Integer) kvs[i], (String) kvs[i + 1]);
        return java.util.Collections.unmodifiableMap(m);
    }
}
