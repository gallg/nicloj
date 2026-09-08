package nicloj.affine;

import nicloj.header.NiftiError;

/**
 * Unit quaternion helpers for the NIfTI qform, in {@code (w, x, y, z)} order.
 *
 * <p>Algorithms follow the public-domain {@code nifti1.h} reference
 * implementation and nibabel's {@code nibabel.quaternions} (MIT).
 */
public final class Quaternions {
    private Quaternions() {}

    private static final double EPS = 2.220446049250313e-16;

    /**
     * Recover the non-negative real part of a unit quaternion from its vector
     * part, as the NIfTI qform stores only {@code b, c, d}.
     */
    public static double[] fillPositive(double[] xyz, double w2Thresh) {
        if (xyz.length != 3) throw new IllegalArgumentException("expected 3 components");
        double w2 = 1.0 - (xyz[0] * xyz[0] + xyz[1] * xyz[1] + xyz[2] * xyz[2]);
        double w;
        if (Math.abs(w2) < Math.abs(w2Thresh)) {
            w = 0.0;
        } else if (w2 < 0) {
            throw new NiftiError("quaternion b,c,d do not describe a rotation (w^2 = " + w2 + ")");
        } else {
            w = Math.sqrt(w2);
        }
        return new double[] {w, xyz[0], xyz[1], xyz[2]};
    }

    /** Rotation matrix for a quaternion; non-unit input is normalised. */
    public static double[][] toMatrix(double[] q) {
        double w = q[0], x = q[1], y = q[2], z = q[3];
        double nq = w * w + x * x + y * y + z * z;
        if (nq < EPS) return Mat.identity(3);
        double s = 2.0 / nq;
        double bigX = x * s, bigY = y * s, bigZ = z * s;
        double wX = w * bigX, wY = w * bigY, wZ = w * bigZ;
        double xX = x * bigX, xY = x * bigY, xZ = x * bigZ;
        double yY = y * bigY, yZ = y * bigZ, zZ = z * bigZ;
        return new double[][] {
                {1.0 - (yY + zZ), xY - wZ, xZ + wY},
                {xY + wZ, 1.0 - (xX + zZ), yZ - wX},
                {xZ - wY, yZ + wX, 1.0 - (xX + yY)}};
    }

    /**
     * Quaternion for a rotation matrix, by Shepperd's method, normalised to a
     * non-negative real part. The input must be orthogonal with determinant
     * {@code +1}; see {@link Svd#polar} to project a general matrix first.
     */
    public static double[] fromMatrix(double[][] m) {
        double trace = m[0][0] + m[1][1] + m[2][2];
        double w, x, y, z;
        if (trace > 0) {
            double s = Math.sqrt(trace + 1.0) * 2.0;
            w = 0.25 * s;
            x = (m[2][1] - m[1][2]) / s;
            y = (m[0][2] - m[2][0]) / s;
            z = (m[1][0] - m[0][1]) / s;
        } else if (m[0][0] > m[1][1] && m[0][0] > m[2][2]) {
            double s = Math.sqrt(1.0 + m[0][0] - m[1][1] - m[2][2]) * 2.0;
            w = (m[2][1] - m[1][2]) / s;
            x = 0.25 * s;
            y = (m[0][1] + m[1][0]) / s;
            z = (m[0][2] + m[2][0]) / s;
        } else if (m[1][1] > m[2][2]) {
            double s = Math.sqrt(1.0 + m[1][1] - m[0][0] - m[2][2]) * 2.0;
            w = (m[0][2] - m[2][0]) / s;
            x = (m[0][1] + m[1][0]) / s;
            y = 0.25 * s;
            z = (m[1][2] + m[2][1]) / s;
        } else {
            double s = Math.sqrt(1.0 + m[2][2] - m[0][0] - m[1][1]) * 2.0;
            w = (m[1][0] - m[0][1]) / s;
            x = (m[0][2] + m[2][0]) / s;
            y = (m[1][2] + m[2][1]) / s;
            z = 0.25 * s;
        }
        return w < 0 ? new double[] {-w, -x, -y, -z} : new double[] {w, x, y, z};
    }
}
