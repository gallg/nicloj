package nicloj.affine;

import nicloj.header.NiftiError;

/** Conversions between NIfTI header geometry fields and 4x4 affine matrices. */
public final class Affines {
    private Affines() {}

    /** Precision of the float32-stored qform quaternion, as used by nibabel. */
    public static final double QUATERNION_THRESHOLD = 3.5762787e-07;

    /** The qform fields describing one affine. */
    public static final class Qform {
        public final double quaternB;
        public final double quaternC;
        public final double quaternD;
        public final double[] qoffset;
        public final double[] zooms;
        public final double qfac;
        /** True when the source affine had non-orthogonal axes, which the qform drops. */
        public final boolean sheared;

        Qform(double[] quat, double[] qoffset, double[] zooms, double qfac, boolean sheared) {
            this.quaternB = quat[1];
            this.quaternC = quat[2];
            this.quaternD = quat[3];
            this.qoffset = qoffset;
            this.zooms = zooms;
            this.qfac = qfac;
            this.sheared = sheared;
        }
    }

    /** Build the affine encoded by a qform. {@code zooms} is {@code pixdim[1..3]}. */
    public static double[][] fromQform(double quaternB, double quaternC, double quaternD,
                                       double[] qoffset, double[] zooms, double qfac) {
        for (double z : zooms) {
            if (z < 0) throw new NiftiError("pixdim[1,2,3] should be positive");
        }
        if (qfac != 1.0 && qfac != -1.0) {
            throw new NiftiError("qfac (pixdim[0]) should be 1 or -1, got " + qfac);
        }
        double[] quat = Quaternions.fillPositive(
                new double[] {quaternB, quaternC, quaternD}, QUATERNION_THRESHOLD);
        double[][] r = Quaternions.toMatrix(quat);
        double[] vox = {zooms[0], zooms[1], zooms[2] * qfac};
        double[][] out = Mat.identity(4);
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) out[i][j] = r[i][j] * vox[j];
            out[i][3] = qoffset[i];
        }
        return out;
    }

    /**
     * Decompose a 4x4 affine into qform fields, dropping any shear by taking
     * the closest orthogonal rotation.
     */
    public static Qform toQform(double[][] affine) {
        if (affine.length != 4 || affine[0].length != 4) {
            throw new NiftiError("need a 4x4 affine");
        }
        double[][] rzs = Mat.block(affine, 3, 3);
        double[] zooms = Mat.columnNorms(rzs);
        for (double z : zooms) {
            if (z == 0.0) throw new NiftiError("cannot decompose a singular affine");
        }
        double[][] r = new double[3][3];
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) r[i][j] = rzs[i][j] / zooms[j];
        }
        double qfac = 1.0;
        if (Mat.det(r) <= 0) {
            qfac = -1.0;
            for (int i = 0; i < 3; i++) r[i][2] = -r[i][2];
        }
        double[][] pr = Svd.of(r).polar(0.0);
        boolean sheared = !Mat.close(pr, r, 1e-8);
        double[] quat = Quaternions.fromMatrix(pr);
        double[] qoffset = {affine[0][3], affine[1][3], affine[2][3]};
        return new Qform(quat, qoffset, zooms, qfac, sheared);
    }

    /** Assemble an affine from the three sform rows. */
    public static double[][] fromSrow(double[][] srow) {
        double[][] out = Mat.identity(4);
        for (int i = 0; i < 3; i++) System.arraycopy(srow[i], 0, out[i], 0, 4);
        return out;
    }

    /**
     * Fallback affine for headers with no qform or sform: a diagonal affine
     * centred on the image, optionally with the Analyze x-axis flip.
     */
    public static double[][] shapeZoomAffine(int[] shape, double[] zooms, boolean xFlip) {
        if (shape.length != zooms.length) {
            throw new NiftiError("shape and zooms must have the same length");
        }
        double[] s = {1, 1, 1};
        double[] z = {1, 1, 1};
        for (int i = 0; i < Math.min(3, shape.length); i++) {
            s[i] = shape[i];
            z[i] = zooms[i];
        }
        if (xFlip) z[0] = -z[0];
        double[][] aff = Mat.identity(4);
        for (int i = 0; i < 3; i++) {
            aff[i][i] = z[i];
            aff[i][3] = -(s[i] - 1) / 2.0 * z[i];
        }
        return aff;
    }

    /** Map one point through an affine. */
    public static double[] apply(double[][] affine, double[] point) {
        int n = affine.length - 1;
        double[] out = new double[n];
        for (int i = 0; i < n; i++) {
            double acc = affine[i][n];
            for (int j = 0; j < point.length; j++) acc += affine[i][j] * point[j];
            out[i] = acc;
        }
        return out;
    }

    /** Map a batch of points, one per row. */
    public static double[][] applyAll(double[][] affine, double[][] points) {
        double[][] out = new double[points.length][];
        for (int i = 0; i < points.length; i++) out[i] = apply(affine, points[i]);
        return out;
    }

    /** Voxel sizes implied by an affine: the column norms of its rotation block. */
    public static double[] voxelSizes(double[][] affine) {
        return Mat.columnNorms(Mat.block(affine, affine.length - 1, affine[0].length - 1));
    }

    /** Compose an affine from a linear part and a translation. */
    public static double[][] fromMatVec(double[][] linear, double[] translation) {
        int rows = linear.length;
        int cols = linear[0].length;
        double[][] out = new double[rows + 1][cols + 1];
        for (int i = 0; i < rows; i++) {
            System.arraycopy(linear[i], 0, out[i], 0, cols);
            out[i][cols] = translation[i];
        }
        out[rows][cols] = 1.0;
        return out;
    }
}
