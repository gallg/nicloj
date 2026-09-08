package nicloj.affine;

/**
 * Thin singular value decomposition for the small matrices found in image
 * affines, via one-sided Jacobi rotations.
 *
 * <p>Yields {@code A == U * diag(s) * transpose(V)} with {@code U} of shape
 * {@code (m, k)}, {@code V} of shape {@code (n, k)} and {@code k = min(m, n)},
 * singular values in descending order. This is the "economy" decomposition,
 * equivalent to numpy's {@code svd(a, full_matrices=False)} up to column signs.
 */
public final class Svd {
    private static final int MAX_SWEEPS = 60;
    private static final double EPS = 2.220446049250313e-16;

    public final double[][] u;
    public final double[] s;
    public final double[][] v;

    private Svd(double[][] u, double[] s, double[][] v) {
        this.u = u;
        this.s = s;
        this.v = v;
    }

    public static Svd of(double[][] a) {
        if (a.length >= a[0].length) return jacobi(a);
        Svd t = jacobi(Mat.transpose(a));
        return new Svd(t.v, t.s, t.u);
    }

    /** One-sided Jacobi, requires {@code m >= n}. */
    private static Svd jacobi(double[][] a) {
        int m = a.length;
        int n = a[0].length;
        double[][] b = Mat.copy(a);
        double[][] v = Mat.identity(n);
        for (int sweep = 0; sweep < MAX_SWEEPS; sweep++) {
            boolean rotated = false;
            for (int p = 0; p < n - 1; p++) {
                for (int q = p + 1; q < n; q++) {
                    double alpha = 0, beta = 0, gamma = 0;
                    for (int i = 0; i < m; i++) {
                        alpha += b[i][p] * b[i][p];
                        beta += b[i][q] * b[i][q];
                        gamma += b[i][p] * b[i][q];
                    }
                    if (gamma == 0.0 || Math.abs(gamma) <= EPS * Math.sqrt(alpha * beta)) continue;
                    rotated = true;
                    double zeta = (beta - alpha) / (2.0 * gamma);
                    double t = Math.signum(zeta) / (Math.abs(zeta) + Math.sqrt(1.0 + zeta * zeta));
                    if (zeta == 0.0) t = 1.0;
                    double c = 1.0 / Math.sqrt(1.0 + t * t);
                    double sn = c * t;
                    rotate(b, m, p, q, c, sn);
                    rotate(v, n, p, q, c, sn);
                }
            }
            if (!rotated) break;
        }
        double[] s = Mat.columnNorms(b);
        double[][] u = new double[m][n];
        for (int j = 0; j < n; j++) {
            if (s[j] == 0.0) continue;
            for (int i = 0; i < m; i++) u[i][j] = b[i][j] / s[j];
        }
        sortDescending(s, u, v);
        return new Svd(u, s, v);
    }

    private static void rotate(double[][] mtx, int rows, int p, int q, double c, double s) {
        for (int i = 0; i < rows; i++) {
            double xp = mtx[i][p];
            double xq = mtx[i][q];
            mtx[i][p] = c * xp - s * xq;
            mtx[i][q] = s * xp + c * xq;
        }
    }

    private static void sortDescending(double[] s, double[][] u, double[][] v) {
        for (int i = 0; i < s.length - 1; i++) {
            int best = i;
            for (int j = i + 1; j < s.length; j++) if (s[j] > s[best]) best = j;
            if (best == i) continue;
            double t = s[i]; s[i] = s[best]; s[best] = t;
            swapColumns(u, i, best);
            swapColumns(v, i, best);
        }
    }

    private static void swapColumns(double[][] mtx, int i, int j) {
        for (double[] row : mtx) {
            double t = row[i]; row[i] = row[j]; row[j] = t;
        }
    }

    /**
     * Orthogonal factor of the polar decomposition, i.e. the orthogonal matrix
     * closest to {@code a}. Singular values at or below {@code tol} are dropped,
     * so a rank-deficient input yields a rank-deficient result.
     */
    public double[][] polar(double tol) {
        int m = u.length;
        int n = v.length;
        double[][] r = new double[m][n];
        for (int k = 0; k < s.length; k++) {
            if (s[k] <= tol) continue;
            for (int i = 0; i < m; i++) {
                for (int j = 0; j < n; j++) r[i][j] += u[i][k] * v[j][k];
            }
        }
        return r;
    }

    /** Default rank tolerance, matching numpy's convention. */
    public double defaultTol(int m, int n) {
        return s.length == 0 ? 0 : s[0] * Math.max(m, n) * EPS;
    }
}
