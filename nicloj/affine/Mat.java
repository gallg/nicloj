package nicloj.affine;

/** Small dense matrix algebra on row-major {@code double[][]}. */
public final class Mat {
    private Mat() {}

    public static double[][] identity(int n) {
        double[][] m = new double[n][n];
        for (int i = 0; i < n; i++) m[i][i] = 1.0;
        return m;
    }

    public static double[][] diag(double[] d) {
        double[][] m = new double[d.length][d.length];
        for (int i = 0; i < d.length; i++) m[i][i] = d[i];
        return m;
    }

    public static double[][] copy(double[][] a) {
        double[][] m = new double[a.length][];
        for (int i = 0; i < a.length; i++) m[i] = a[i].clone();
        return m;
    }

    public static double[][] transpose(double[][] a) {
        double[][] t = new double[a[0].length][a.length];
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a[0].length; j++) t[j][i] = a[i][j];
        }
        return t;
    }

    public static double[][] mul(double[][] a, double[][] b) {
        int n = a.length;
        int k = b.length;
        int m = b[0].length;
        if (a[0].length != k) {
            throw new IllegalArgumentException("shape mismatch: " + a[0].length + " vs " + k);
        }
        double[][] c = new double[n][m];
        for (int i = 0; i < n; i++) {
            for (int p = 0; p < k; p++) {
                double aip = a[i][p];
                if (aip == 0.0) continue;
                for (int j = 0; j < m; j++) c[i][j] += aip * b[p][j];
            }
        }
        return c;
    }


    /** Top-left {@code rows x cols} block. */
    public static double[][] block(double[][] a, int rows, int cols) {
        double[][] m = new double[rows][cols];
        for (int i = 0; i < rows; i++) System.arraycopy(a[i], 0, m[i], 0, cols);
        return m;
    }

    public static double[] column(double[][] a, int j) {
        double[] c = new double[a.length];
        for (int i = 0; i < a.length; i++) c[i] = a[i][j];
        return c;
    }

    /** Euclidean norms of each column. */
    public static double[] columnNorms(double[][] a) {
        double[] out = new double[a[0].length];
        for (double[] row : a) {
            for (int j = 0; j < out.length; j++) out[j] += row[j] * row[j];
        }
        for (int j = 0; j < out.length; j++) out[j] = Math.sqrt(out[j]);
        return out;
    }

    /** Gauss-Jordan inverse with partial pivoting. */
    public static double[][] inverse(double[][] a) {
        int n = a.length;
        if (a[0].length != n) throw new IllegalArgumentException("matrix must be square");
        double[][] w = copy(a);
        double[][] inv = identity(n);
        for (int col = 0; col < n; col++) {
            int pivot = col;
            for (int r = col + 1; r < n; r++) {
                if (Math.abs(w[r][col]) > Math.abs(w[pivot][col])) pivot = r;
            }
            if (w[pivot][col] == 0.0) throw new ArithmeticException("singular matrix");
            double[] t = w[col]; w[col] = w[pivot]; w[pivot] = t;
            t = inv[col]; inv[col] = inv[pivot]; inv[pivot] = t;
            double d = w[col][col];
            for (int j = 0; j < n; j++) { w[col][j] /= d; inv[col][j] /= d; }
            for (int r = 0; r < n; r++) {
                if (r == col) continue;
                double f = w[r][col];
                if (f == 0.0) continue;
                for (int j = 0; j < n; j++) {
                    w[r][j] -= f * w[col][j];
                    inv[r][j] -= f * inv[col][j];
                }
            }
        }
        return inv;
    }

    /** Determinant via LU decomposition with partial pivoting. */
    public static double det(double[][] a) {
        int n = a.length;
        double[][] w = copy(a);
        double d = 1.0;
        for (int col = 0; col < n; col++) {
            int pivot = col;
            for (int r = col + 1; r < n; r++) {
                if (Math.abs(w[r][col]) > Math.abs(w[pivot][col])) pivot = r;
            }
            if (w[pivot][col] == 0.0) return 0.0;
            if (pivot != col) {
                double[] t = w[col]; w[col] = w[pivot]; w[pivot] = t;
                d = -d;
            }
            d *= w[col][col];
            for (int r = col + 1; r < n; r++) {
                double f = w[r][col] / w[col][col];
                for (int j = col; j < n; j++) w[r][j] -= f * w[col][j];
            }
        }
        return d;
    }

    public static boolean close(double[][] a, double[][] b, double atol) {
        if (a.length != b.length || a[0].length != b[0].length) return false;
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a[0].length; j++) {
                if (!(Math.abs(a[i][j] - b[i][j]) <= atol)) return false;
            }
        }
        return true;
    }

    public static String format(double[][] a) {
        StringBuilder sb = new StringBuilder();
        for (double[] row : a) {
            for (double v : row) sb.append(String.format("%12.6g", v));
            sb.append('\n');
        }
        return sb.toString();
    }
}
