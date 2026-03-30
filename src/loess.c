/*
 * loess.c — Standalone 1D LOESS (local polynomial regression).
 *
 * Two modes:
 *   loess_fit():        Direct evaluation at each point (more accurate)
 *   loess_fit_interp(): Evaluate at kdtree vertices + Hermite interpolation
 *                       (approximates R's default loess behavior)
 *
 * Copyright: Cleveland's algorithm is public domain (AT&T 1989, 1992).
 * This implementation is original C code following the same algorithm.
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>

static int cmp_double_idx(const void *a, const void *b) {
    double da = ((const double *)a)[0], db = ((const double *)b)[0];
    return (da > db) - (da < db);
}

/*
 * Core: evaluate 1D weighted local polynomial at a single point.
 * Returns value in *val and derivative in *deriv.
 */
static void loess_eval_one(
    const double *x, const double *y, const double *w, int n,
    double xi, int h, int degree,
    double *val, double *deriv)
{
    int d = degree + 1;  /* number of coefficients */
    double *dist_idx = (double *)malloc(n * 2 * sizeof(double));

    for (int j = 0; j < n; j++) {
        dist_idx[j * 2] = fabs(x[j] - xi);
        dist_idx[j * 2 + 1] = (double)j;
    }
    qsort(dist_idx, n, 2 * sizeof(double), cmp_double_idx);

    double max_dist = dist_idx[(h - 1) * 2];
    if (max_dist == 0.0) max_dist = 1.0;

    /* Build normal equations */
    double ata[9] = {0};  /* max 3x3 */
    double aty[3] = {0};

    for (int k = 0; k < h; k++) {
        int j = (int)dist_idx[k * 2 + 1];
        double u = dist_idx[k * 2] / max_dist;
        double kern = 0.0;
        if (u < 1.0) {
            double u3 = u * u * u;
            double t = 1.0 - u3;
            kern = t * t * t;
        }
        double wt = kern * w[j];
        if (wt == 0.0) continue;

        double xc = x[j] - xi;
        double v[3];
        v[0] = 1.0;
        for (int p = 1; p < d; p++) v[p] = v[p - 1] * xc;

        for (int r = 0; r < d; r++) {
            for (int c = 0; c < d; c++)
                ata[r * d + c] += wt * v[r] * v[c];
            aty[r] += wt * v[r] * y[j];
        }
    }

    /* Solve via Gaussian elimination with partial pivoting */
    double A[9], b[3];
    memcpy(A, ata, d * d * sizeof(double));
    memcpy(b, aty, d * sizeof(double));

    int ok = 1;
    for (int col = 0; col < d; col++) {
        int piv = col;
        for (int row = col + 1; row < d; row++)
            if (fabs(A[row * d + col]) > fabs(A[piv * d + col])) piv = row;
        if (piv != col) {
            for (int c = 0; c < d; c++) {
                double tmp = A[col*d+c]; A[col*d+c] = A[piv*d+c]; A[piv*d+c] = tmp;
            }
            double tmp = b[col]; b[col] = b[piv]; b[piv] = tmp;
        }
        if (fabs(A[col * d + col]) < 1e-15) { ok = 0; break; }
        for (int row = col + 1; row < d; row++) {
            double factor = A[row * d + col] / A[col * d + col];
            for (int c = col; c < d; c++)
                A[row * d + c] -= factor * A[col * d + c];
            b[row] -= factor * b[col];
        }
    }

    if (ok) {
        for (int row = d - 1; row >= 0; row--) {
            for (int c = row + 1; c < d; c++)
                b[row] -= A[row * d + c] * b[c];
            b[row] /= A[row * d + row];
        }
        *val = b[0];
        *deriv = (d > 1) ? b[1] : 0.0;
    } else {
        *val = 0.0;
        *deriv = 0.0;
    }

    free(dist_idx);
}

/*
 * Direct LOESS: evaluate at each output point independently.
 * Most accurate, matches R's loess(control=loess.control(surface="direct")).
 */
void loess_fit(const double *x, const double *y, const double *w,
               double *out, int n, double span, int degree)
{
    int h = (int)floor(span * n);
    if (h < degree + 1) h = degree + 1;
    if (h > n) h = n;

    double dummy_deriv;
    for (int i = 0; i < n; i++)
        loess_eval_one(x, y, w, n, x[i], h, degree, &out[i], &dummy_deriv);
}

/*
 * Interpolated LOESS: evaluate at nv evenly-spaced vertices, then use
 * Hermite cubic interpolation between them. Approximates R's default
 * loess(control=loess.control(surface="interpolate")) behavior.
 *
 * nv: number of vertices (R uses ~9 for n~41)
 * If nv <= 0, auto-select based on n.
 */
void loess_fit_interp(const double *x, const double *y, const double *w,
                      double *out, int n, double span, int degree, int nv)
{
    int h = (int)floor(span * n);
    if (h < degree + 1) h = degree + 1;
    if (h > n) h = n;

    /* Auto-select nv to match R's behavior */
    if (nv <= 0) {
        /* R uses nvmax = max(200, n), and the kdtree produces ~ceil(log2(n))+1 leaves */
        /* For n=41, nv=9. For larger n, more vertices. */
        nv = 9;  /* Good default for typical quality score ranges */
        if (n > 100) nv = 17;
        if (n > 200) nv = 33;
    }
    if (nv > n) nv = n;

    /* Find data range with 0.5% margin (matching R) */
    double xmin = x[0], xmax = x[0];
    for (int i = 1; i < n; i++) {
        if (x[i] < xmin) xmin = x[i];
        if (x[i] > xmax) xmax = x[i];
    }
    double margin = 0.005 * (xmax - xmin);
    double lo = xmin - margin;
    double hi = xmax + margin;

    /* Evaluate at vertices */
    double *vx = (double *)malloc(nv * sizeof(double));
    double *vvals = (double *)malloc(nv * sizeof(double));
    double *vderivs = (double *)malloc(nv * sizeof(double));

    for (int i = 0; i < nv; i++) {
        vx[i] = lo + (hi - lo) * i / (nv - 1);
        loess_eval_one(x, y, w, n, vx[i], h, degree, &vvals[i], &vderivs[i]);
    }

    /* Hermite cubic interpolation at each output point */
    for (int i = 0; i < n; i++) {
        double xi = x[i];
        /* Find enclosing interval */
        int j = 0;
        while (j < nv - 2 && vx[j + 1] < xi) j++;

        double x0 = vx[j], x1 = vx[j + 1];
        double H = x1 - x0;
        if (H == 0.0) { out[i] = vvals[j]; continue; }
        double t = (xi - x0) / H;

        /* Hermite basis functions */
        double phi0 = (1-t)*(1-t) * (1 + 2*t);
        double phi1 = t*t * (3 - 2*t);
        double psi0 = t * (1-t)*(1-t);
        double psi1 = t*t * (t - 1);

        out[i] = phi0*vvals[j] + phi1*vvals[j+1]
               + (psi0*vderivs[j] + psi1*vderivs[j+1]) * H;
    }

    free(vx);
    free(vvals);
    free(vderivs);
}
