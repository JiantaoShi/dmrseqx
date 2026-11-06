/*
 * Local quadratic regression with evaluation on an adaptive tree, the default
 * evaluation structure of local regression described by Loader (1999, Local
 * Regression and Likelihood, Springer, chapter 12).
 *
 * The first vertices are the smallest and the largest position. A cell [a, b] is
 * split at its midpoint when (b - a) / min(h(a), h(b)) > cut, where h is the
 * bandwidth at a vertex. At a vertex x0 the bandwidth is max(h, d_k), with d_k
 * the distance to the k-th nearest point and k = n * nn (truncated). Points get
 * the weight w * (1 - |u|^3)^3, u = (x - x0) / bandwidth, and a quadratic in u
 * is fitted, which gives a value and a slope at x0. Between vertices the fit is
 * the cubic Hermite interpolant of the values and slopes; positions outside the
 * outermost vertices extrapolate linearly from the outermost vertex.
 *
 * The local fit starts from the weighted mean of the points inside the bandwidth,
 * and the normal equations for the correction are solved after scaling them to a
 * unit diagonal. When they are nearly singular (few points inside the bandwidth),
 * directions with negligible eigenvalues are dropped, so that the fit stays at the
 * start in those directions.
 */
#include "dmrseqx.h"
#include <math.h>

typedef struct {
    const double *x, *y, *w;
    long n, k;
    double h;
} fit_data_t;

typedef struct { double x, h, val, slope; } vertex_t;

static long lower_bound(const double *x, long n, double v)
{
    long lo = 0, hi = n;
    while (lo < hi) {
        const long mid = lo + (hi - lo) / 2;
        if (x[mid] < v) lo = mid + 1; else hi = mid;
    }
    return lo;
}

/* Eigenvalues and vectors of a symmetric 3 x 3 matrix by cyclic Jacobi rotations:
 * a is overwritten, its diagonal holding the eigenvalues, the columns of q the
 * eigenvectors. */
static void jacobi3(double a[3][3], double q[3][3])
{
    int sweep, p, r, i;
    for (p = 0; p < 3; p++) for (r = 0; r < 3; r++) q[p][r] = (p == r);
    for (sweep = 0; sweep < 30; sweep++) {
        const double off = fabs(a[0][1]) + fabs(a[0][2]) + fabs(a[1][2]);
        if (off < 1e-300) return;
        for (p = 0; p < 2; p++) {
            for (r = p + 1; r < 3; r++) {
                double theta, t, c, s;
                if (a[p][r] == 0) continue;
                theta = (a[r][r] - a[p][p]) / (2 * a[p][r]);
                t = (theta >= 0 ? 1.0 : -1.0) / (fabs(theta) + sqrt(theta * theta + 1));
                c = 1 / sqrt(t * t + 1);
                s = t * c;
                for (i = 0; i < 3; i++) {
                    const double aip = a[i][p], air = a[i][r];
                    a[i][p] = c * aip - s * air;
                    a[i][r] = s * aip + c * air;
                }
                for (i = 0; i < 3; i++) {
                    const double api = a[p][i], ari = a[r][i];
                    a[p][i] = c * api - s * ari;
                    a[r][i] = s * api + c * ari;
                }
                for (i = 0; i < 3; i++) {
                    const double qip = q[i][p], qir = q[i][r];
                    q[i][p] = c * qip - s * qir;
                    q[i][r] = s * qip + c * qir;
                }
            }
        }
    }
}

/* Coefficients of the local quadratic from the moments of the kernel weights. */
static void local_solve(const double mom[5], const double rhs[3], double beta[3])
{
    const double ybar = rhs[0] / mom[0];
    double dg[3], z[3][3], b[3], chol[3][3];
    int i, j, k, ok = 1;
    for (i = 0; i < 3; i++) dg[i] = sqrt(mom[2 * i]);
    for (i = 0; i < 3; i++) {
        b[i] = dg[i] > 0 ? (rhs[i] - ybar * mom[i]) / dg[i] : 0.0;
        for (j = 0; j < 3; j++)
            z[i][j] = (dg[i] > 0 && dg[j] > 0) ? mom[i + j] / (dg[i] * dg[j]) : 0.0;
    }
    /* Well conditioned: Cholesky of the scaled matrix */
    for (i = 0; i < 3 && ok; i++) {
        for (j = 0; j <= i; j++) {
            double sum = z[i][j];
            for (k = 0; k < j; k++) sum -= chol[i][k] * chol[j][k];
            if (i == j) {
                if (!(sum > 1e-4)) { ok = 0; break; }
                chol[i][i] = sqrt(sum);
            } else {
                chol[i][j] = sum / chol[j][j];
            }
        }
    }
    if (ok) {
        double e[3];
        for (i = 0; i < 3; i++) {
            double sum = b[i];
            for (k = 0; k < i; k++) sum -= chol[i][k] * e[k];
            e[i] = sum / chol[i][i];
        }
        for (i = 2; i >= 0; i--) {
            double sum = e[i];
            for (k = i + 1; k < 3; k++) sum -= chol[k][i] * beta[k];
            beta[i] = sum / chol[i][i];
        }
    } else {
        double q[3][3], w[3], tol = 0;
        jacobi3(z, q);
        for (i = 0; i < 3; i++) if (z[i][i] > tol) tol = z[i][i];
        tol *= 1e-8;
        for (i = 0; i < 3; i++) {
            double sum = 0;
            for (j = 0; j < 3; j++) sum += q[j][i] * b[j];
            w[i] = z[i][i] > tol ? sum / z[i][i] : 0.0;
        }
        for (i = 0; i < 3; i++) {
            double sum = 0;
            for (j = 0; j < 3; j++) sum += q[i][j] * w[j];
            beta[i] = sum;
        }
    }
    for (i = 0; i < 3; i++) beta[i] = dg[i] > 0 ? beta[i] / dg[i] : 0.0;
    beta[0] += ybar;
}

static void fit_vertex(const fit_data_t *d, double x0, vertex_t *out)
{
    long right = lower_bound(d->x, d->n, x0), left = right - 1, i, first, last;
    double dk = 0, bw, mom[5] = {0, 0, 0, 0, 0}, rhs[3] = {0, 0, 0}, beta[3];
    for (i = 0; i < d->k; i++) {
        if (left >= 0 && (right >= d->n || x0 - d->x[left] <= d->x[right] - x0)) {
            dk = x0 - d->x[left];
            left--;
        } else {
            dk = d->x[right] - x0;
            right++;
        }
    }
    bw = dk > d->h ? dk : d->h;
    first = lower_bound(d->x, d->n, x0 - bw);
    last = lower_bound(d->x, d->n, x0 + bw);
    for (i = first; i < last; i++) {
        const double u = (d->x[i] - x0) / bw, au = fabs(u);
        double ker, wu, u2;
        if (au >= 1.0) continue;
        ker = 1.0 - au * au * au;
        wu = d->w[i] * ker * ker * ker;
        u2 = u * u;
        mom[0] += wu;
        mom[1] += wu * u;
        mom[2] += wu * u2;
        mom[3] += wu * u2 * u;
        mom[4] += wu * u2 * u2;
        rhs[0] += wu * d->y[i];
        rhs[1] += wu * d->y[i] * u;
        rhs[2] += wu * d->y[i] * u2;
    }
    out->x = x0;
    out->h = bw;
    out->val = 0;
    out->slope = 0;
    if (!(mom[0] > 0)) return;
    local_solve(mom, rhs, beta);
    out->val = beta[0];
    out->slope = beta[1] / bw;
}

static void grow(const fit_data_t *d, vertex_t a, vertex_t b, double cut,
                 const double *xEval, double *fitted, long lo, long hi)
{
    const double hmin = a.h < b.h ? a.h : b.h;
    if ((b.x - a.x) / hmin > cut) {
        const double mid = (a.x + b.x) / 2;
        const long m = lower_bound(xEval + lo, hi - lo, mid) + lo;
        vertex_t v;
        fit_vertex(d, mid, &v);
        grow(d, a, v, cut, xEval, fitted, lo, m);
        grow(d, v, b, cut, xEval, fitted, m, hi);
    } else {
        const double len = b.x - a.x;
        long i;
        for (i = lo; i < hi; i++) {
            const double u = (xEval[i] - a.x) / len, v = 1 - u;
            if (u < 0) {
                fitted[i] = a.val + (xEval[i] - a.x) * a.slope;
            } else if (u > 1) {
                fitted[i] = b.val + (xEval[i] - b.x) * b.slope;
            } else {
                const double phi0 = v * v * (3 - 2 * v), phi1 = u * u * (3 - 2 * u);
                const double psi0 = v * v * u, psi1 = -u * u * v;
                fitted[i] = phi0 * a.val + phi1 * b.val + len * (psi0 * a.slope + psi1 * b.slope);
            }
        }
    }
}

int smooth_block(const double *xFit, const double *yFit, const double *wFit, long nFit,
                 double nn, double h, double cut,
                 const double *xEval, long nEval, double *fitted)
{
    fit_data_t d;
    vertex_t a, b;
    d.x = xFit;
    d.y = yFit;
    d.w = wFit;
    d.n = nFit;
    d.h = h;
    /* tolerance so that, for example, nn = 30 / n gives k = 30 */
    d.k = (long) (nFit * nn + 1e-8);
    if (d.k < 1) d.k = 1;
    if (d.k > d.n) d.k = d.n;
    fit_vertex(&d, xFit[0], &a);
    fit_vertex(&d, xFit[nFit - 1], &b);
    grow(&d, a, b, cut, xEval, fitted, 0, nEval);
    return 1;
}
