/*
 * The region level model of dmrseq (Korthauer et al., Biostatistics 2019),
 *
 *   Z ~ covariate of interest + factor(locus) + adjustment covariates,
 *   Z = asin(2 (M + 0.05) / (N + 0.1) - 1),
 *
 * with variance proportional to 1 / MedCov (MedCov is the median coverage of a
 * locus) and errors correlated within samples: AR(1) over rows for regions of at
 * least 40 loci, and continuous AR(1) over base pairs (rho^distance) for shorter
 * regions. It is fitted by restricted maximum likelihood.
 *
 * The inverse of the correlation matrix of a sample only connects neighbouring
 * rows, u' R^-1 v = u_1 v_1 + sum_t (u_t - r_t u_t-1) (v_t - r_t v_t-1) / (1 - r_t^2),
 * so that the likelihood is cheap at any value of the correlation parameter:
 *
 * - Every sample covers every locus: the locus effects can be profiled out in
 *   closed form, and an evaluation costs O(L), or O(1) for AR(1).
 * - Otherwise: the normal equations of the locus effects are sparse, and are
 *   eliminated with a Cholesky factorization on their envelope.
 */
#include "dmrseqx.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/*
 * Stop rule of fit = "dmrseq". gls() starts corCAR1 at rho = 0.2, and nlminb stops
 * at once where the objective is flat there, which happens when CpGs are about 8
 * bp or more apart. Comparing with gls() on real regions, the outcome is
 * determined by g0^2 / (2 |f0|), f0 and g0 being the objective and its derivative
 * at the start: below 1e-11 gls() always stops, above 1e-10 it nearly always
 * moves, and between the two the rounding errors of the objective decide. The
 * value below gives as many gained as lost calls.
 */
#define STOP_RULE_TOL 3e-11
#define STOP_RULE_STEP 1e-2

typedef struct {
    int car1, n, nLoci, nSamples, p, dfResid, balanced;
    double sumLogW, constant;
    fit_ws_t *ws;
    /* results of the last evaluation */
    double rss, cQc;
} model_t;

struct fit_ws {
    /* a region */
    size_t capL, capN, capRows, capEnv;
    int capP;
    int *cv, *mt, *rowLocus, *rowSample, *locusMap, *envFirst, *gapIdx;
    size_t *envPtr;
    unsigned char *rowFirst;
    double *medCov, *medTmp, *z, *w, *rowGap, *gapU, *phiTab, *sTab, *logTab;
    double *c, *gapL, *y0, *A, *Bl, *Wk, *env, *Bm, *ybar;
    int nGapU;
    /* covariates */
    double *G, *Ug, *v, *r, *Sm, *gRow, *cholS, *xy, *xbar, *e, *Fm, *Fu, *Fg, *Fprev;
    int capF;
    /* AR(1) sums of the balanced case */
    double cc0, ccCur, ccPrev, ccCross, A0, ACur, APrev, BSum;
    double *v0, *vCur, *vPrev, *vCross;
    const double *X;
};

fit_ws_t *fit_ws_new(void)
{
    return (fit_ws_t *) calloc(1, sizeof(fit_ws_t));
}

void fit_ws_free(fit_ws_t *ws)
{
    if (!ws) return;
    free(ws->cv); free(ws->mt); free(ws->rowLocus); free(ws->rowSample); free(ws->locusMap);
    free(ws->envFirst); free(ws->gapIdx); free(ws->envPtr); free(ws->rowFirst);
    free(ws->medCov); free(ws->medTmp); free(ws->z); free(ws->w); free(ws->rowGap);
    free(ws->gapU); free(ws->phiTab); free(ws->sTab); free(ws->logTab);
    free(ws->c); free(ws->gapL); free(ws->y0); free(ws->A); free(ws->Bl); free(ws->Wk);
    free(ws->env); free(ws->Bm); free(ws->ybar);
    free(ws->G); free(ws->Ug); free(ws->v); free(ws->r); free(ws->Sm); free(ws->gRow);
    free(ws->cholS); free(ws->xy); free(ws->xbar); free(ws->e);
    free(ws->Fm); free(ws->Fu); free(ws->Fg); free(ws->Fprev);
    free(ws->v0); free(ws->vCur); free(ws->vPrev); free(ws->vCross);
    free(ws);
}

#define GROW(ptr, type, count) do { \
        void *tmp_ = realloc((ptr), sizeof(type) * (size_t) (count)); \
        if (!tmp_) return 0; \
        (ptr) = (type *) tmp_; \
    } while (0)

static int ws_reserve(fit_ws_t *ws, int L, int N, int p)
{
    const size_t rows = (size_t) L * N;
    const int m = p + 1;
    if ((size_t) L > ws->capL || (size_t) N > ws->capN || p > ws->capP) {
        const size_t capL = (size_t) L > ws->capL ? (size_t) L * 2 : ws->capL;
        const size_t capN = (size_t) N > ws->capN ? (size_t) N : ws->capN;
        const int capP = p > ws->capP ? p : ws->capP;
        const int capM = capP + 1;
        ws->capL = capL;
        ws->capN = capN;
        ws->capP = capP;
        GROW(ws->medTmp, double, capN);
        GROW(ws->medCov, double, capL);
        GROW(ws->locusMap, int, capL);
        GROW(ws->envFirst, int, capL);
        GROW(ws->envPtr, size_t, capL + 1);
        GROW(ws->c, double, capL);
        GROW(ws->gapL, double, capL);
        GROW(ws->A, double, capL);
        GROW(ws->Bl, double, capL);
        GROW(ws->ybar, double, capL);
        GROW(ws->Wk, double, capL * capP);
        GROW(ws->Bm, double, capL * capM);
        GROW(ws->G, double, capN * capP);
        GROW(ws->Ug, double, (size_t) capP * capP);
        GROW(ws->cholS, double, (size_t) capP * capP);
        GROW(ws->v, double, capP);
        GROW(ws->r, double, capP);
        GROW(ws->xy, double, capP);
        GROW(ws->xbar, double, capP);
        GROW(ws->e, double, capP);
        GROW(ws->v0, double, capP);
        GROW(ws->vCur, double, capP);
        GROW(ws->vPrev, double, capP);
        GROW(ws->vCross, double, capP);
        GROW(ws->gRow, double, capM);
        GROW(ws->Sm, double, (size_t) capM * capM);
    }
    (void) m;
    if (rows > ws->capRows) {
        const size_t cap = rows * 2;
        ws->capRows = cap;
        GROW(ws->cv, int, cap);
        GROW(ws->mt, int, cap);
        GROW(ws->rowLocus, int, cap);
        GROW(ws->rowSample, int, cap);
        GROW(ws->gapIdx, int, cap);
        GROW(ws->rowFirst, unsigned char, cap);
        GROW(ws->z, double, cap);
        GROW(ws->w, double, cap);
        GROW(ws->rowGap, double, cap);
        GROW(ws->y0, double, cap);
        GROW(ws->gapU, double, cap);
        GROW(ws->phiTab, double, cap);
        GROW(ws->sTab, double, cap);
        GROW(ws->logTab, double, cap);
    }
    return 1;
}

/* ---------------------------------------------------------------- likelihood */

/* Upper Cholesky factor of the p x p matrix A (column-major, upper part read);
 * returns log determinant through *logDet, or 0 if singular. */
static int chol_upper(const double *A, int p, int ld, double *U, double *logDet)
{
    int k, q, j;
    double ld_sum = 0;
    for (k = 0; k < p; k++) {
        double d = A[k + k * ld];
        const double d0 = d;
        for (q = 0; q < k; q++) {
            double e = A[q + k * ld];
            for (j = 0; j < q; j++) e -= U[j + q * p] * U[j + k * p];
            e /= U[q + q * p];
            U[q + k * p] = e;
            d -= e * e;
        }
        if (!(d > 1e-10 * d0) || !(d > 0)) return 0;
        d = sqrt(d);
        U[k + k * p] = d;
        ld_sum += 2 * log(d);
    }
    *logDet = ld_sum;
    return 1;
}

static double lik_balanced(model_t *m, double phi)
{
    fit_ws_t *ws = m->ws;
    const int L = m->nLoci, p = m->p;
    const double *c = ws->c, *A = ws->A, *Bl = ws->Bl, *W = ws->Wk;
    double *v = ws->v, *r = ws->r;
    double cQc, rss0, sumLog = 0, rr = 0;
    int k, q, i;

    if (m->car1) {
        cQc = c[0] * c[0];
        rss0 = A[0];
        for (k = 0; k < p; k++) v[k] = c[0] * W[k * L];
        for (i = 1; i < L; i++) {
            const double ph = pow(phi, ws->gapL[i]);
            const double d = 1 - ph * ph;
            double ci;
            if (!(d > 0)) return -HUGE_VAL;
            sumLog += log(d);
            ci = c[i] - ph * c[i - 1];
            cQc += ci * ci / d;
            rss0 += (A[i] - 2 * ph * Bl[i] + ph * ph * A[i - 1]) / d;
            for (k = 0; k < p; k++) v[k] += ci * (W[i + k * L] - ph * W[i - 1 + k * L]) / d;
        }
    } else {
        const double d = 1 - phi * phi;
        if (!(d > 0)) return -HUGE_VAL;
        sumLog = (L - 1) * log(d);
        cQc = ws->cc0 + (ws->ccCur - 2 * phi * ws->ccCross + phi * phi * ws->ccPrev) / d;
        rss0 = ws->A0 + (ws->ACur - 2 * phi * ws->BSum + phi * phi * ws->APrev) / d;
        for (k = 0; k < p; k++)
            v[k] = ws->v0[k] + (ws->vCur[k] - phi * ws->vCross[k] + phi * phi * ws->vPrev[k]) / d;
    }
    for (k = 0; k < p; k++) {
        double e = v[k];
        for (q = 0; q < k; q++) e -= ws->Ug[q + k * p] * r[q];
        e /= ws->Ug[k + k * p];
        r[k] = e;
        rr += e * e;
    }
    m->rss = rss0 - rr / cQc;
    m->cQc = cQc;
    if (!(m->rss > 0) || !(cQc > 0)) return -HUGE_VAL;
    return -0.5 * m->dfResid * log(m->rss) - 0.5 * (m->nSamples - 1) * sumLog
        - 0.5 * p * log(cQc) + m->constant;
}

static double lik_general(model_t *m, double phi)
{
    fit_ws_t *ws = m->ws;
    const int n = m->n, L = m->nLoci, p = m->p, mm = p + 1, N = m->nSamples;
    const double *z = ws->z, *w = ws->w, *X = ws->X;
    const int *samp = ws->rowSample, *locus = ws->rowLocus;
    const int *first = ws->envFirst;
    const size_t *ptr = ws->envPtr;
    double *env = ws->env, *B = ws->Bm, *S = ws->Sm, *g = ws->gRow;
    double logR = 0, logDet = 0;
    int t, k, q, i, j;

    for (i = 0; i < ws->nGapU; i++) {
        const double ph = m->car1 ? pow(phi, ws->gapU[i]) : phi;
        const double d = 1 - ph * ph;
        if (!(d > 0)) return -HUGE_VAL;
        ws->phiTab[i] = ph;
        ws->sTab[i] = sqrt(d);
        ws->logTab[i] = log(d);
    }
    memset(env, 0, sizeof(double) * ptr[L]);
    memset(B, 0, sizeof(double) * (size_t) L * mm);
    memset(S, 0, sizeof(double) * (size_t) mm * mm);

    for (t = 0; t < n; t++) {
        const int li = locus[t], si = samp[t];
        double a, c = 0;
        if (ws->rowFirst[t]) {
            a = w[t];
            for (k = 0; k < p; k++) g[k] = w[t] * X[si + k * N];
            g[p] = w[t] * z[t];
        } else {
            const int gi = ws->gapIdx[t];
            const double ph = ws->phiTab[gi], s = ws->sTab[gi];
            double d;
            logR += ws->logTab[gi];
            a = w[t] / s;
            c = ph * w[t - 1] / s;
            d = a - c;
            for (k = 0; k < p; k++) g[k] = d * X[si + k * N];
            g[p] = (w[t] * z[t] - ph * w[t - 1] * z[t - 1]) / s;
        }
        env[ptr[li] + (li - first[li])] += a * a;
        for (k = 0; k < mm; k++) B[li + k * L] += a * g[k];
        if (!ws->rowFirst[t]) {
            const int h = locus[t - 1];
            env[ptr[h] + (h - first[h])] += c * c;
            env[ptr[li] + (h - first[li])] -= a * c;
            for (k = 0; k < mm; k++) B[h + k * L] -= c * g[k];
        }
        for (k = 0; k < mm; k++)
            for (q = 0; q <= k; q++) S[q + k * mm] += g[q] * g[k];
    }

    /* Cholesky on the envelope, by rows */
    for (i = 0; i < L; i++) {
        double *ri = env + ptr[i];
        const int fi = first[i];
        double d = ri[i - fi];
        for (j = fi; j < i; j++) {
            const double *rj = env + ptr[j];
            const int fj = first[j];
            const int from = fi > fj ? fi : fj;
            double e = ri[j - fi];
            int h;
            for (h = from; h < j; h++) e -= ri[h - fi] * rj[h - fj];
            e /= rj[j - fj];
            ri[j - fi] = e;
            d -= e * e;
        }
        if (!(d > 0)) return -HUGE_VAL;
        d = sqrt(d);
        ri[i - fi] = d;
        logDet += 2 * log(d);
    }
    /* forward solves */
    for (k = 0; k < mm; k++) {
        double *col = B + k * L;
        for (i = 0; i < L; i++) {
            const double *ri = env + ptr[i];
            const int fi = first[i];
            double e = col[i];
            for (j = fi; j < i; j++) e -= ri[j - fi] * col[j];
            col[i] = e / ri[i - fi];
        }
    }
    for (k = 0; k < mm; k++)
        for (q = 0; q <= k; q++) {
            double s = 0;
            for (j = 0; j < L; j++) s += B[j + q * L] * B[j + k * L];
            S[q + k * mm] -= s;
        }
    /* Schur complement */
    {
        double ld = 0;
        if (!chol_upper(S, p, mm, ws->cholS, &ld)) return -HUGE_VAL;
        logDet += ld;
    }
    m->rss = S[p + p * mm];
    for (k = 0; k < p; k++) {
        double e = S[k + p * mm];
        for (q = 0; q < k; q++) e -= ws->cholS[q + k * p] * ws->xy[q];
        e /= ws->cholS[k + k * p];
        ws->xy[k] = e;
        m->rss -= e * e;
    }
    if (!(m->rss > 0)) return -HUGE_VAL;
    return -0.5 * m->dfResid * log(m->rss) - 0.5 * logR - 0.5 * logDet;
}

static double likelihood(model_t *m, double phi)
{
    return m->balanced ? lik_balanced(m, phi) : lik_general(m, phi);
}

/* ------------------------------------------------------------- optimization */

/* Parameterization of nlme: logit(phi) for CAR(1), log((1 + rho) / (1 - rho)) for AR(1) */
static double theta_to_phi(const model_t *m, double theta)
{
    const double e = exp(theta);
    return m->car1 ? e / (1 + e) : (e - 1) / (e + 1);
}

static double lik_theta(model_t *m, double theta)
{
    const double v = likelihood(m, theta_to_phi(m, theta));
    return isfinite(v) ? v : -1e300;
}

/* Brent's method for the maximum on [ax, bx] */
static double brent_max(model_t *m, double ax, double bx, double tol)
{
    const double golden = (3.0 - sqrt(5.0)) * 0.5, eps = 1.490116e-08;
    double a = ax, b = bx, v, w, x, d = 0, e = 0, u, fu, fv, fw, fx, xm, tol1, t2;
    int iter;
    v = w = x = a + golden * (b - a);
    fv = fw = fx = -lik_theta(m, x);
    for (iter = 0; iter < 200; iter++) {
        double p = 0, q = 0, r = 0;
        xm = (a + b) * 0.5;
        tol1 = eps * fabs(x) + tol / 3;
        t2 = tol1 * 2;
        if (fabs(x - xm) <= t2 - (b - a) * 0.5) break;
        if (fabs(e) > tol1) {
            r = (x - w) * (fx - fv);
            q = (x - v) * (fx - fw);
            p = (x - v) * q - (x - w) * r;
            q = (q - r) * 2;
            if (q > 0) p = -p; else q = -q;
            r = e;
            e = d;
        }
        if (fabs(p) >= fabs(q * 0.5 * r) || p <= q * (a - x) || p >= q * (b - x)) {
            e = (x < xm) ? b - x : a - x;
            d = golden * e;
        } else {
            d = p / q;
            u = x + d;
            if (u - a < t2 || b - u < t2) d = (x >= xm) ? -tol1 : tol1;
        }
        u = (fabs(d) >= tol1) ? x + d : ((d > 0) ? x + tol1 : x - tol1);
        fu = -lik_theta(m, u);
        if (fu <= fx) {
            if (u < x) b = x; else a = x;
            v = w; w = x; x = u;
            fv = fw; fw = fx; fx = fu;
        } else {
            if (u < x) a = u; else b = u;
            if (fu <= fw || w == x) {
                v = w; fv = fw;
                w = u; fw = fu;
            } else if (fu <= fv || v == x || v == w) {
                v = u; fv = fu;
            }
        }
    }
    return x;
}

static double search_global(model_t *m)
{
    const double lo = m->car1 ? -6.0 : -9.0, hi = m->car1 ? 14.0 : 9.0, step = m->car1 ? 1.0 : 1.5;
    double th, bestTh = lo, best = -HUGE_VAL, opt, v;
    for (th = lo; th <= hi + 1e-9; th += step) {
        v = lik_theta(m, th);
        if (v > best) { best = v; bestTh = th; }
    }
    opt = brent_max(m, bestTh - step, bestTh + step, 1e-6);
    v = lik_theta(m, opt);
    return v >= best ? opt : bestTh;
}

/* What gls() does for short regions: stop at the start where the objective is
 * flat, otherwise find the local optimum on the uphill side of the start. */
static double search_from_start(model_t *m)
{
    const double th0 = log(0.25), step = STOP_RULE_STEP;
    const double l0 = lik_theta(m, th0), lp = lik_theta(m, th0 + step), lm = lik_theta(m, th0 - step);
    const double g0 = (lm - lp) / (2 * step), f0 = -(l0 + m->sumLogW);
    double prev = th0, best = th0, v = l0, delta = 0.5, next, opt, dir, lo, hi;

    if (g0 * g0 / (2 * fabs(f0)) < STOP_RULE_TOL) return th0;
    dir = lp >= lm ? 1.0 : -1.0;
    for (;;) {
        double nv;
        next = best + dir * delta;
        nv = lik_theta(m, next);
        if (!(nv > v) || fabs(next) > 30) break;
        prev = best;
        best = next;
        v = nv;
        delta *= 2;
    }
    lo = prev < next ? prev : next;
    hi = prev < next ? next : prev;
    opt = brent_max(m, lo, hi, 1e-6);
    return lik_theta(m, opt) >= v ? opt : best;
}

/* ------------------------------------------------------------------- set up */

static int setup_covariates(model_t *m, const double *X)
{
    fit_ws_t *ws = m->ws;
    const int N = m->nSamples, p = m->p;
    int s, k, q;
    double *G = ws->G, ld = 0;
    double *GG = ws->Sm;
    for (k = 0; k < p; k++) {
        double sum = 0;
        for (s = 0; s < N; s++) sum += X[s + k * N];
        sum /= N;
        ws->xbar[k] = sum;
        for (s = 0; s < N; s++) G[s + k * N] = X[s + k * N] - sum;
    }
    for (k = 0; k < p; k++)
        for (q = 0; q <= k; q++) {
            double sum = 0;
            for (s = 0; s < N; s++) sum += G[s + q * N] * G[s + k * N];
            GG[q + k * p] = sum;
        }
    if (!chol_upper(GG, p, p, ws->Ug, &ld)) return 0;
    m->constant = -0.5 * ld;
    return 1;
}

/* z (L x N) and c (L) are set */
static int setup_balanced(model_t *m, const double *X)
{
    fit_ws_t *ws = m->ws;
    const int L = m->nLoci, N = m->nSamples, p = m->p;
    const double *c = ws->c, *z = ws->z;
    double *y0 = ws->y0, *ybar = ws->ybar, *A = ws->A, *Bl = ws->Bl, *W = ws->Wk;
    double sumLogC = 0;
    int i, s, k;

    if (!setup_covariates(m, X)) return 0;
    for (i = 0; i < L; i++) {
        double sum = 0;
        for (s = 0; s < N; s++) sum += z[i + (size_t) s * L];
        ybar[i] = c[i] * (sum / N);
        sumLogC += log(c[i]);
    }
    memset(A, 0, sizeof(double) * L);
    memset(Bl, 0, sizeof(double) * L);
    memset(W, 0, sizeof(double) * (size_t) L * p);
    for (s = 0; s < N; s++) {
        const double *zs = z + (size_t) s * L;
        double *ys = y0 + (size_t) s * L;
        for (i = 0; i < L; i++) ys[i] = c[i] * zs[i] - ybar[i];
        A[0] += ys[0] * ys[0];
        for (i = 1; i < L; i++) {
            A[i] += ys[i] * ys[i];
            Bl[i] += ys[i] * ys[i - 1];
        }
        for (k = 0; k < p; k++) {
            const double gsk = ws->G[s + k * N];
            double *Wk = W + (size_t) k * L;
            for (i = 0; i < L; i++) Wk[i] += gsk * ys[i];
        }
    }
    m->constant += -0.5 * L * log((double) N) - sumLogC;
    m->sumLogW = N * sumLogC;

    if (!m->car1) {
        ws->cc0 = c[0] * c[0];
        ws->A0 = A[0];
        ws->ccCur = ws->ccPrev = ws->ccCross = ws->ACur = ws->APrev = ws->BSum = 0;
        for (k = 0; k < p; k++) {
            const double *Wk = W + (size_t) k * L;
            ws->v0[k] = c[0] * Wk[0];
            ws->vCur[k] = ws->vPrev[k] = ws->vCross[k] = 0;
            for (i = 1; i < L; i++) {
                ws->vCur[k] += c[i] * Wk[i];
                ws->vPrev[k] += c[i - 1] * Wk[i - 1];
                ws->vCross[k] += c[i] * Wk[i - 1] + c[i - 1] * Wk[i];
            }
        }
        for (i = 1; i < L; i++) {
            ws->ccCur += c[i] * c[i];
            ws->ccPrev += c[i - 1] * c[i - 1];
            ws->ccCross += c[i] * c[i - 1];
            ws->ACur += A[i];
            ws->APrev += A[i - 1];
            ws->BSum += Bl[i];
        }
    }
    return 1;
}

static int cmp_double(const void *a, const void *b)
{
    const double x = *(const double *) a, y = *(const double *) b;
    return (x > y) - (x < y);
}

static int find_gap(const double *u, int n, double v)
{
    int lo = 0, hi = n - 1;
    while (lo < hi) {
        const int mid = lo + (hi - lo) / 2;
        if (u[mid] < v) lo = mid + 1; else hi = mid;
    }
    return lo;
}

/* The rows z, w, rowSample, rowLocus, rowFirst, rowGap are set */
static int setup_general(model_t *m)
{
    fit_ws_t *ws = m->ws;
    const int n = m->n, L = m->nLoci;
    int t, i;
    size_t total = 0;

    /* unique distances between neighbouring rows */
    if (m->car1) {
        int nu = 0;
        for (t = 0; t < n; t++) if (!ws->rowFirst[t]) ws->gapU[nu++] = ws->rowGap[t];
        if (nu > 0) {
            int w = 1;
            qsort(ws->gapU, nu, sizeof(double), cmp_double);
            for (i = 1; i < nu; i++) if (ws->gapU[i] != ws->gapU[w - 1]) ws->gapU[w++] = ws->gapU[i];
            nu = w;
        }
        ws->nGapU = nu;
        for (t = 0; t < n; t++)
            ws->gapIdx[t] = ws->rowFirst[t] ? 0 : find_gap(ws->gapU, nu, ws->rowGap[t]);
    } else {
        ws->nGapU = 1;
        ws->gapU[0] = 1;
        for (t = 0; t < n; t++) ws->gapIdx[t] = 0;
    }
    /* envelope of the locus effects */
    for (i = 0; i < L; i++) ws->envFirst[i] = i;
    for (t = 0; t < n; t++)
        if (!ws->rowFirst[t] && ws->rowLocus[t - 1] < ws->envFirst[ws->rowLocus[t]])
            ws->envFirst[ws->rowLocus[t]] = ws->rowLocus[t - 1];
    for (i = 0; i < L; i++) {
        ws->envPtr[i] = total;
        total += (size_t) (i - ws->envFirst[i] + 1);
    }
    ws->envPtr[L] = total;
    if (total > ws->capEnv) {
        GROW(ws->env, double, total * 2);
        ws->capEnv = total * 2;
    }
    m->constant = 0;
    m->sumLogW = 0;
    for (t = 0; t < n; t++) m->sumLogW += log(ws->w[t]);
    return 1;
}

/* Coefficient of the first covariate and its standard error at the last
 * evaluation of the likelihood */
static void estimate(model_t *m, double *beta, double *se)
{
    fit_ws_t *ws = m->ws;
    const int p = m->p;
    const double sigma2 = m->rss / m->dfResid;
    const double *U = m->balanced ? ws->Ug : ws->cholS;
    const double *rhs = m->balanced ? ws->r : ws->xy;
    const double scale = m->balanced ? m->cQc : 1.0;
    double *e = ws->e, sum = 0, b = 0;
    int k, q;

    /* first row of the inverse of U: forward substitution on U' x = e_1 */
    for (k = 0; k < p; k++) {
        double v = (k == 0) ? 1.0 : 0.0;
        for (q = 0; q < k; q++) v -= U[q + k * p] * e[q];
        e[k] = v / U[k + k * p];
        sum += e[k] * e[k];
        b += e[k] * rhs[k];
    }
    *beta = b / scale;
    *se = sqrt(sigma2 * sum / scale);
}

/* --------------------------------------------------------------------- data */

static void sort_small(double *x, int n)
{
    int i, j;
    if (n > 40) {
        qsort(x, n, sizeof(double), cmp_double);
        return;
    }
    for (i = 1; i < n; i++) {
        const double v = x[i];
        for (j = i - 1; j >= 0 && x[j] > v; j--) x[j + 1] = x[j];
        x[j + 1] = v;
    }
}

/* Critical values at p = 0.01 of the test for two outliers on one side (Grubbs
 * 1950), for n = 4, ..., 10: small values indicate outliers */
static const double GRUBBS_TWO[7] = {1e-05, 0.0035, 0.0186, 0.044, 0.075, 0.1082, 0.1415};

/* Grubbs tests for one and for two outliers at p < 0.01 */
static int grubbs_outliers(const double *values, int n, const params_t *P, double *x)
{
    long double total = 0;
    double mean, ss = 0;
    int i, lowSide;
    if (n < 3 || n > 10) return 0;
    for (i = 0; i < n; i++) x[i] = values[i];
    sort_small(x, n);
    if (x[0] == x[n - 1]) return 0;
    for (i = 0; i < n; i++) total += x[i];
    mean = (double) (total / n);
    for (i = 0; i < n; i++) ss += (x[i] - mean) * (x[i] - mean);
    lowSide = (x[n - 1] - mean) < (mean - x[0]);
    {
        const double o = lowSide ? x[0] : x[n - 1];
        const double g2 = (o - mean) * (o - mean) / (ss / (n - 1));
        const double s = (g2 * n * (2 - n)) / (g2 * n - (double) (n - 1) * (n - 1));
        const double t = sqrt(s);
        if (isnan(t) || t > P->grubbsT[n]) return 1;
    }
    if (n > 3) {
        /* the variance without the two extreme values */
        const int from = lowSide ? 2 : 0, to = lowSide ? n : n - 2;
        long double t2 = 0;
        double m2, ss2 = 0, u;
        for (i = from; i < to; i++) t2 += x[i];
        m2 = (double) (t2 / (n - 2));
        for (i = from; i < to; i++) ss2 += (x[i] - m2) * (x[i] - m2);
        u = (ss2 / (n - 3)) / (ss / (n - 1)) * (n - 3) / (n - 1);
        if (u < GRUBBS_TWO[n - 4]) return 1;
    }
    return 0;
}

/* The coefficients of the first nTest covariates at the last evaluation of the
 * likelihood */
static void coefficients(model_t *m, int nTest, double *beta)
{
    fit_ws_t *ws = m->ws;
    const int p = m->p;
    const double *U = m->balanced ? ws->Ug : ws->cholS;
    const double *rhs = m->balanced ? ws->r : ws->xy;
    const double scale = m->balanced ? m->cQc : 1.0;
    double *e = ws->e;
    int k, q;
    for (k = p - 1; k >= 0; k--) {
        double v = rhs[k];
        for (q = k + 1; q < p; q++) v -= U[k + q * p] * e[q];
        e[k] = v / U[k + k * p];
    }
    for (k = 0; k < nTest; k++) beta[k] = e[k] / scale;
}

/*
 * The sequential F statistic of the first nTest covariates, as anova() of gls gives
 * for a first term: the reduction of the sum of squares from the intercept alone, on
 * the whitened scale and ignoring all later terms, over nTest and sigma^2.
 */
static double sequential_f(model_t *m, const double *X, int nTest, double phi)
{
    fit_ws_t *ws = m->ws;
    const int n = m->n, N = m->nSamples, d = nTest + 1, dm = d + 1;
    double *C = ws->Fm, *U = ws->Fu, *g = ws->Fg, *prev = ws->Fprev;
    double ss = 0, ld = 0, e[64];
    int t, j, k;
    memset(C, 0, sizeof(double) * (size_t) dm * dm);
    for (t = 0; t < n; t++) {
        const int si = ws->rowSample[t];
        const double wt = ws->w[t];
        double cur[64];
        cur[0] = wt;
        for (j = 0; j < nTest; j++) cur[j + 1] = wt * X[si + j * N];
        cur[d] = wt * ws->z[t];
        if (ws->rowFirst[t]) {
            for (j = 0; j < dm; j++) g[j] = cur[j];
        } else {
            const double ph = m->car1 ? pow(phi, ws->rowGap[t]) : phi;
            const double s = sqrt(1 - ph * ph);
            for (j = 0; j < dm; j++) g[j] = (cur[j] - ph * prev[j]) / s;
        }
        for (j = 0; j < dm; j++) prev[j] = cur[j];
        for (j = 0; j < dm; j++)
            for (k = 0; k <= j; k++) C[k + j * dm] += g[k] * g[j];
    }
    if (!chol_upper(C, d, dm, U, &ld)) return NAN;
    for (j = 0; j < d; j++) {
        double v = C[j + d * dm];
        for (k = 0; k < j; k++) v -= U[k + j * d] * e[k];
        e[j] = v / U[j + j * d];
        if (j > 0) ss += e[j] * e[j];
    }
    return ss / nTest / (m->rss / m->dfResid);
}

/* One attempt at the fit: returns 1 if it worked */
static int fit_attempt(model_t *m, const double *X, int car1, const params_t *P, int search,
                       double *phiOut, double *betaOut, double *statOut)
{
    double theta, phi, l, beta, se;
    const int nTest = P->nTest;
    m->car1 = car1;
    if (m->balanced) {
        if (!setup_balanced(m, X)) return 0;
    } else {
        if (!setup_general(m)) return 0;
    }
    if (!search) {
        phi = 0;
    } else if (car1 && !P->fitGlobal) {
        theta = search_from_start(m);
        phi = theta_to_phi(m, theta);
    } else {
        theta = search_global(m);
        phi = theta_to_phi(m, theta);
    }
    l = likelihood(m, phi);
    if (!isfinite(l)) return 0;
    estimate(m, &beta, &se);
    if (!isfinite(beta) || !(se > 0)) return 0;
    *phiOut = phi;
    if (nTest == 1) {
        betaOut[0] = beta;
        *statOut = beta / se;
    } else {
        double f;
        int k;
        /* before the sums below overwrite anything that the coefficients need */
        coefficients(m, nTest, betaOut);
        for (k = 0; k < nTest; k++) if (!isfinite(betaOut[k])) return 0;
        f = sequential_f(m, X, nTest, phi);
        if (!isfinite(f) || f < 0) return 0;
        *statOut = sqrt(f);
    }
    return 1;
}

/* Is the matrix with intercept, the covariate of interest and the position of full
 * rank on the rows with coverage? */
static int full_rank_rows(model_t *m, const double *X, int nTest, const double *posRow)
{
    fit_ws_t *ws = m->ws;
    const int n = m->n, N = m->nSamples, d = nTest + 2;
    int t, j, k;
    double *C = ws->Fm, *U = ws->Fu, ld = 0;
    double mean = 0;
    if (n < d + 1) return 0;
    for (t = 0; t < n; t++) mean += posRow[t];
    mean /= n;
    memset(C, 0, sizeof(double) * (size_t) d * d);
    for (t = 0; t < n; t++) {
        double row[64];
        row[0] = 1;
        for (j = 0; j < nTest; j++) row[j + 1] = X[ws->rowSample[t] + j * N];
        row[nTest + 1] = posRow[t] - mean;
        for (j = 0; j < d; j++)
            for (k = 0; k <= j; k++) C[k + j * d] += row[k] * row[j];
    }
    return chol_upper(C, d, d, U, &ld);
}

int fit_region(const dataset_t *ds, const params_t *P, size_t offset, const int *ix, int L,
               const double *X, int p, int weighted, rng_t *rng, fit_ws_t *ws, fit_out_t *out)
{
    const int N = ds->nSamples;
    const size_t nTotal = ds->nTotal;
    model_t m;
    int s, i, t, n = 0, anyMeth = 0, anyUnmeth = 0, medOK = 1, grubbs = 0, jitter = 0;
    int nPresentLoci = 0;

    out->stat = out->phi = NAN;
    for (i = 0; i < P->nTest; i++) out->beta[i] = NAN;
    out->constant = 0;
    if (!ws_reserve(ws, L, N, p)) return 0;
    if (P->nTest + 4 > ws->capF) {
        const int cap = P->nTest + 4;
        GROW(ws->Fm, double, (size_t) cap * cap);
        GROW(ws->Fu, double, (size_t) cap * cap);
        GROW(ws->Fg, double, cap);
        GROW(ws->Fprev, double, cap);
        ws->capF = cap;
    }
    ws->X = X;

    /* counts */
    for (s = 0; s < N; s++) {
        const int *Ms = ds->M + (size_t) s * nTotal + offset;
        const int *Cs = ds->Cov + (size_t) s * nTotal + offset;
        int *cv = ws->cv + (size_t) s * L, *mt = ws->mt + (size_t) s * L;
        for (i = 0; i < L; i++) {
            const int c = Cs[ix[i]], me = Ms[ix[i]];
            cv[i] = c;
            mt[i] = me;
            if (me != 0) anyMeth = 1;
            if (c != me) anyUnmeth = 1;
        }
    }
    /* constant methylation over the region */
    if (!anyMeth || !anyUnmeth) {
        out->constant = 1;
        return 1;
    }

    /* median coverage of each locus, including samples without coverage */
    for (i = 0; i < L; i++) {
        for (s = 0; s < N; s++) ws->medTmp[s] = ws->cv[i + (size_t) s * L];
        sort_small(ws->medTmp, N);
        ws->medCov[i] = (N % 2) ? ws->medTmp[N / 2] : (ws->medTmp[N / 2 - 1] + ws->medTmp[N / 2]) / 2;
    }

    /* rows with coverage */
    for (i = 0; i < L; i++) ws->locusMap[i] = -1;
    for (s = 0; s < N; s++) {
        int prevLocus = -1;
        for (i = 0; i < L; i++) {
            const int c = ws->cv[i + (size_t) s * L];
            if (c > 0) {
                const double me = ws->mt[i + (size_t) s * L];
                ws->rowSample[n] = s;
                ws->rowLocus[n] = i;
                ws->rowFirst[n] = prevLocus < 0;
                ws->rowGap[n] = prevLocus < 0 ? 0.0 : (double) ds->pos[offset + ix[i]] - ds->pos[offset + ix[prevLocus]];
                ws->locusMap[i] = 1;
                ws->z[n] = asin(2 * (me + 0.05) / (c + 0.1) - 1);
                prevLocus = i;
                n++;
            }
        }
    }
    m.ws = ws;
    m.n = n;
    m.p = p;
    m.nSamples = N;
    m.balanced = (n == L * N);
    m.rss = m.cQc = 0;

    /* Designs that are not of full rank on the rows with coverage are skipped. With
     * other covariates the factor of the loci replaces the position, and the attempts
     * below fail if the model is not of full rank. */
    if (!m.balanced && p == P->nTest) {
        int li;
        double *posRow = ws->y0;
        for (t = 0; t < n; t++) {
            li = ws->rowLocus[t];
            posRow[t] = ds->pos[offset + ix[li]];
        }
        if (!full_rank_rows(&m, X, P->nTest, posRow)) return 1;
    }

    /* compress the loci that have rows */
    for (i = 0; i < L; i++)
        if (ws->locusMap[i] == 1) ws->locusMap[i] = nPresentLoci++;
    m.nLoci = nPresentLoci;
    m.dfResid = n - nPresentLoci - p;
    if (m.dfResid <= 0) return 1;

    /* jitter for numerically constant values within a sample: a value shared by
     * at least L - 1 rows of a sample */
    for (s = 0, t = 0; s < N && !jitter; s++) {
        int from = t, to = t, cand = -1, count = 0, k;
        while (to < n && ws->rowSample[to] == s) to++;
        t = to;
        if (to - from < L - 1) continue;
        for (k = from; k < to; k++) {
            if (count == 0) { cand = k; count = 1; }
            else if (fabs(ws->z[k] - ws->z[cand]) <= 1e-14) count++;
            else count--;
        }
        if (cand >= 0) {
            count = 0;
            for (k = from; k < to; k++) if (fabs(ws->z[k] - ws->z[cand]) <= 1e-14) count++;
            if (count >= L - 1) jitter = 1;
        }
    }
    if (jitter)
        for (t = 0; t < n; t++) ws->z[t] += -0.001 + 0.002 * rng_unif(rng);

    /* coverage outliers in short regions: no variance weights */
    if (L < 40 && L <= 10) {
        double first[10], tmp[10];
        const int nf = n < L ? n : L;
        for (t = 0; t < nf; t++) first[t] = ws->medCov[ws->rowLocus[t]];
        grubbs = grubbs_outliers(first, nf, P, tmp);
    }
    for (t = 0; t < n; t++) if (!(ws->medCov[ws->rowLocus[t]] > 0)) medOK = 0;

    /* Map the loci, and bring the balanced z to L x N order (it already is) */
    for (t = 0; t < n; t++) ws->rowLocus[t] = ws->locusMap[ws->rowLocus[t]];

    {
        /* weights, then no weights, then neither weights nor correlation */
        const int car1 = L < 40;
        int attempt;
        for (attempt = 0; attempt < 3; attempt++) {
            const int useW = attempt == 0 && weighted && !grubbs && medOK;
            const int search = attempt < 2;
            double stat, phi;
            if (attempt == 1 && !(weighted && !grubbs && medOK)) continue;
            if (m.balanced) {
                int li;
                for (li = 0; li < L; li++) {
                    ws->c[li] = useW ? sqrt(ws->medCov[li]) : 1.0;
                    ws->gapL[li] = li == 0 ? 0.0 : (double) ds->pos[offset + ix[li]] - ds->pos[offset + ix[li - 1]];
                }
                if (P->nTest > 1)
                    for (t = 0; t < n; t++) ws->w[t] = ws->c[ws->rowLocus[t]];
            } else {
                int li = 0, raw;
                /* weights of the rows from the median coverage of the original loci */
                double *mapped = ws->c;
                for (raw = 0; raw < L; raw++)
                    if (ws->locusMap[raw] >= 0) mapped[li++] = ws->medCov[raw];
                for (t = 0; t < n; t++) ws->w[t] = useW ? sqrt(mapped[ws->rowLocus[t]]) : 1.0;
            }
            if (fit_attempt(&m, X, car1, P, search, &phi, out->beta, &stat)) {
                out->stat = stat;
                out->phi = phi;
                return 1;
            }
        }
    }
    return 1;
}
