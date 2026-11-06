/*
 * One chromosome under one labelling of the samples: the effect of each CpG, its
 * smoothing, the candidate regions, and the fit of each region. The steps follow
 * dmrseq (Korthauer et al., Biostatistics 2019, functions bumphunt, estim, getEstimate,
 * smoother, refineEdges, trimEdges and regionScanner).
 */
#include "dmrseqx.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------ random */

static uint64_t splitmix64(uint64_t *x)
{
    uint64_t z = (*x += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

void rng_seed(rng_t *r, uint64_t seed, uint64_t a, uint64_t b)
{
    uint64_t x = seed;
    x ^= splitmix64(&a);
    x ^= splitmix64(&b) << 1;
    r->s[0] = splitmix64(&x);
    r->s[1] = splitmix64(&x);
}

/* xorshift128+ */
double rng_unif(rng_t *r)
{
    uint64_t s1 = r->s[0];
    const uint64_t s0 = r->s[1];
    const uint64_t result = s0 + s1;
    r->s[0] = s0;
    s1 ^= s1 << 23;
    r->s[1] = s1 ^ s0 ^ (s1 >> 18) ^ (s0 >> 5);
    return (double) (result >> 11) * (1.0 / 9007199254740992.0);
}

/* ---------------------------------------------------------------- workspace */

typedef struct { int a, b, sig; } candidate_t;

struct chrom_ws {
    double *rawBeta, *sdRaw, *weight, *fitted, *y, *x, *xd, *xFit, *yFit, *wFit, *tmp, *v1, *v0;
    unsigned char *smoothed;
    int *ind, *cluster, *g1, *g0, *level;
    candidate_t *cand;
    size_t capCand;
    fit_ws_t *fit;
    /* regression at each CpG */
    double *Q, *Qsub, *R, *Rinv, *Pfull, *Dsub, *ysub, *bhat, *Ulev, *betaFit;
    unsigned char *hasLevel;
};

chrom_ws_t *chrom_ws_new(int nSamples, int maxCpG, int nDesign, int nTest)
{
    chrom_ws_t *w = (chrom_ws_t *) calloc(1, sizeof(chrom_ws_t));
    const size_t n = (size_t) maxCpG + 1, N = (size_t) nSamples + 1, q = (size_t) nDesign;
    if (!w) return NULL;
    w->Q = (double *) malloc(N * q * sizeof(double));
    w->Qsub = (double *) malloc(N * q * sizeof(double));
    w->R = (double *) malloc(q * q * sizeof(double));
    w->Rinv = (double *) malloc(q * q * sizeof(double));
    w->Pfull = (double *) malloc(N * q * sizeof(double));
    w->Dsub = (double *) malloc(N * q * sizeof(double));
    w->ysub = (double *) malloc(N * sizeof(double));
    w->bhat = (double *) malloc(q * sizeof(double));
    w->Ulev = (double *) malloc(N * q * sizeof(double));
    w->betaFit = (double *) malloc((size_t) nTest * sizeof(double));
    w->level = (int *) malloc(N * sizeof(int));
    w->hasLevel = (unsigned char *) malloc(N);
    if (!w->Q || !w->Qsub || !w->R || !w->Rinv || !w->Pfull || !w->Dsub || !w->ysub || !w->bhat || !w->Ulev ||
        !w->betaFit || !w->level || !w->hasLevel) {
        chrom_ws_free(w);
        return NULL;
    }
    w->rawBeta = (double *) malloc(n * sizeof(double));
    w->sdRaw = (double *) malloc(n * sizeof(double));
    w->weight = (double *) malloc(n * sizeof(double));
    w->fitted = (double *) malloc(n * sizeof(double));
    w->y = (double *) malloc(n * sizeof(double));
    w->x = (double *) malloc(n * sizeof(double));
    w->xd = (double *) malloc(n * sizeof(double));
    w->xFit = (double *) malloc(n * sizeof(double));
    w->yFit = (double *) malloc(n * sizeof(double));
    w->wFit = (double *) malloc(n * sizeof(double));
    w->tmp = (double *) malloc(n * sizeof(double));
    w->v1 = (double *) malloc((size_t) (nSamples + 1) * sizeof(double));
    w->v0 = (double *) malloc((size_t) (nSamples + 1) * sizeof(double));
    w->smoothed = (unsigned char *) malloc(n);
    w->ind = (int *) malloc(n * sizeof(int));
    w->cluster = (int *) malloc(n * sizeof(int));
    w->g1 = (int *) malloc((size_t) (nSamples + 1) * sizeof(int));
    w->g0 = (int *) malloc((size_t) (nSamples + 1) * sizeof(int));
    w->fit = fit_ws_new();
    w->capCand = 1024;
    w->cand = (candidate_t *) malloc(w->capCand * sizeof(candidate_t));
    if (!w->rawBeta || !w->sdRaw || !w->weight || !w->fitted || !w->y || !w->x || !w->xd ||
        !w->xFit || !w->yFit || !w->wFit || !w->tmp || !w->v1 || !w->v0 || !w->smoothed ||
        !w->ind || !w->cluster || !w->g1 || !w->g0 || !w->fit || !w->cand) {
        chrom_ws_free(w);
        return NULL;
    }
    return w;
}

void chrom_ws_free(chrom_ws_t *w)
{
    if (!w) return;
    free(w->rawBeta); free(w->sdRaw); free(w->weight); free(w->fitted); free(w->y); free(w->x);
    free(w->xd); free(w->xFit); free(w->yFit); free(w->wFit); free(w->tmp); free(w->v1);
    free(w->v0); free(w->smoothed); free(w->ind); free(w->cluster); free(w->g1); free(w->g0);
    free(w->cand);
    free(w->Q); free(w->Qsub); free(w->R); free(w->Rinv); free(w->Pfull); free(w->Dsub); free(w->ysub);
    free(w->bhat); free(w->Ulev); free(w->betaFit); free(w->level); free(w->hasLevel);
    fit_ws_free(w->fit);
    free(w);
}

/* ------------------------------------------------------------ effect of a CpG */

static int cmp_double(const void *a, const void *b)
{
    const double x = *(const double *) a, y = *(const double *) b;
    return (x > y) - (x < y);
}

static void sort_values(double *x, int n)
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

static double median_sorted(const double *x, int n)
{
    return (n % 2) ? x[n / 2] : (x[n / 2 - 1] + x[n / 2]) / 2.0;
}

/* Median absolute deviation without constant (x is overwritten) */
static double mad_of(double *x, int n)
{
    double center;
    int i;
    sort_values(x, n);
    center = median_sorted(x, n);
    for (i = 0; i < n; i++) x[i] = fabs(x[i] - center);
    sort_values(x, n);
    return median_sorted(x, n);
}

/* The pooled difference in methylation, its robust variation between samples
 * and the smoothing weight of each CpG. The factor 1.4826 enters twice as in
 * dmrseq, once in the MAD and once in the standard deviation. */
static void estimate_effects(const dataset_t *ds, const params_t *P, size_t offset, int n,
                             const int *g1, int n1All, const int *g0, int n0All, chrom_ws_t *w)
{
    const int N = ds->nSamples;
    const size_t nTotal = ds->nTotal;
    const double scale = 2.0 / sqrt((double) N);
    int i, s;
    for (i = 0; i < n; i++) {
        double s1 = 0, t1 = 0, s0 = 0, t0 = 0, total = 0, sq1, sq0, sd, effect, covMean, lower;
        int n1 = 0, n0 = 0;
        const size_t row = offset + (size_t) i;
        for (s = 0; s < n1All; s++) {
            const size_t col = (size_t) g1[s] * nTotal + row;
            const int c = ds->Cov[col], m = ds->M[col];
            s1 += m;
            t1 += c;
            if (c > 0) w->v1[n1++] = (double) m / c;
        }
        for (s = 0; s < n0All; s++) {
            const size_t col = (size_t) g0[s] * nTotal + row;
            const int c = ds->Cov[col], m = ds->M[col];
            s0 += m;
            t0 += c;
            if (c > 0) w->v0[n0++] = (double) m / c;
        }
        total = t1 + t0;
        effect = s1 / t1 - s0 / t0;
        if (n1 == 1 && n0 == 1) {
            sq1 = sq0 = 1;
        } else if (n1 == 0 || n0 == 0) {
            sq1 = sq0 = NAN;
        } else {
            double mad;
            mad = 1.4826 * mad_of(w->v1, n1);
            sq1 = mad * mad;
            mad = 1.4826 * mad_of(w->v0, n0);
            sq0 = mad * mad;
            if (n1 == 1) sq1 = sq0;
            if (n0 == 1) sq0 = sq1;
        }
        sd = 1.4826 * sqrt(sq1 + sq0);
        if (sd < 1e-05) sd = 1e-05;
        covMean = total / N;
        lower = 1.0 / (covMean > 5 ? covMean : 5);
        w->rawBeta[i] = effect;
        w->sdRaw[i] = sd;
        w->weight[i] = isnan(sd) ? NAN : (covMean < P->covQ75 ? covMean : P->covQ75) / (sd > lower ? sd : lower);
        w->y[i] = effect / (sd * scale);
    }
}

/* ------------------------------------------- effect from a regression at each CpG */

/* QR by classical Gram-Schmidt with reorthogonalization: D, Q are m x q and R is upper
 * triangular q x q, all column-major. Returns 0 if not of full rank. */
static int qr_small(const double *D, int m, int q, double *Q, double *R)
{
    int i, j, k, pass;
    for (j = 0; j < q; j++) {
        double *v = Q + (size_t) j * m;
        double norm0 = 0, norm = 0;
        for (k = 0; k < m; k++) {
            v[k] = D[k + (size_t) j * m];
            norm0 += v[k] * v[k];
        }
        for (i = 0; i < q; i++) R[i + j * q] = 0;
        for (pass = 0; pass < 2; pass++) {
            for (i = 0; i < j; i++) {
                const double *u = Q + (size_t) i * m;
                double r = 0;
                for (k = 0; k < m; k++) r += u[k] * v[k];
                for (k = 0; k < m; k++) v[k] -= r * u[k];
                R[i + j * q] += r;
            }
        }
        for (k = 0; k < m; k++) norm += v[k] * v[k];
        if (!(norm > 1e-14 * norm0) || !(norm > 0)) return 0;
        norm = sqrt(norm);
        R[j + j * q] = norm;
        for (k = 0; k < m; k++) v[k] /= norm;
    }
    return 1;
}

static void invert_upper(const double *R, int q, double *Rinv)
{
    int i, j, k;
    for (j = 0; j < q; j++) {
        for (i = j + 1; i < q; i++) Rinv[i + j * q] = 0;
        for (i = j; i >= 0; i--) {
            double v = (i == j) ? 1.0 : 0.0;
            for (k = i + 1; k <= j; k++) v -= R[i + k * q] * Rinv[k + j * q];
            Rinv[i + j * q] = v / R[i + i * q];
        }
    }
}

/* (R'R)^-1 at the covariate of interest, its mean over several columns */
static double variance_factor(const double *Rinv, int q, int nTest)
{
    double total = 0;
    int c, k;
    for (c = 1; c <= nTest; c++)
        for (k = c; k < q; k++) total += Rinv[c + k * q] * Rinv[c + k * q];
    return total / nTest;
}

static double expit(double v)
{
    const double e = exp(v);
    return e / (1 + e);
}

/* The largest and smallest proportion that are not 0 and 1 give the value that replaces
 * 0 and 1 before the logit */
static double logit_offset(const dataset_t *ds, size_t offset, int n)
{
    double minNonZero = HUGE_VAL, maxNonOne = -HUGE_VAL, eps;
    int s, i;
    for (s = 0; s < ds->nSamples; s++) {
        const int *Ms = ds->M + (size_t) s * ds->nTotal + offset;
        const int *Cs = ds->Cov + (size_t) s * ds->nTotal + offset;
        for (i = 0; i < n; i++) {
            if (Cs[i] > 0) {
                const double p = (double) Ms[i] / Cs[i];
                if (p != 0 && p < minNonZero) minNonZero = p;
                if (p != 1 && p > maxNonOne) maxNonOne = p;
            }
        }
    }
    if (minNonZero < HUGE_VAL && maxNonOne > -HUGE_VAL)
        eps = minNonZero < 1 - maxNonOne ? minNonZero : 1 - maxNonOne;
    else if (minNonZero < HUGE_VAL)
        eps = minNonZero;
    else
        eps = 1 - maxNonOne;
    return eps == 1 ? 0.1 : eps;
}

/* The levels of a factor of interest, from the rows of the design */
static int find_levels(const double *design, int N, const params_t *P, int *level, double *Ulev)
{
    const int nCols = P->nTest + 1;
    int s, l, j, nLev = 0;
    for (s = 0; s < N; s++) {
        for (l = 0; l < nLev; l++) {
            for (j = 0; j < nCols; j++)
                if (design[s + j * N] != Ulev[l + j * N]) break;
            if (j == nCols) break;
        }
        if (l == nLev) {
            for (j = 0; j < nCols; j++) Ulev[l + j * N] = design[s + j * N];
            nLev++;
        }
        level[s] = l;
    }
    return nLev;
}

/* A regression of the logit of the methylation on the design at each CpG, by the rows
 * with coverage. The effect of a continuous covariate is the difference between its
 * values 1 and 0 on the scale of proportions, and for several groups the largest
 * difference between two groups. */
static void estimate_regression(const dataset_t *ds, const params_t *P, size_t offset, int n,
                                const double *design, int nLev, chrom_ws_t *w)
{
    const int N = ds->nSamples, q = P->nDesign, nTest = P->nTest, df = N - q;
    const size_t nTotal = ds->nTotal;
    const double eps = logit_offset(ds, offset, n);
    const int fullRank = qr_small(design, N, q, w->Q, w->R);
    double vbFull = 0;
    double *y = w->v1;
    int i, s, j, k;

    if (fullRank) {
        invert_upper(w->R, q, w->Rinv);
        vbFull = variance_factor(w->Rinv, q, nTest);
        for (s = 0; s < N; s++)
            for (j = 0; j < q; j++) {
                double v = 0;
                for (k = j; k < q; k++) v += w->Rinv[j + k * q] * w->Q[s + (size_t) k * N];
                w->Pfull[j + (size_t) s * q] = v;
            }
    }
    for (i = 0; i < n; i++) {
        const size_t row = offset + (size_t) i;
        double total = 0, rss = 0, vb = NAN, effect = NAN, sd, covMean, lower;
        int present = 0, ok = 0;
        double *b = w->bhat;
        for (s = 0; s < N; s++) {
            const size_t col = (size_t) s * nTotal + row;
            const int c = ds->Cov[col];
            total += c;
            if (c > 0) {
                double p = (double) ds->M[col] / c;
                if (p == 0) p = eps; else if (p == 1) p = 1 - eps;
                y[s] = log(p / (1 - p));
                present++;
            } else {
                y[s] = NAN;
            }
        }
        if (present == N) {
            if (fullRank) {
                for (j = 0; j < q; j++) {
                    const double *Pj = w->Pfull;
                    double v = 0;
                    for (s = 0; s < N; s++) v += Pj[j + (size_t) s * q] * y[s];
                    b[j] = v;
                }
                vb = vbFull;
                ok = 1;
            }
        } else if (present >= q) {
            /* The rows with coverage */
            int m = 0;
            for (s = 0; s < N; s++) {
                if (!isnan(y[s])) {
                    for (j = 0; j < q; j++) w->Dsub[m + (size_t) j * present] = design[s + j * N];
                    w->ysub[m] = y[s];
                    m++;
                }
            }
            if (qr_small(w->Dsub, present, q, w->Qsub, w->R)) {
                const double *Qt = w->Qsub;
                invert_upper(w->R, q, w->Rinv);
                for (j = 0; j < q; j++) {
                    double v = 0;
                    for (k = j; k < q; k++) {
                        double qy = 0;
                        for (m = 0; m < present; m++) qy += Qt[m + (size_t) k * present] * w->ysub[m];
                        v += w->Rinv[j + k * q] * qy;
                    }
                    b[j] = v;
                }
                vb = variance_factor(w->Rinv, q, nTest);
                ok = 1;
            }
        }
        if (ok) {
            for (s = 0; s < N; s++) {
                if (!isnan(y[s])) {
                    double fit = 0, res;
                    for (j = 0; j < q; j++) fit += design[s + j * N] * b[j];
                    res = y[s] - fit;
                    rss += res * res;
                }
            }
            if (nTest == 1) {
                /* without other covariates, at the values 1 and 0 of the covariate */
                effect = expit(b[0] + b[1]) - expit(b[0]);
            } else {
                double lo = HUGE_VAL, hi = -HUGE_VAL;
                int l;
                for (l = 0; l < nLev; l++) {
                    double g = 0, e;
                    for (j = 0; j <= nTest; j++) g += w->Ulev[l + j * N] * b[j];
                    e = expit(g);
                    if (e < lo) lo = e;
                    if (e > hi) hi = e;
                }
                effect = hi - lo;
            }
            sd = sqrt(vb * (rss / df));
            if (sd < 1e-05) sd = 1e-05;
        } else {
            sd = NAN;
        }
        covMean = total / N;
        lower = 1.0 / (covMean > 5 ? covMean : 5);
        w->rawBeta[i] = effect;
        w->sdRaw[i] = sd;
        w->weight[i] = isnan(sd) ? NAN : (covMean < P->covQ75 ? covMean : P->covQ75) / (sd > lower ? sd : lower);
        w->y[i] = effect / sd;
    }
}

/* ---------------------------------------------------------------- smoothing */

static int smooth_chromosome(const params_t *P, const int *pos, int n, chrom_ws_t *w)
{
    int a = 0, i;
    for (i = 0; i < n; i++) w->xd[i] = pos[i];
    while (a < n) {
        int b = a + 1;
        while (b < n && (double) pos[b] - pos[b - 1] <= P->maxGapSmooth) b++;
        {
            const int size = b - a;
            int nFit = 0;
            for (i = a; i < b; i++)
                if (!isnan(w->rawBeta[i]) && !isnan(w->weight[i])) nFit++;
            if (size >= P->minNumRegion && nFit >= P->minNumRegion) {
                double nn = P->bpSpan2 / ((double) pos[b - 1] - pos[a] + 1);
                if (P->minInSpan / size < nn) nn = P->minInSpan / size;
                if (0.75 < nn) nn = 0.75;
                if (nFit == size) {
                    smooth_block(w->xd + a, w->rawBeta + a, w->weight + a, nFit, nn, P->bpSpan, 0.8,
                                 w->xd + a, size, w->fitted + a);
                } else {
                    int k = 0;
                    for (i = a; i < b; i++) {
                        if (!isnan(w->rawBeta[i]) && !isnan(w->weight[i])) {
                            w->xFit[k] = w->xd[i];
                            w->yFit[k] = w->rawBeta[i];
                            w->wFit[k] = w->weight[i];
                            k++;
                        }
                    }
                    smooth_block(w->xFit, w->yFit, w->wFit, nFit, nn, P->bpSpan, 0.8,
                                 w->xd + a, size, w->fitted + a);
                }
                memset(w->smoothed + a, 1, (size_t) size);
            } else {
                for (i = a; i < b; i++) w->fitted[i] = NAN;
                memset(w->smoothed + a, 0, (size_t) size);
            }
        }
        a = b;
    }
    return 1;
}

/* --------------------------------------------------------------- candidates */

#define GE_PRECISION 1.490116119384765625e-08  /* sqrt(DBL_EPSILON) */

static int direction_of(double v, double cutoff)
{
    int d = (v >= cutoff || fabs(v - cutoff) <= GE_PRECISION) ? 1 : 0;
    if (v <= -cutoff) d = -1;
    return d;
}

/* Is the slope of x[first..last] on the index (0-based, inclusive) significant at
 * p < 0.01 and of the given sign? */
static int slope_significant(const double *x, int first, int last, int sign, const params_t *P)
{
    const int n = last - first + 1, df = n - 2;
    long double meanI = 0, meanX = 0, sxx = 0, sxy = 0, rss = 0;
    double slope, intercept, se, tval, crit;
    int i;
    if (n <= 4) return 0;
    for (i = first; i <= last; i++) {
        meanI += i;
        meanX += x[i];
    }
    meanI /= n;
    meanX /= n;
    for (i = first; i <= last; i++) {
        sxx += (i - meanI) * (i - meanI);
        sxy += (i - meanI) * (x[i] - meanX);
    }
    slope = (double) (sxy / sxx);
    intercept = (double) (meanX - slope * meanI);
    for (i = first; i <= last; i++) {
        const long double res = x[i] - (intercept + slope * i);
        rss += res * res;
    }
    se = sqrt((double) (rss / df) / (double) sxx);
    if (!(slope * sign > 0)) return 0;
    tval = fabs(slope) / se;
    if (isnan(tval)) return 0;
    if (df <= P->nTCrit) {
        crit = P->tCrit[df - 1];
    } else {
        /* Cornish-Fisher expansion of qt(0.995, df) */
        const double z = 2.5758293035489008, z3 = z * z * z, z5 = z3 * z * z;
        crit = z + (z3 + z) / (4.0 * df) + (5 * z5 + 16 * z3 + 3 * z) / (96.0 * df * (double) df);
    }
    return tval > crit;
}

/* For larger regions, remove the portions at the beginning and the end where the
 * effect differs drastically from the rest. x has the sign of the region removed,
 * *first and *last (0-based, inclusive) are the part kept. */
static void trim_region(const double *x, int n, int minNumRegion, const params_t *P, int *first, int *last)
{
    int mid = 0, i, newFirst = 0, newLast = n - 1;
    double maximum = x[0], minimum = x[0], cut, sumD;
    long double sum = 0, sum2 = 0;
    *first = 0;
    *last = n - 1;
    for (i = 0; i < n; i++) {
        if (x[i] > maximum) { maximum = x[i]; mid = i; }
        if (x[i] < minimum) minimum = x[i];
        sum += x[i];
    }
    if (!(maximum / minimum > 4.0 / 3.0)) return;
    /* mean() as in R: long double sum and one refinement */
    sum /= n;
    for (i = 0; i < n; i++) sum2 += (x[i] - sum);
    sumD = (double) (sum + sum2 / n);
    cut = (0.5 * (maximum - minimum) + minimum + 0.75 * sumD) / 2;
    /* 1-based mid below, as the rounding of 0.125 * n refers to it */
    if (slope_significant(x, 0, mid, 1, P)) {
        const int mid1 = mid + 1;
        double a = rint(mid1 - 0.125 * n), b = mid1 - 2, c = HUGE_VAL, v;
        if (a < 1) a = 1;
        if (b < 1) b = 1;
        for (i = 0; i <= mid; i++) if (x[i] >= cut) { c = i + 1; break; }
        v = a < b ? a : b;
        if (c < v) v = c;
        newFirst = (int) v - 1;
    }
    if (slope_significant(x, mid, n - 1, -1, P)) {
        const int mid1 = mid + 1;
        double a = rint(mid1 + 0.125 * n), b = mid1 + 2, c = -HUGE_VAL, v;
        if (a > n) a = n;
        if (b > n) b = n;
        for (i = n - 1; i >= mid; i--) if (x[i] >= cut) { c = i + 1; break; }
        v = a > b ? a : b;
        if (c > v) v = c;
        newLast = (int) v - 1;
    }
    if (newLast - newFirst + 1 >= minNumRegion) {
        *first = newFirst;
        *last = newLast;
    }
}

/* ------------------------------------------------------------------- output */

static int add_region(regions_t *out, const region_t *r, const double *beta, int nTest)
{
    if (out->n == out->cap) {
        const size_t cap = out->cap ? out->cap * 2 : 1024;
        region_t *rows = (region_t *) realloc(out->rows, cap * sizeof(region_t));
        double *betas;
        if (!rows) return 0;
        out->rows = rows;
        betas = (double *) realloc(out->beta, cap * (size_t) nTest * sizeof(double));
        if (!betas) return 0;
        out->beta = betas;
        out->cap = cap;
    }
    memcpy(out->beta + out->n * (size_t) nTest, beta, (size_t) nTest * sizeof(double));
    out->rows[out->n++] = *r;
    return 1;
}

/* ---------------------------------------------------------------- the task */

int process_chromosome(const dataset_t *ds, const params_t *P, size_t offset, int n,
                       const double *design, const double *X, int p, uint64_t seed, uint64_t pass,
                       uint64_t chrom, chrom_ws_t *w, regions_t *out)
{
    const int N = ds->nSamples;
    const size_t nTotal = ds->nTotal;
    const int *pos = ds->pos + offset;
    const int minNum = P->minNumRegion;
    int n1 = 0, n0 = 0, s, i, m = 0, nCand = 0, c, cl = 0, nLev = 0;
    rng_t rng;

    if (n < minNum) return 1;
    rng_seed(&rng, seed, pass, chrom);
    if (P->factor) nLev = find_levels(design, N, P, w->level, w->Ulev);
    if (P->twoGroup) {
        for (s = 0; s < N; s++) {
            if (design[s + N] == 1) w->g1[n1++] = s;
            else w->g0[n0++] = s;
        }
        estimate_effects(ds, P, offset, n, w->g1, n1, w->g0, n0, w);
    } else {
        estimate_regression(ds, P, offset, n, design, nLev, w);
    }
    if (P->smooth) {
        smooth_chromosome(P, pos, n, w);
        for (i = 0; i < n; i++) w->x[i] = w->smoothed[i] ? w->fitted[i] : w->y[i];
    } else {
        for (i = 0; i < n; i++) w->x[i] = w->rawBeta[i];
    }

    /* the CpGs with a value, and the clusters of all CpGs */
    for (i = 0; i < n; i++) {
        if (i > 0 && (double) pos[i] - pos[i - 1] > P->maxGap) cl++;
        w->cluster[i] = cl;
        if (!isnan(w->x[i])) w->ind[m++] = i;
    }
    if (m == 0) return 1;

    /* runs of the same direction within clusters */
    {
        int t = 0;
        while (t < m) {
            const int dir = direction_of(w->x[w->ind[t]], P->cutoff);
            const int clust = w->cluster[w->ind[t]];
            int e = t + 1;
            while (e < m && w->cluster[w->ind[e]] == clust && direction_of(w->x[w->ind[e]], P->cutoff) == dir) e++;
            if (dir != 0) {
                if ((size_t) nCand == w->capCand) {
                    candidate_t *cc = (candidate_t *) realloc(w->cand, w->capCand * 2 * sizeof(candidate_t));
                    if (!cc) return 0;
                    w->cand = cc;
                    w->capCand *= 2;
                }
                w->cand[nCand].a = t;
                w->cand[nCand].b = e - 1;
                w->cand[nCand].sig = dir;
                nCand++;
            }
            t = e;
        }
    }

    for (c = 0; c < nCand; c++) {
        candidate_t *cd = &w->cand[c];
        fit_out_t fo;
        int a = cd->a, b = cd->b, r, L;

        /* first refinement: from the first to the last CpG whose raw standardized
         * effect passes the cutoff in the direction of the region */
        if (b - a + 1 > minNum) {
            int lo = -1, hi = -1, t;
            for (t = a; t <= b; t++) {
                const double yv = w->y[w->ind[t]];
                if (!isnan(yv) && direction_of(yv, P->cutoff) == cd->sig) {
                    if (lo < 0) lo = t;
                    hi = t;
                }
            }
            if (lo >= 0 && hi - lo + 1 >= minNum) {
                a = lo;
                b = hi;
            }
        }
        /* second refinement */
        if (b - a + 1 > minNum) {
            const int len = b - a + 1;
            int t, f, l;
            for (t = 0; t < len; t++) w->tmp[t] = w->x[w->ind[a + t]] * cd->sig;
            trim_region(w->tmp, len, minNum, P, &f, &l);
            b = a + l;
            a = a + f;
        }
        L = b - a + 1;
        if (L < minNum) continue;

        /* for factors, at least two CpGs with coverage in all groups. As in dmrseq the
         * row numbers count the CpGs with a value. */
        if (P->factor) {
            int count = 0;
            for (r = a; r <= b && count < 2; r++) {
                const size_t row = offset + (size_t) r;
                int found = 0;
                memset(w->hasLevel, 0, (size_t) nLev);
                for (s = 0; s < N && found < nLev; s++) {
                    if (ds->Cov[(size_t) s * nTotal + row] > 0 && !w->hasLevel[w->level[s]]) {
                        w->hasLevel[w->level[s]] = 1;
                        found++;
                    }
                }
                if (found == nLev) count++;
            }
            if (count < 2) continue;
        }

        fo.beta = w->betaFit;
        if (!fit_region(ds, P, offset, w->ind + a, L, X, p, 1, &rng, w->fit, &fo)) return 0;
        if (fo.constant) continue;

        /* two groups, extreme coefficients: refit without variance weights */
        if (P->twoGroup && !isnan(fo.beta[0])) {
            long double total = 0;
            int count = 0, refit = 0;
            for (r = a; r <= b; r++) {
                const size_t row = offset + (size_t) r;
                double m1 = 0, m0 = 0;
                int c1 = 0, c0 = 0;
                for (s = 0; s < n1; s++) {
                    const size_t col = (size_t) w->g1[s] * nTotal + row;
                    if (ds->Cov[col] > 0) { m1 += (double) ds->M[col] / ds->Cov[col]; c1++; }
                }
                for (s = 0; s < n0; s++) {
                    const size_t col = (size_t) w->g0[s] * nTotal + row;
                    if (ds->Cov[col] > 0) { m0 += (double) ds->M[col] / ds->Cov[col]; c0++; }
                }
                if (c1 > 0 && c0 > 0) {
                    total += m1 / c1 - m0 / c0;
                    count++;
                }
            }
            if (fabs(fo.beta[0]) > M_PI) {
                refit = 1;
            } else if (count > 0) {
                const double diff = (double) (total / count);
                const int sd = (diff > 0) - (diff < 0), sb = (fo.beta[0] > 0) - (fo.beta[0] < 0);
                if (fabs(diff - fo.beta[0] / M_PI) > 1.0 / 3.0 && sd != sb) refit = 1;
            }
            if (refit && !fit_region(ds, P, offset, w->ind + a, L, X, p, 0, &rng, w->fit, &fo)) return 0;
        }

        {
            region_t reg;
            long double area = 0;
            int t;
            for (t = a; t <= b; t++) area += w->x[w->ind[t]];
            reg.start = pos[w->ind[a]];
            reg.end = pos[w->ind[b]];
            reg.indexStart = w->ind[a];
            reg.indexEnd = w->ind[b];
            reg.L = L;
            reg.area = fabs((double) area);
            reg.stat = fo.stat;
            if (!add_region(out, &reg, fo.beta, P->nTest)) return 0;
        }
    }
    return 1;
}
