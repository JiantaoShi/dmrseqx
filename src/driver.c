/*
 * Runs the tasks (one labelling of the samples on one chromosome) in worker threads
 * and returns the regions to R. The workers don't touch the R API; the main thread
 * waits for them, watches for user interrupts and reports progress.
 */
#include "dmrseqx.h"
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <pthread.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef struct {
    const dataset_t *ds;
    const params_t *P;
    int nPass, nChr, p, maxCpG;
    const double *chrStart;
    const double *design;
    const double *X;
    uint64_t seed;
    regions_t *results;
    int *passLeft;
    pthread_mutex_t lock;
    pthread_cond_t cond;
    int next, done, failed, stop, active;
} pool_t;

static void *worker(void *arg)
{
    pool_t *pool = (pool_t *) arg;
    const int N = pool->ds->nSamples, total = pool->nPass * pool->nChr;
    chrom_ws_t *ws = chrom_ws_new(N, pool->maxCpG, pool->P->nDesign, pool->P->nTest);

    if (!ws) {
        pthread_mutex_lock(&pool->lock);
        pool->failed = 1;
        pool->active--;
        pthread_cond_broadcast(&pool->cond);
        pthread_mutex_unlock(&pool->lock);
        return NULL;
    }
    for (;;) {
        int task, pass, chrom, ok;
        pthread_mutex_lock(&pool->lock);
        if (pool->stop || pool->failed || pool->next >= total) {
            pool->active--;
            pthread_cond_broadcast(&pool->cond);
            pthread_mutex_unlock(&pool->lock);
            break;
        }
        task = pool->next++;
        pthread_mutex_unlock(&pool->lock);

        pass = task / pool->nChr;
        chrom = task % pool->nChr;
        ok = process_chromosome(pool->ds, pool->P, (size_t) pool->chrStart[chrom],
                                (int) (pool->chrStart[chrom + 1] - pool->chrStart[chrom]),
                                pool->design + (size_t) pass * N * pool->P->nDesign,
                                pool->X + (size_t) pass * N * pool->p, pool->p, pool->seed,
                                (uint64_t) pass, (uint64_t) chrom, ws, &pool->results[task]);

        pthread_mutex_lock(&pool->lock);
        if (!ok) pool->failed = 1;
        pool->done++;
        pool->passLeft[pass]--;
        pthread_cond_broadcast(&pool->cond);
        pthread_mutex_unlock(&pool->lock);
    }
    chrom_ws_free(ws);
    return NULL;
}

static void check_interrupt_fn(void *dummy)
{
    (void) dummy;
    R_CheckUserInterrupt();
}

static int interrupted(void)
{
    return R_ToplevelExec(check_interrupt_fn, NULL) == FALSE;
}

static void read_params(SEXP params, SEXP tCrit, SEXP grubbsT, params_t *P)
{
    const double *v = REAL(params);
    int i;
    if (XLENGTH(params) != 14) error("Invalid parameters.");
    P->cutoff = v[0];
    P->bpSpan = v[1];
    P->minInSpan = v[2];
    P->bpSpan2 = v[3];
    P->maxGapSmooth = v[4];
    P->maxGap = v[5];
    P->covQ75 = v[6];
    P->minNumRegion = (int) v[7];
    P->smooth = (int) v[8];
    P->fitGlobal = (int) v[9];
    P->twoGroup = (int) v[10];
    P->factor = (int) v[11];
    P->nTest = (int) v[12];
    P->nDesign = (int) v[13];
    P->tCrit = REAL(tCrit);
    P->nTCrit = (int) XLENGTH(tCrit);
    if (XLENGTH(grubbsT) != 8) error("Invalid critical values.");
    for (i = 0; i < 8; i++) P->grubbsT[i + 3] = REAL(grubbsT)[i];
}

SEXP C_dmrseqx_run(SEXP M, SEXP Cov, SEXP pos, SEXP chrStart, SEXP design, SEXP X, SEXP params,
                   SEXP tCrit, SEXP grubbsT, SEXP threads, SEXP seed, SEXP verbose)
{
    const int *dimM = INTEGER(getAttrib(M, R_DimSymbol));
    const int nThreads = asInteger(threads), verb = asLogical(verbose);
    dataset_t ds;
    params_t P;
    pool_t pool;
    pthread_t *tid;
    int i, started = 0, nextReport = 0, total, wasInterrupted = 0;
    size_t nOut = 0;
    SEXP out, names, col[10];
    const char *colNames[] = {"pass", "chrom", "start", "end", "indexStart", "indexEnd", "L", "area", "stat", "beta"};

    if (TYPEOF(M) != INTSXP || TYPEOF(Cov) != INTSXP || TYPEOF(pos) != INTSXP)
        error("Counts and positions must be integers.");
    if (XLENGTH(pos) != dimM[0] || XLENGTH(M) != XLENGTH(Cov)) error("Inconsistent dimensions.");
    ds.M = INTEGER(M);
    ds.Cov = INTEGER(Cov);
    ds.pos = INTEGER(pos);
    ds.nTotal = (size_t) dimM[0];
    ds.nSamples = dimM[1];
    read_params(params, tCrit, grubbsT, &P);

    memset(&pool, 0, sizeof(pool));
    pool.ds = &ds;
    pool.P = &P;
    pool.nChr = (int) XLENGTH(chrStart) - 1;
    pool.nPass = (int) (XLENGTH(design) / ((R_xlen_t) ds.nSamples * P.nDesign));
    pool.p = (int) (XLENGTH(X) / ((R_xlen_t) ds.nSamples * pool.nPass));
    pool.chrStart = REAL(chrStart);
    pool.design = REAL(design);
    pool.X = REAL(X);
    pool.seed = (uint64_t) asReal(seed);
    if (pool.nPass < 1 || pool.p < P.nTest || P.nTest < 1 || P.nTest > 60) error("Inconsistent design.");
    for (i = 0; i < pool.nChr; i++) {
        const int n = (int) (pool.chrStart[i + 1] - pool.chrStart[i]);
        if (n > pool.maxCpG) pool.maxCpG = n;
    }
    total = pool.nPass * pool.nChr;
    pool.results = (regions_t *) R_alloc((size_t) total, sizeof(regions_t));
    memset(pool.results, 0, sizeof(regions_t) * (size_t) total);
    pool.passLeft = (int *) R_alloc((size_t) pool.nPass, sizeof(int));
    for (i = 0; i < pool.nPass; i++) pool.passLeft[i] = pool.nChr;
    tid = (pthread_t *) R_alloc((size_t) nThreads, sizeof(pthread_t));
    pthread_mutex_init(&pool.lock, NULL);
    pthread_cond_init(&pool.cond, NULL);

    pthread_mutex_lock(&pool.lock);
    pool.active = nThreads;
    pthread_mutex_unlock(&pool.lock);
    for (i = 0; i < nThreads; i++) {
        if (pthread_create(&tid[i], NULL, worker, &pool) != 0) {
            pthread_mutex_lock(&pool.lock);
            pool.active -= nThreads - i;
            if (i == 0) pool.failed = 1;
            pthread_mutex_unlock(&pool.lock);
            break;
        }
        started++;
    }

    pthread_mutex_lock(&pool.lock);
    while (pool.active > 0) {
        struct timespec ts;
        clock_gettime(CLOCK_REALTIME, &ts);
        ts.tv_nsec += 200000000L;
        if (ts.tv_nsec >= 1000000000L) {
            ts.tv_nsec -= 1000000000L;
            ts.tv_sec += 1;
        }
        pthread_cond_timedwait(&pool.cond, &pool.lock, &ts);
        pthread_mutex_unlock(&pool.lock);
        if (!wasInterrupted && interrupted()) {
            wasInterrupted = 1;
            pthread_mutex_lock(&pool.lock);
            pool.stop = 1;
            pthread_mutex_unlock(&pool.lock);
        }
        pthread_mutex_lock(&pool.lock);
        while (verb && nextReport < pool.nPass && pool.passLeft[nextReport] == 0) {
            size_t count = 0;
            int c;
            for (c = 0; c < pool.nChr; c++) count += pool.results[nextReport * pool.nChr + c].n;
            if (nextReport == 0)
                REprintf("* %lu candidates detected\n", (unsigned long) count);
            else
                REprintf("* %d out of %d permutations completed (%lu null candidates)\n",
                         nextReport, pool.nPass - 1, (unsigned long) count);
            nextReport++;
        }
    }
    pthread_mutex_unlock(&pool.lock);
    for (i = 0; i < started; i++) pthread_join(tid[i], NULL);
    pthread_mutex_destroy(&pool.lock);
    pthread_cond_destroy(&pool.cond);

    if (wasInterrupted || pool.failed) {
        for (i = 0; i < total; i++) {
            free(pool.results[i].rows);
            free(pool.results[i].beta);
        }
        if (wasInterrupted) error("Interrupted.");
        error("Out of memory.");
    }

    for (i = 0; i < total; i++) nOut += pool.results[i].n;
    PROTECT(out = allocVector(VECSXP, 10));
    PROTECT(names = allocVector(STRSXP, 10));
    for (i = 0; i < 10; i++) {
        SET_STRING_ELT(names, i, mkChar(colNames[i]));
        col[i] = allocVector(i < 7 ? INTSXP : REALSXP, (R_xlen_t) (i == 9 ? nOut * (size_t) P.nTest : nOut));
        SET_VECTOR_ELT(out, i, col[i]);
    }
    setAttrib(out, R_NamesSymbol, names);
    {
        size_t k = 0;
        for (i = 0; i < total; i++) {
            const int pass = i / pool.nChr, chrom = i % pool.nChr;
            const int first = (int) pool.chrStart[chrom];
            size_t j;
            int c;
            for (j = 0; j < pool.results[i].n; j++, k++) {
                const region_t *r = &pool.results[i].rows[j];
                INTEGER(col[0])[k] = pass;
                INTEGER(col[1])[k] = chrom + 1;
                INTEGER(col[2])[k] = r->start;
                INTEGER(col[3])[k] = r->end;
                INTEGER(col[4])[k] = first + r->indexStart + 1;
                INTEGER(col[5])[k] = first + r->indexEnd + 1;
                INTEGER(col[6])[k] = r->L;
                REAL(col[7])[k] = r->area;
                REAL(col[8])[k] = isnan(r->stat) ? NA_REAL : r->stat;
                for (c = 0; c < P.nTest; c++) {
                    const double b = pool.results[i].beta[j * (size_t) P.nTest + (size_t) c];
                    REAL(col[9])[k + nOut * (size_t) c] = isnan(b) ? NA_REAL : b;
                }
            }
            free(pool.results[i].rows);
            free(pool.results[i].beta);
        }
    }
    UNPROTECT(2);
    return out;
}

/* ------------------------------------------------- access for checks and tests */

/* One region: M and Cov are L x samples, X samples x p */
SEXP C_fit_one(SEXP M, SEXP Cov, SEXP pos, SEXP X, SEXP weighted, SEXP params, SEXP tCrit,
               SEXP grubbsT, SEXP seed)
{
    const int *dim = INTEGER(getAttrib(M, R_DimSymbol));
    const int L = dim[0];
    int *ix;
    dataset_t ds;
    params_t P;
    fit_ws_t *ws;
    fit_out_t fo;
    rng_t rng;
    int i, ok;
    SEXP out;

    if (TYPEOF(M) != INTSXP || TYPEOF(Cov) != INTSXP || TYPEOF(pos) != INTSXP)
        error("Counts and positions must be integers.");
    ds.M = INTEGER(M);
    ds.Cov = INTEGER(Cov);
    ds.pos = INTEGER(pos);
    ds.nTotal = (size_t) L;
    ds.nSamples = dim[1];
    read_params(params, tCrit, grubbsT, &P);
    ix = (int *) R_alloc((size_t) L, sizeof(int));
    for (i = 0; i < L; i++) ix[i] = i;
    rng_seed(&rng, (uint64_t) asReal(seed), 0, 0);
    fo.beta = (double *) R_alloc((size_t) P.nTest, sizeof(double));
    ws = fit_ws_new();
    if (!ws) error("Out of memory.");
    ok = fit_region(&ds, &P, 0, ix, L, REAL(X), (int) (XLENGTH(X) / dim[1]), asLogical(weighted), &rng, ws, &fo);
    fit_ws_free(ws);
    if (!ok) error("Out of memory.");
    PROTECT(out = allocVector(REALSXP, 3 + P.nTest));
    REAL(out)[0] = isnan(fo.stat) ? NA_REAL : fo.stat;
    REAL(out)[1] = isnan(fo.phi) ? NA_REAL : fo.phi;
    REAL(out)[2] = fo.constant;
    for (i = 0; i < P.nTest; i++) REAL(out)[3 + i] = isnan(fo.beta[i]) ? NA_REAL : fo.beta[i];
    UNPROTECT(1);
    return out;
}

/* One smoothing block */
SEXP C_smooth_one(SEXP xFit, SEXP yFit, SEXP wFit, SEXP nn, SEXP h, SEXP cut, SEXP xEval)
{
    SEXP out;
    PROTECT(out = allocVector(REALSXP, XLENGTH(xEval)));
    smooth_block(REAL(xFit), REAL(yFit), REAL(wFit), (long) XLENGTH(xFit), asReal(nn), asReal(h),
                 asReal(cut), REAL(xEval), (long) XLENGTH(xEval), REAL(out));
    UNPROTECT(1);
    return out;
}

static const R_CallMethodDef call_methods[] = {
    {"C_dmrseqx_run", (DL_FUNC) &C_dmrseqx_run, 12},
    {"C_fit_one", (DL_FUNC) &C_fit_one, 9},
    {"C_smooth_one", (DL_FUNC) &C_smooth_one, 7},
    {NULL, NULL, 0}
};

void R_init_dmrseqx(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, call_methods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
