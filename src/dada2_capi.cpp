/*
 * dada2_capi.cpp - Pure C API wrapper for DADA2 core algorithm.
 * Replaces Rmain.cpp for standalone (non-R) builds.
 * Compiled only with -DNO_RCPP.
 */
#ifndef NO_RCPP
/* This file is only for standalone builds. When compiled as R package,
 * Rmain.cpp provides the entry point instead. */
#else

#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <chrono>
#include "dada.h"
#include "dada2_capi.h"

/* Forward declaration of the standalone run_dada */
static B *run_dada_c(Raw **raws, int nraw, double *err_mat, int ncol_err,
                     int match, int mismatch, int gap_pen, int homo_gap_pen,
                     bool use_kmers, double kdist_cutoff, int band_size,
                     double omegaA, double omegaP, bool detect_singletons,
                     int max_clust, double min_fold, int min_hamming, int min_abund,
                     bool use_quals, bool vectorized_alignment,
                     bool multithread, bool verbose, int SSE, bool gapless, bool greedy);

/* Build transition matrix from clustering (mirrors b_make_transition_by_quality_matrix) */
static void fill_trans_matrix(B *b, Sub **subs, bool has_quals, int ncol_err,
                              int **out_trans, int *out_ncol) {
    int ncol = ncol_err;
    int *mat = (int *)calloc(16 * ncol, sizeof(int));
    if (!mat) { *out_trans = NULL; *out_ncol = 0; return; }

    for (unsigned int i = 0; i < b->nclust; i++) {
        for (unsigned int r = 0; r < b->bi[i]->nraw; r++) {
            Raw *raw = b->bi[i]->raw[r];
            if (!raw->correct) continue;
            Sub *sub = subs[raw->index];
            if (!sub) continue;

            /* For each position in the center (reference) sequence */
            Raw *center = b->bi[i]->center;
            for (unsigned int pos0 = 0; pos0 < center->length; pos0++) {
                uint16_t pos1 = sub->map[pos0];
                if (pos1 == GAP_GLYPH || pos1 >= raw->length) continue;

                int nti0 = ((int)b->bi[i]->center->seq[pos0]) - 1;
                int nti1 = ((int)raw->seq[pos1]) - 1;
                if (nti0 < 0 || nti0 > 3 || nti1 < 0 || nti1 > 3) continue;

                int qind = has_quals ? (int)raw->qual[pos1] : 0;
                if (qind >= ncol) qind = ncol - 1;

                int row = nti0 * 4 + nti1;
                mat[row * ncol + qind] += raw->reads;
            }
        }
    }
    *out_trans = mat;
    *out_ncol = ncol;
}

extern "C" DadaResult* dada2_run(
    const char **seqs,
    const int  *abundances,
    const int  *priors,
    int         nraw,
    const double *err_mat,
    int         ncol_err,
    const double *quals,
    int         maxlen,
    int match, int mismatch, int gap_pen,
    int use_kmers, double kdist_cutoff,
    int band_size,
    double omegaA, double omegaP, double omegaC,
    int detect_singletons,
    int max_clust,
    double min_fold, int min_hamming, int min_abund,
    int use_quals,
    int vectorized_alignment,
    int homo_gap_pen,
    int multithread,
    int verbose,
    int SSE, int gapless, int greedy)
{
    unsigned int i, r, index, pos;
    Raw *raw;
    auto t_total0 = std::chrono::steady_clock::now();

    if (nraw == 0) { fprintf(stderr, "dada2_run: zero input sequences.\n"); return NULL; }
    if (ncol_err < 1) { fprintf(stderr, "dada2_run: invalid error matrix.\n"); return NULL; }

    bool has_quals = (quals != NULL);

    /* Determine actual maxlen from sequences */
    unsigned int actual_maxlen = 0;
    for (index = 0; index < (unsigned)nraw; index++) {
        unsigned int slen = strlen(seqs[index]);
        if (slen > actual_maxlen) actual_maxlen = slen;
    }
    if (maxlen == 0) maxlen = actual_maxlen;

    /* Construct Raw structures */
    char seq_buf[SEQLEN];
    double qual_buf[SEQLEN];
    Raw **raws = (Raw **)malloc(nraw * sizeof(Raw *));
    if (!raws) { fprintf(stderr, "Memory allocation failed.\n"); return NULL; }

    for (index = 0; index < (unsigned)nraw; index++) {
        strncpy(seq_buf, seqs[index], SEQLEN - 1);
        seq_buf[SEQLEN - 1] = '\0';
        nt2int(seq_buf, seq_buf);
        if (has_quals) {
            for (pos = 0; pos < strlen(seqs[index]); pos++) {
                qual_buf[pos] = quals[(size_t)index * maxlen + pos];
            }
            raws[index] = raw_new(seq_buf, qual_buf, abundances[index], priors ? priors[index] : 0);
        } else {
            raws[index] = raw_new(seq_buf, NULL, abundances[index], priors ? priors[index] : 0);
        }
        raws[index]->index = index;
    }

    /* Build kmer structures */
    uint8_t *k8 = NULL;
    uint16_t *k16 = NULL;
    uint16_t *kord = NULL;
    if (use_kmers) {
        size_t n_kmer = 1 << (2 * KMER_SIZE);
        k8 = (uint8_t *)malloc(nraw * n_kmer * sizeof(uint8_t));
        k16 = (uint16_t *)malloc(nraw * n_kmer * sizeof(uint16_t));
        kord = (uint16_t *)malloc(nraw * actual_maxlen * sizeof(uint16_t));
        if (!k8 || !k16 || !kord) {
            free(k8); free(k16); free(kord);
            for (index = 0; index < (unsigned)nraw; index++) raw_free(raws[index]);
            free(raws);
            fprintf(stderr, "Memory allocation failed (kmer structures).\n");
            return NULL;
        }
        for (index = 0; index < (unsigned)nraw; index++) {
            raw = raws[index];
            raw->kmer8 = &k8[index * n_kmer];
            assign_kmer8(raw->kmer8, raw->seq, KMER_SIZE);
            raw->kmer = &k16[index * n_kmer];
            assign_kmer(raw->kmer, raw->seq, KMER_SIZE);
            raw->kord = &kord[index * actual_maxlen];
            assign_kmer_order(raw->kord, raw->seq, KMER_SIZE);
        }
    } else {
        for (index = 0; index < (unsigned)nraw; index++) {
            raws[index]->kmer8 = NULL;
            raws[index]->kmer = NULL;
            raws[index]->kord = NULL;
        }
    }

    /* Convert error matrix to mutable copy for internal use */
    double *err_mat_c = (double *)malloc(16 * ncol_err * sizeof(double));
    if (!err_mat_c) {
        for (index = 0; index < (unsigned)nraw; index++) raw_free(raws[index]);
        free(raws);
        if (use_kmers) { free(k8); free(k16); free(kord); }
        fprintf(stderr, "Memory allocation failed (error matrix).\n");
        return NULL;
    }
    memcpy(err_mat_c, err_mat, 16 * ncol_err * sizeof(double));

    /* Run DADA algorithm */
    auto t_core0 = std::chrono::steady_clock::now();
    B *bb = run_dada_c(raws, nraw, err_mat_c, ncol_err,
                       match, mismatch, gap_pen, homo_gap_pen,
                       use_kmers, kdist_cutoff, band_size,
                       omegaA, omegaP, detect_singletons,
                       max_clust, min_fold, min_hamming, min_abund,
                       use_quals, vectorized_alignment,
                       multithread, verbose, SSE, gapless, greedy);
    auto t_core1 = std::chrono::steady_clock::now();

    /* Build final alignments for output */
    auto t_subs0 = std::chrono::steady_clock::now();
    Sub **subs = (Sub **)malloc(bb->nraw * sizeof(Sub *));
    if (!subs) { fprintf(stderr, "Memory allocation failed.\n"); return NULL; }
    for (i = 0; i < bb->nclust; i++) {
        for (r = 0; r < bb->bi[i]->nraw; r++) {
            raw = bb->bi[i]->raw[r];
            subs[raw->index] = sub_new(bb->bi[i]->center, raw, match, mismatch, gap_pen, homo_gap_pen,
                                       false, 1.0, band_size, vectorized_alignment, SSE, gapless);
        }
    }
    auto t_subs1 = std::chrono::steady_clock::now();

    /* Assign final p-values */
    auto t_pval0 = std::chrono::steady_clock::now();
    for (i = 0; i < bb->nclust; i++) {
        for (r = 0; r < bb->bi[i]->nraw; r++) {
            raw = bb->bi[i]->raw[r];
            if (bb->bi[i]->center == raw) {
                raw->p = 1.0;
            } else {
                raw->p = calc_pA(raw->reads, raw->comp.lambda * bb->bi[i]->reads, true);
                if (raw->p < omegaC) { raw->correct = false; }
            }
        }
    }
    auto t_pval1 = std::chrono::steady_clock::now();

    /* Construct DadaResult */
    DadaResult *res = (DadaResult *)calloc(1, sizeof(DadaResult));
    if (!res) { fprintf(stderr, "Memory allocation failed.\n"); return NULL; }

    res->nclust = bb->nclust;
    res->nraw = bb->nraw;
    res->maxlen = actual_maxlen;

    /* Cluster sequences and abundances */
    res->cluster_seqs = (char **)malloc(bb->nclust * sizeof(char *));
    res->cluster_abunds = (int *)malloc(bb->nclust * sizeof(int));
    res->cluster_n0 = (int *)calloc(bb->nclust, sizeof(int));
    res->cluster_n1 = (int *)calloc(bb->nclust, sizeof(int));
    res->cluster_nunq = (int *)malloc(bb->nclust * sizeof(int));
    res->cluster_pval = (double *)malloc(bb->nclust * sizeof(double));

    char oseq[SEQLEN];
    for (i = 0; i < bb->nclust; i++) {
        /* Match R's clustering data.frame construction exactly:
         * representative sequence is the most abundant raw in the cluster,
         * but abundance/n0/n1/nunq count only raws that remain correct. */
        Raw *max_raw = NULL;
        unsigned int max_reads = 0;
        unsigned int max0 = 0, max1 = 0;
        int abund = 0;
        int nunq = 0;

        for (r = 0; r < bb->bi[i]->nraw; r++) {
            raw = bb->bi[i]->raw[r];
            if (raw->reads > max_reads) {
                max_raw = raw;
                max_reads = raw->reads;
            }
            if (!raw->correct) continue;

            abund += raw->reads;
            nunq++;
            if (subs[raw->index]) {
                if (subs[raw->index]->nsubs == 0) {
                    max0 += raw->reads;
                }
                if (subs[raw->index]->nsubs == 1) {
                    max1 += raw->reads;
                }
            }
        }

        if (!max_raw) {
            res->cluster_seqs[i] = strdup("");
        } else {
            int2nt(oseq, max_raw->seq);
            res->cluster_seqs[i] = strdup(oseq);
        }
        res->cluster_abunds[i] = abund;
        res->cluster_n0[i] = (int)max0;
        res->cluster_n1[i] = (int)max1;
        res->cluster_nunq[i] = nunq;
        res->cluster_pval[i] = bb->bi[i]->birth_pval;
    }

    /* Transition matrix */
    auto t_trans0 = std::chrono::steady_clock::now();
    fill_trans_matrix(bb, subs, has_quals, ncol_err, &res->trans, &res->ncol_trans);
    auto t_trans1 = std::chrono::steady_clock::now();

    /* Map: raw -> cluster (0-indexed, -1 for uncorrected) */
    res->map = (int *)malloc(bb->nraw * sizeof(int));
    res->pval = (double *)malloc(bb->nraw * sizeof(double));
    for (i = 0; i < bb->nclust; i++) {
        for (r = 0; r < bb->bi[i]->nraw; r++) {
            raw = bb->bi[i]->raw[r];
            res->map[raw->index] = raw->correct ? (int)i : -1;
            res->pval[raw->index] = raw->p;
        }
    }

    /* Cleanup */
    for (index = 0; index < bb->nraw; index++) sub_free(subs[index]);
    free(subs);
    b_free(bb);
    for (index = 0; index < (unsigned)nraw; index++) raw_free(raws[index]);
    free(raws);
    if (use_kmers) { free(k8); free(k16); free(kord); }
    free(err_mat_c);

    if (verbose) {
        auto t_total1 = std::chrono::steady_clock::now();
        double core_ms = std::chrono::duration<double, std::milli>(t_core1 - t_core0).count();
        double subs_ms = std::chrono::duration<double, std::milli>(t_subs1 - t_subs0).count();
        double pval_ms = std::chrono::duration<double, std::milli>(t_pval1 - t_pval0).count();
        double trans_ms = std::chrono::duration<double, std::milli>(t_trans1 - t_trans0).count();
        double total_ms = std::chrono::duration<double, std::milli>(t_total1 - t_total0).count();
        printf("CPU PROFILE: core=%.3f ms final_subs=%.3f ms final_pvals=%.3f ms trans=%.3f ms total=%.3f ms\n",
               core_ms, subs_ms, pval_ms, trans_ms, total_ms);
    }

    return res;
}

extern "C" void dada2_result_free(DadaResult *res) {
    if (!res) return;
    if (res->cluster_seqs) {
        for (int i = 0; i < res->nclust; i++) free(res->cluster_seqs[i]);
        free(res->cluster_seqs);
    }
    free(res->cluster_abunds);
    free(res->cluster_n0);
    free(res->cluster_n1);
    free(res->cluster_nunq);
    free(res->cluster_pval);
    free(res->trans);
    free(res->map);
    free(res->pval);
    free(res);
}

/* Standalone run_dada (mirrors Rmain.cpp::run_dada but takes flat arrays) */
static B *run_dada_c(Raw **raws, int nraw, double *err_mat, int ncol_err,
                     int match, int mismatch, int gap_pen, int homo_gap_pen,
                     bool use_kmers, double kdist_cutoff, int band_size,
                     double omegaA, double omegaP, bool detect_singletons,
                     int max_clust, double min_fold, int min_hamming, int min_abund,
                     bool use_quals, bool vectorized_alignment,
                     bool multithread, bool verbose, int SSE, bool gapless, bool greedy) {
    int newi = 0, nshuffle = 0;
    bool shuffled = false;
    auto t_compare = std::chrono::nanoseconds::zero();
    auto t_pupdate = std::chrono::nanoseconds::zero();
    auto t_shuffle = std::chrono::nanoseconds::zero();
    int n_compare_calls = 0;
    int n_pupdate_calls = 0;

    B *bb = b_new(raws, nraw, omegaA, omegaP, use_quals);

    /* Initial comparison - all raws vs cluster 0, no kmer screen */
    auto tc0 = std::chrono::steady_clock::now();
    b_compare_omp(bb, 0, err_mat, ncol_err, match, mismatch, gap_pen, homo_gap_pen,
                  use_kmers, 1.0, band_size, vectorized_alignment, SSE, gapless, greedy, verbose);
    auto tc1 = std::chrono::steady_clock::now();
    t_compare += (tc1 - tc0);
    n_compare_calls++;

    auto tp0 = std::chrono::steady_clock::now();
    b_p_update(bb, greedy, detect_singletons);
    auto tp1 = std::chrono::steady_clock::now();
    t_pupdate += (tp1 - tp0);
    n_pupdate_calls++;

    if (max_clust < 1) max_clust = bb->nraw;

    while ((bb->nclust < (unsigned)max_clust) && (newi = b_bud(bb, min_fold, min_hamming, min_abund, verbose))) {
        if (verbose) printf("\nNew Cluster C%d:", newi);

        auto tci0 = std::chrono::steady_clock::now();
        b_compare_omp(bb, newi, err_mat, ncol_err, match, mismatch, gap_pen, homo_gap_pen,
                      use_kmers, kdist_cutoff, band_size, vectorized_alignment, SSE, gapless, greedy, verbose);
        auto tci1 = std::chrono::steady_clock::now();
        t_compare += (tci1 - tci0);
        n_compare_calls++;

        nshuffle = 0;
        auto ts0 = std::chrono::steady_clock::now();
        do {
            shuffled = b_shuffle2(bb);
            if (verbose) printf("S");
        } while (shuffled && ++nshuffle < MAX_SHUFFLE);
        auto ts1 = std::chrono::steady_clock::now();
        t_shuffle += (ts1 - ts0);
        if (verbose && nshuffle >= MAX_SHUFFLE) printf("Warning: Reached maximum (%d) shuffles.\n", MAX_SHUFFLE);

        auto tpi0 = std::chrono::steady_clock::now();
        b_p_update(bb, greedy, detect_singletons);
        auto tpi1 = std::chrono::steady_clock::now();
        t_pupdate += (tpi1 - tpi0);
        n_pupdate_calls++;
    }

    if (verbose) {
        double compare_ms = std::chrono::duration<double, std::milli>(t_compare).count();
        double pupdate_ms = std::chrono::duration<double, std::milli>(t_pupdate).count();
        double shuffle_ms = std::chrono::duration<double, std::milli>(t_shuffle).count();
        printf("\nALIGN: %d aligns, %d shrouded (%d raw).\n", bb->nalign, bb->nshroud, bb->nraw);
        printf("CPU CORE PROFILE: compares=%d compare=%.3f ms p_updates=%d p_update=%.3f ms shuffle=%.3f ms clusters=%u\n",
               n_compare_calls, compare_ms, n_pupdate_calls, pupdate_ms, shuffle_ms, bb->nclust);
    }

    return bb;
}
#endif /* NO_RCPP */
