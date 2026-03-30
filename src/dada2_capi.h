#ifndef DADA2_CAPI_H
#define DADA2_CAPI_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int      nclust;
    int      nraw;
    int      maxlen;

    /* Per-cluster arrays (nclust entries) */
    char   **cluster_seqs;      /* ACGT sequences, null-terminated */
    int     *cluster_abunds;
    int     *cluster_n0;
    int     *cluster_n1;
    int     *cluster_nunq;
    double  *cluster_pval;

    /* Transition matrix: 16 x ncol_trans, row-major */
    int      ncol_trans;
    int     *trans;

    /* Map: nraw entries, 0-indexed cluster assignment (-1 for uncorrected) */
    int     *map;

    /* Per-raw p-values */
    double  *pval;
} DadaResult;

DadaResult* dada2_run(
    const char **seqs,
    const int  *abundances,
    const int  *priors,
    int         nraw,
    const double *err_mat,      /* 16 x ncol_err, row-major */
    int         ncol_err,
    const double *quals,        /* nraw x maxlen, row-major (avg Q per pos), or NULL */
    int         maxlen,
    /* alignment params */
    int match, int mismatch, int gap_pen,
    int use_kmers, double kdist_cutoff,
    int band_size,
    /* statistical params */
    double omegaA, double omegaP, double omegaC,
    int detect_singletons,
    int max_clust,
    double min_fold, int min_hamming, int min_abund,
    int use_quals,
    int vectorized_alignment,
    int homo_gap_pen,
    int multithread,
    int verbose,
    int SSE, int gapless, int greedy
);

void dada2_result_free(DadaResult *res);

/* Taxonomy assignment result */
typedef struct {
    int      nseq;      /* number of query sequences */
    int      nlevel;    /* number of taxonomy levels */
    int     *rval;      /* best genus index per query (nseq), 1-indexed, 0=NA */
    int     *rboot;     /* bootstrap counts (nseq x nlevel), row-major */
} TaxResult;

/* Assign taxonomy using naive Bayesian kmer classifier.
 *
 * seqs:         query sequences (nseq)
 * nseq:         number of query sequences
 * refs:         reference sequences (nref)
 * nref:         number of references
 * ref_to_genus: 0-indexed genus ID per reference (nref)
 * genusmat:     genus-to-level assignment matrix (ngenus x nlevel), row-major
 * ngenus:       number of unique genera
 * nlevel:       number of taxonomy levels
 * verbose:      print progress
 *
 * Returns a TaxResult* that must be freed with dada2_tax_result_free().
 */
TaxResult* dada2_assign_taxonomy(
    const char **seqs,
    int nseq,
    const char **refs,
    int nref,
    const int *ref_to_genus,
    const int *genusmat,
    int ngenus,
    int nlevel,
    int verbose
);

void dada2_tax_result_free(TaxResult *res);

/* Paired-read merging functions */

/* NW ends-free alignment of two ACGT strings.
 * Returns 0 on success, -1 on error.
 * Caller must free *al1_out and *al2_out with dada2_free_string(). */
int dada2_nwalign(const char *s1, const char *s2,
                  int match, int mismatch, int gap_p, int band,
                  char **al1_out, char **al2_out);

/* Evaluate an alignment: count matches, mismatches, indels (skipping end gaps). */
void dada2_eval_pair(const char *al1, const char *al2,
                     int *out_match, int *out_mismatch, int *out_indel);

/* Build consensus from two aligned strings.
 * prefer=1: al1 wins mismatches; prefer=2: al2 wins.
 * Caller must free result with dada2_free_string(). */
char *dada2_pair_consensus(const char *al1, const char *al2,
                           int prefer, int trim_overhang);

/* Reverse complement an ACGT string.
 * Caller must free result with dada2_free_string(). */
char *dada2_rc(const char *seq);

/* Free a string returned by the paired-read functions. */
void dada2_free_string(char *s);

/* ---- Chimera detection ---- */

typedef struct {
    int  n_seqs;
    int *nflag;    /* per-ASV: number of samples flagging as chimeric */
    int *nsam;     /* per-ASV: number of samples where ASV is present */
} ChimeraResult;

/* Check if seq is a bimera of the given parent sequences.
 * Returns 1 if bimera, 0 if not. */
int dada2_is_bimera(
    const char *seq,
    const char **parents,
    int n_parents,
    int allow_one_off,
    int min_one_off_par_dist,
    int match, int mismatch, int gap_p, int max_shift
);

/* Table-level consensus chimera detection.
 * mat: count matrix, column-major (nrow x ncol), rows=samples, cols=ASVs.
 * seqs: ASV sequences (ncol entries).
 * Returns a ChimeraResult* that must be freed with dada2_chimera_result_free(). */
ChimeraResult* dada2_table_bimera(
    const int *mat,
    int nrow, int ncol,
    const char **seqs,
    double min_fold,
    int min_abund,
    int allow_one_off,
    int min_one_off_par_dist,
    int match, int mismatch, int gap_p, int max_shift
);

void dada2_chimera_result_free(ChimeraResult *res);

#ifdef __cplusplus
}
#endif

#endif /* DADA2_CAPI_H */
