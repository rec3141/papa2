/*
 * taxonomy_capi.cpp — Pure C API for dada2 taxonomy assignment
 *
 * Ports C_assign_taxonomy2 from taxonomy.cpp to the standalone C API,
 * removing Rcpp/RcppParallel dependencies. Uses OpenMP for parallelism.
 *
 * Algorithm (Wang et al. 2007, as implemented in dada2):
 *   1. Group references by genus, build per-genus kmer probability table
 *   2. For each query: find best genus via log-probability scoring
 *   3. Bootstrap (100x): subsample 1/8 of query kmers, find best genus
 *   4. Count bootstrap agreement at each taxonomic level
 */

#include "dada2_capi.h"
#include "dada.h"
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cfloat>
#include <cstdio>
#include <random>
#include <algorithm>

#define TAX_KMER_SIZE 8
#define TAX_N_KMERS (1 << (2 * TAX_KMER_SIZE))  /* 65536 */
#define TAX_NBOOT 100

/* Compute kmer integer from sequence position */
static int tax_kmer_at(const char *seq, unsigned int k) {
    unsigned int j;
    int kmer = 0;
    for (j = 0; j < k; j++) {
        int nti;
        switch (seq[j]) {
            case 'A': nti = 0; break;
            case 'C': nti = 1; break;
            case 'G': nti = 2; break;
            case 'T': nti = 3; break;
            default:  return -1;
        }
        kmer = 4 * kmer + nti;
    }
    return kmer;
}

/* Build boolean kmer presence vector for a sequence */
static void tax_kvec_build(const char *seq, unsigned int k, unsigned char *kvec) {
    size_t len = strlen(seq);
    size_t n_kmers = (size_t)(1 << (2 * k));
    memset(kvec, 0, n_kmers);
    if (len < k) return;
    size_t klen = len - k + 1;
    for (size_t i = 0; i < klen; i++) {
        int kmer = tax_kmer_at(&seq[i], k);
        if (kmer >= 0 && kmer < (int)n_kmers) {
            kvec[kmer] = 1;
        }
    }
}

/* Build sorted kmer array for a query sequence. Returns length written. */
static unsigned int tax_karray_build(const char *seq, unsigned int k, int *karray) {
    size_t len = strlen(seq);
    if (len < k) return 0;
    size_t klen = len - k + 1;
    unsigned int j = 0;
    for (size_t i = 0; i < klen; i++) {
        int kmer = tax_kmer_at(&seq[i], k);
        if (kmer >= 0) {
            karray[j++] = kmer;
        }
    }
    std::sort(karray, karray + j);
    return j;
}

/* Find the genus with highest log-probability for a given kmer array */
static int get_best_genus(int *karray, float *out_logp, unsigned int arraylen,
                          unsigned int n_kmers, unsigned int ngenus,
                          float *lgk_probability) {
    int max_g = -1;
    float max_logp = -FLT_MAX;
    unsigned int nmax = 0;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> cunif(0.0, 1.0);

    for (unsigned int g = 0; g < ngenus; g++) {
        float *lgk_v = &lgk_probability[g * n_kmers];
        float logp = 0.0f;

        for (unsigned int pos = 0; pos < arraylen; pos++) {
            logp += lgk_v[karray[pos]];
            if (logp < max_logp) break;  /* early termination */
        }

        if (max_logp < -FLT_MAX + 1 || logp > max_logp) {
            max_logp = logp;
            max_g = (int)g;
            nmax = 1;
        } else if (logp == max_logp) {
            nmax++;
            if (cunif(gen) < 1.0 / nmax) {
                max_g = (int)g;
            }
        }
    }
    *out_logp = max_logp;
    return max_g;
}

extern "C"
TaxResult* dada2_assign_taxonomy(
    const char **seqs, int nseq,
    const char **refs, int nref,
    const int *ref_to_genus,
    const int *genusmat,
    int ngenus, int nlevel,
    int verbose)
{
    unsigned int k = TAX_KMER_SIZE;
    size_t n_kmers = TAX_N_KMERS;

    /* Count references per genus */
    float *genus_num_plus1 = (float *)calloc(ngenus, sizeof(float));
    for (int i = 0; i < nref; i++) {
        genus_num_plus1[ref_to_genus[i]]++;
    }
    for (int g = 0; g < ngenus; g++) {
        genus_num_plus1[g]++;
    }

    /* Build per-genus kmer probability table */
    float *kmer_prior = (float *)calloc(n_kmers, sizeof(float));
    float *lgk_probability = (float *)calloc(ngenus * n_kmers, sizeof(float));
    unsigned char *ref_kv = (unsigned char *)malloc(n_kmers);

    if (!genus_num_plus1 || !kmer_prior || !lgk_probability || !ref_kv) {
        fprintf(stderr, "[ERROR] Memory allocation failed in taxonomy.\n");
        return NULL;
    }

    for (int i = 0; i < nref; i++) {
        if (verbose && i > 0 && i % 100000 == 0) {
            fprintf(stderr, "[INFO] Processed %d/%d references\n", i, nref);
        }
        tax_kvec_build(refs[i], k, ref_kv);
        int g = ref_to_genus[i];
        float *lgk_v = &lgk_probability[g * n_kmers];
        for (size_t km = 0; km < n_kmers; km++) {
            if (ref_kv[km]) {
                lgk_v[km]++;
                kmer_prior[km]++;
            }
        }
    }
    free(ref_kv);

    /* Compute kmer priors and log-probabilities */
    for (size_t km = 0; km < n_kmers; km++) {
        kmer_prior[km] = (kmer_prior[km] + 0.5f) / (1.0f + nref);
    }
    for (int g = 0; g < ngenus; g++) {
        float *lgk_v = &lgk_probability[g * n_kmers];
        for (size_t km = 0; km < n_kmers; km++) {
            lgk_v[km] = logf((lgk_v[km] + kmer_prior[km]) / genus_num_plus1[g]);
        }
    }
    free(kmer_prior);
    free(genus_num_plus1);

    if (verbose) {
        fprintf(stderr, "[INFO] Reference kmer table built (%d genera x %zu kmers = %.1f MB)\n",
                ngenus, n_kmers, (double)(ngenus * n_kmers * sizeof(float)) / 1e6);
    }

    /* Get max query length for array allocation */
    unsigned int max_arraylen = 0;
    for (int i = 0; i < nseq; i++) {
        unsigned int slen = strlen(seqs[i]);
        if (slen > k && (slen - k + 1) > max_arraylen) {
            max_arraylen = slen - k + 1;
        }
    }

    /* Allocate result */
    TaxResult *result = (TaxResult *)malloc(sizeof(TaxResult));
    result->nseq = nseq;
    result->nlevel = nlevel;
    result->rval = (int *)malloc(nseq * sizeof(int));
    result->rboot = (int *)calloc(nseq * nlevel, sizeof(int));

    /* Generate random numbers for bootstrapping */
    size_t n_unifs = (size_t)nseq * TAX_NBOOT * (max_arraylen / 8 + 1);
    double *unifs = (double *)malloc(n_unifs * sizeof(double));
    std::mt19937 rng(42);
    std::uniform_real_distribution<double> udist(0.0, 1.0);
    for (size_t i = 0; i < n_unifs; i++) {
        unifs[i] = udist(rng);
    }

    /* Classify each query */
    if (verbose) {
        fprintf(stderr, "[INFO] Classifying %d query sequences...\n", nseq);
    }

    #pragma omp parallel for schedule(dynamic, 1)
    for (int j = 0; j < nseq; j++) {
        int karray[10000];
        int bootarray[10000 / 8 + 1];
        float logp;

        size_t seqlen = strlen(seqs[j]);
        if (seqlen < 50) {
            result->rval[j] = 0;  /* NA */
            for (int lev = 0; lev < nlevel; lev++) {
                result->rboot[j * nlevel + lev] = 0;
            }
            continue;
        }

        unsigned int arraylen = tax_karray_build(seqs[j], k, karray);
        int max_g = get_best_genus(karray, &logp, arraylen, n_kmers, ngenus, lgk_probability);
        result->rval[j] = max_g + 1;  /* 1-indexed */

        /* Bootstrap */
        size_t unif_offset = (size_t)j * TAX_NBOOT * (max_arraylen / 8 + 1);
        size_t booti = 0;
        unsigned int sub_len = arraylen / 8;
        if (sub_len == 0) sub_len = 1;

        for (int boot = 0; boot < TAX_NBOOT; boot++) {
            for (unsigned int bi = 0; bi < sub_len; bi++, booti++) {
                int idx = (int)(arraylen * unifs[unif_offset + booti]);
                if (idx >= (int)arraylen) idx = arraylen - 1;
                bootarray[bi] = karray[idx];
            }
            int boot_g = get_best_genus(bootarray, &logp, sub_len, n_kmers, ngenus, lgk_probability);

            /* Count agreement at each level */
            for (int lev = 0; lev < nlevel; lev++) {
                if (genusmat[boot_g * nlevel + lev] == genusmat[max_g * nlevel + lev]) {
                    result->rboot[j * nlevel + lev]++;
                } else {
                    break;
                }
            }
        }
    }

    free(unifs);
    free(lgk_probability);

    if (verbose) {
        fprintf(stderr, "[INFO] Taxonomy classification complete.\n");
    }

    return result;
}

extern "C"
void dada2_tax_result_free(TaxResult *res) {
    if (res) {
        free(res->rval);
        free(res->rboot);
        free(res);
    }
}
