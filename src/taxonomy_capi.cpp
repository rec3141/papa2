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
#include <cstdint>
#include <cstring>
#include <cmath>
#include <cfloat>
#include <cstdio>
#include <random>
#include <algorithm>

extern "C" void r_rng_runif_fill(uint32_t seed, double *out, long long n);

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
    /* Tie-breaking RNG, constructed lazily: exact logp ties are rare and
     * std::random_device construction is costly in this hot path.  Ties
     * are broken non-deterministically, as in R dada2. */
    std::mt19937 *gen = nullptr;
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
            if (!gen) {
                std::random_device rd;
                gen = new std::mt19937(rd());
            }
            if (cunif(*gen) < 1.0 / nmax) {
                max_g = (int)g;
            }
        }
    }
    delete gen;
    *out_logp = max_logp;
    return max_g;
}


/* Score genera from a pre-gathered probability matrix S (ngenus x
 * arraylen, row-major): S[g*arraylen + i] == lgk_probability[g*n_kmers +
 * karray[i]].  Iterates the same values in the same order with the same
 * early break as get_best_genus, so results are bit-identical — only the
 * memory layout changes (contiguous rows instead of 65536-wide strides).
 * positions == NULL scores the full array; otherwise scores
 * S[g*arraylen + positions[p]] for p in 0..npos. */
static int get_best_genus_gathered(const float *S, float *out_logp,
                                   unsigned int arraylen,
                                   const int *positions, unsigned int npos,
                                   unsigned int ngenus) {
    int max_g = -1;
    float max_logp = -FLT_MAX;
    unsigned int nmax = 0;
    std::mt19937 *gen = nullptr;
    std::uniform_real_distribution<> cunif(0.0, 1.0);

    for (unsigned int g = 0; g < ngenus; g++) {
        const float *row = S + (size_t)g * arraylen;
        float logp = 0.0f;
        if (positions == NULL) {
            for (unsigned int i = 0; i < arraylen; i++) {
                logp += row[i];
                if (logp < max_logp) break;
            }
        } else {
            for (unsigned int p = 0; p < npos; p++) {
                logp += row[positions[p]];
                if (logp < max_logp) break;
            }
        }

        if (max_logp < -FLT_MAX + 1 || logp > max_logp) {
            max_logp = logp;
            max_g = (int)g;
            nmax = 1;
        } else if (logp == max_logp) {
            nmax++;
            if (!gen) {
                std::random_device rd;
                gen = new std::mt19937(rd());
            }
            if (cunif(*gen) < 1.0 / nmax) {
                max_g = (int)g;
            }
        }
    }
    delete gen;
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
    int verbose, long long seed)
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

    /* Build per-genus kmer probability table.  Parallelised by genus so
     * no two threads touch the same lgk row; kmer priors accumulate into
     * per-thread buffers.  All accumulated values are integer counts
     * (exactly representable in float), so the result is bit-identical
     * to the serial order. */
    float *kmer_prior = (float *)calloc(n_kmers, sizeof(float));
    float *lgk_probability = (float *)calloc(ngenus * n_kmers, sizeof(float));

    if (!genus_num_plus1 || !kmer_prior || !lgk_probability) {
        fprintf(stderr, "[ERROR] Memory allocation failed in taxonomy.\n");
        return NULL;
    }

    {
        /* Bucket references by genus (counting sort, preserves order) */
        int *genus_start = (int *)calloc(ngenus + 1, sizeof(int));
        int *ref_by_genus = (int *)malloc(nref * sizeof(int));
        for (int i = 0; i < nref; i++) genus_start[ref_to_genus[i] + 1]++;
        for (int g = 0; g < ngenus; g++) genus_start[g + 1] += genus_start[g];
        {
            int *fill = (int *)malloc(ngenus * sizeof(int));
            memcpy(fill, genus_start, ngenus * sizeof(int));
            for (int i = 0; i < nref; i++) {
                ref_by_genus[fill[ref_to_genus[i]]++] = i;
            }
            free(fill);
        }

        #pragma omp parallel
        {
            unsigned char *ref_kv = (unsigned char *)malloc(n_kmers);
            float *prior_local = (float *)calloc(n_kmers, sizeof(float));

            #pragma omp for schedule(dynamic, 16)
            for (int g = 0; g < ngenus; g++) {
                float *lgk_v = &lgk_probability[(size_t)g * n_kmers];
                for (int p = genus_start[g]; p < genus_start[g + 1]; p++) {
                    tax_kvec_build(refs[ref_by_genus[p]], k, ref_kv);
                    for (size_t km = 0; km < n_kmers; km++) {
                        if (ref_kv[km]) {
                            lgk_v[km]++;
                            prior_local[km]++;
                        }
                    }
                }
            }

            #pragma omp critical
            for (size_t km = 0; km < n_kmers; km++) {
                kmer_prior[km] += prior_local[km];
            }

            free(prior_local);
            free(ref_kv);
        }
        free(genus_start);
        free(ref_by_genus);
    }

    /* Compute kmer priors and log-probabilities.  The double-precision
     * literals match upstream dada2: the expression promotes to double
     * before truncating back to float, and the bootstrap counts are
     * sensitive to that last-ulp difference. */
    for (size_t km = 0; km < n_kmers; km++) {
        kmer_prior[km] = (kmer_prior[km] + 0.5) / (1.0 + nref);
    }
    #pragma omp parallel for schedule(static)
    for (int g = 0; g < ngenus; g++) {
        float *lgk_v = &lgk_probability[(size_t)g * n_kmers];
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

    /* Random numbers for bootstrapping.  Count and layout replicate R
     * dada2 1.40: one flat runif(nseq*NBOOT*(max_arraylen/8)) vector,
     * indexed per query with a stride of max_arraylen.  With seed >= 0
     * the values come from R's own RNG stream (set.seed-compatible), so
     * bootstrap values are identical to a seeded R session. */
    size_t n_unifs = (size_t)nseq * TAX_NBOOT * (max_arraylen / 8);
    if (n_unifs == 0) n_unifs = 1;
    double *unifs = (double *)malloc(n_unifs * sizeof(double));
    if (seed >= 0) {
        r_rng_runif_fill((uint32_t)seed, unifs, (long long)n_unifs);
    } else {
        std::random_device rd;
        std::mt19937 rng(rd());
        std::uniform_real_distribution<double> udist(0.0, 1.0);
        for (size_t i = 0; i < n_unifs; i++) {
            unifs[i] = udist(rng);
        }
    }

    /* Classify each query */
    if (verbose) {
        fprintf(stderr, "[INFO] Classifying %d query sequences...\n", nseq);
    }

    #pragma omp parallel
    {
        /* Per-thread gather buffer: S[g][i] = lgk[g][karray[i]] */
        float *S = (float *)malloc((size_t)ngenus * max_arraylen * sizeof(float));

        #pragma omp for schedule(dynamic, 1)
        for (int j = 0; j < nseq; j++) {
            int karray[10000];
            int bootpos[10000 / 8 + 1];
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
            int max_g;
            if (S != NULL && arraylen > 0) {
                for (unsigned int g = 0; g < (unsigned int)ngenus; g++) {
                    const float *lgk_v = &lgk_probability[(size_t)g * n_kmers];
                    float *row = S + (size_t)g * arraylen;
                    for (unsigned int i = 0; i < arraylen; i++) {
                        row[i] = lgk_v[karray[i]];
                    }
                }
                max_g = get_best_genus_gathered(S, &logp, arraylen, NULL,
                                                arraylen, ngenus);
            } else {
                max_g = get_best_genus(karray, &logp, arraylen, n_kmers,
                                       ngenus, lgk_probability);
            }
            result->rval[j] = max_g + 1;  /* 1-indexed */

            /* Bootstrap — indexing replicates R dada2 1.40: per-query offset
             * j*max_arraylen into the flat unifs vector, subsample size
             * arraylen/8 with a running index across bootstrap replicates. */
            const double *unifs_j = unifs + (size_t)j * max_arraylen;
            size_t booti = 0;
            unsigned int sub_len = arraylen / 8;

            for (int boot = 0; boot < TAX_NBOOT; boot++) {
                int boot_g;
                if (S != NULL && arraylen > 0) {
                    for (unsigned int bi = 0; bi < sub_len; bi++, booti++) {
                        bootpos[bi] = (int)(arraylen * unifs_j[booti]);
                    }
                    boot_g = get_best_genus_gathered(S, &logp, arraylen,
                                                     bootpos, sub_len, ngenus);
                } else {
                    int bootarray[10000 / 8 + 1];
                    for (unsigned int bi = 0; bi < sub_len; bi++, booti++) {
                        bootarray[bi] = karray[(int)(arraylen * unifs_j[booti])];
                    }
                    boot_g = get_best_genus(bootarray, &logp, sub_len, n_kmers,
                                            ngenus, lgk_probability);
                }

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

        free(S);
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
