/*
 * filter_capi.cpp — C helpers for FASTQ filtering.
 *
 * dada2_match_ref_counts ports C_matchRef from dada2's filter.cpp: count
 * word_size-mers of each query that occur in the reference, treating the
 * reference as circular, optionally skipping ahead by word_size after each
 * hit (non-overlapping counting).  Must reproduce R's counts exactly —
 * isPhiX() and rm.phix filtering depend on it.
 */

#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <unordered_set>

#ifdef _OPENMP
#include <omp.h>
#endif

/* 2-bit encode one base; -1 otherwise.  Uppercase only: R's C_matchRef
 * compares raw kmer strings, so lowercase never matches the (uppercase)
 * reference — a lowercase base here must likewise produce no match. */
static inline int nt2i(char c) {
    switch (c) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return -1;
    }
}

/* Encode the word_size-mer starting at s; returns false if it contains a
 * non-ACGT base.  word_size <= 32. */
static inline bool encode_kmer(const char *s, int word_size, uint64_t *out) {
    uint64_t k = 0;
    for (int j = 0; j < word_size; j++) {
        int v = nt2i(s[j]);
        if (v < 0) return false;
        k = (k << 2) | (uint64_t)v;
    }
    *out = k;
    return true;
}


/* Shared implementation over explicit per-sequence windows. */
static void match_ref_windows_impl(const char *buf, const int64_t *starts,
                                   const int64_t *ends, int nseq,
                                   const char *ref, int word_size,
                                   int non_overlapping, int *out) {
    size_t ref_len = strlen(ref);
    std::unordered_set<uint64_t> phash;
    phash.reserve(ref_len * 2);

    /* Circular reference: kmers starting at every position, wrapping. */
    {
        char *circ = (char *)malloc(ref_len + word_size + 1);
        if (!circ) { for (int i = 0; i < nseq; i++) out[i] = 0; return; }
        memcpy(circ, ref, ref_len);
        memcpy(circ + ref_len, ref, word_size);
        circ[ref_len + word_size] = '\0';
        for (size_t i = 0; i < ref_len; i++) {
            uint64_t k;
            if (encode_kmer(circ + i, word_size, &k)) phash.insert(k);
        }
        free(circ);
    }

    #pragma omp parallel for schedule(dynamic, 256)
    for (int i = 0; i < nseq; i++) {
        const char *seq = buf + starts[i];
        int64_t len = ends[i] - starts[i];
        int count = 0;
        if (len >= word_size) {
            for (int64_t j = 0; j <= len - word_size; j++) {
                uint64_t k;
                if (encode_kmer(seq + j, word_size, &k) && phash.count(k)) {
                    count++;
                    if (non_overlapping) j += word_size;
                }
            }
        }
        out[i] = count;
    }
}

extern "C" {

/*
 * Count reference kmer hits per query sequence.
 *
 * concat:   all query sequences concatenated (no separators)
 * offsets:  int64 array of length nseq+1; query i is concat[offsets[i]..offsets[i+1])
 * ref:      reference sequence (treated as circular, like C_matchRef)
 * word_size: kmer length (<= 32)
 * non_overlapping: if nonzero, skip word_size positions after each hit
 * out:      int array of length nseq to receive hit counts
 */
void dada2_match_ref_counts(const char *concat, const int64_t *offsets,
                            int nseq, const char *ref, int word_size,
                            int non_overlapping, int *out) {
    /* offsets is starts[0..n-1] with offsets[i+1] as each end */
    match_ref_windows_impl(concat, offsets, offsets + 1, nseq, ref,
                           word_size, non_overlapping, out);
}

/* Same, but with explicit per-sequence start/end windows into buf. */
void dada2_match_ref_windows(const char *buf, const int64_t *starts,
                             const int64_t *ends, int nseq, const char *ref,
                             int word_size, int non_overlapping, int *out) {
    match_ref_windows_impl(buf, starts, ends, nseq, ref, word_size,
                           non_overlapping, out);
}

} /* extern "C" */
