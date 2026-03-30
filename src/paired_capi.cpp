/*
 * paired_capi.cpp - C API for paired-read merging functions.
 * Ports C_nwalign, C_eval_pair, C_pair_consensus from evaluate.cpp
 * and adds reverse complement (rc).
 *
 * Compiled only with -DNO_RCPP.
 */
#ifndef NO_RCPP
/* This file is only for standalone builds. */
#else

#include <cstring>
#include <cstdlib>
#include <cstdio>
#include "dada.h"
#include "dada2_capi.h"

/* ----------------------------------------------------------------
 * dada2_nwalign: Needleman-Wunsch ends-free alignment of two ACGT strings.
 * Caller must free *al1_out and *al2_out with dada2_free_string().
 * Returns 0 on success, -1 on error.
 * ---------------------------------------------------------------- */
extern "C" int dada2_nwalign(const char *s1, const char *s2,
                              int match, int mismatch, int gap_p, int band,
                              char **al1_out, char **al2_out)
{
    if (!s1 || !s2 || !al1_out || !al2_out) return -1;

    size_t len1 = strlen(s1);
    size_t len2 = strlen(s2);

    /* Convert ACGT -> integer encoding */
    char *seq1 = (char *)malloc(len1 + 1);
    char *seq2 = (char *)malloc(len2 + 1);
    if (!seq1 || !seq2) { free(seq1); free(seq2); return -1; }

    nt2int(seq1, s1);
    nt2int(seq2, s2);

    /* Build score matrix */
    int c_score[4][4];
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            c_score[i][j] = (i == j) ? match : mismatch;
        }
    }

    /* Perform ends-free alignment */
    char **al = nwalign_endsfree(seq1, len1, seq2, len2, c_score, gap_p, band);

    /* Convert back to ACGT */
    int2nt(al[0], al[0]);
    int2nt(al[1], al[1]);

    /* Return copies */
    *al1_out = strdup(al[0]);
    *al2_out = strdup(al[1]);

    free(seq1);
    free(seq2);
    free(al[0]);
    free(al[1]);
    free(al);

    return 0;
}

/* ----------------------------------------------------------------
 * dada2_eval_pair: Count matches, mismatches, indels in an alignment.
 * Skips initial and terminal gaps (ends-free evaluation).
 * ---------------------------------------------------------------- */
extern "C" void dada2_eval_pair(const char *al1, const char *al2,
                                 int *out_match, int *out_mismatch, int *out_indel)
{
    *out_match = 0;
    *out_mismatch = 0;
    *out_indel = 0;

    if (!al1 || !al2) return;

    size_t len = strlen(al1);
    if (len != strlen(al2)) {
        fprintf(stderr, "Warning: Aligned strings are not the same length.\n");
        return;
    }
    if (len == 0) return;

    /* Find start of internal alignment (skip initial gaps) */
    bool s1gap = true, s2gap = true;
    int start = -1;
    do {
        start++;
        s1gap = s1gap && (al1[start] == '-');
        s2gap = s2gap && (al2[start] == '-');
    } while ((s1gap || s2gap) && start < (int)len);

    /* Find end of internal alignment (skip terminal gaps) */
    s1gap = s2gap = true;
    int end = (int)len;
    do {
        end--;
        s1gap = s1gap && (al1[end] == '-');
        s2gap = s2gap && (al2[end] == '-');
    } while ((s1gap || s2gap) && end >= start);

    /* Count matches, mismatches, indels in the internal alignment */
    int nmatch = 0, nmismatch = 0, nindel = 0;
    for (int i = start; i <= end; i++) {
        if (al1[i] == '-' || al2[i] == '-') {
            nindel++;
        } else if (al1[i] == al2[i]) {
            nmatch++;
        } else {
            nmismatch++;
        }
    }

    *out_match = nmatch;
    *out_mismatch = nmismatch;
    *out_indel = nindel;
}

/* ----------------------------------------------------------------
 * dada2_pair_consensus: Build consensus from two aligned strings.
 * prefer=1: s1 wins mismatches; prefer=2: s2 wins.
 * trim_overhang: replace overhanging portions with gaps, then strip all gaps.
 * Caller must free result with dada2_free_string().
 * ---------------------------------------------------------------- */
extern "C" char *dada2_pair_consensus(const char *al1, const char *al2,
                                       int prefer, int trim_overhang)
{
    if (!al1 || !al2) return NULL;

    size_t len = strlen(al1);
    if (len != strlen(al2)) {
        fprintf(stderr, "Warning: Aligned strings are not the same length.\n");
        return NULL;
    }

    char *oseq = (char *)malloc(len + 1);
    if (!oseq) return NULL;

    for (size_t i = 0; i < len; i++) {
        if (al1[i] == al2[i]) {
            oseq[i] = al1[i];
        } else if (al2[i] == '-') {
            oseq[i] = al1[i];
        } else if (al1[i] == '-') {
            oseq[i] = al2[i];
        } else {
            if (prefer == 1) {
                oseq[i] = al1[i];
            } else if (prefer == 2) {
                oseq[i] = al2[i];
            } else {
                oseq[i] = 'N';
            }
        }
    }

    /* Trim overhangs */
    if (trim_overhang) {
        /* Mark leading overhang (where al1 has gaps) */
        for (size_t i = 0; i < len; i++) {
            if (al1[i] != '-') break;
            oseq[i] = '-';
        }
        /* Mark trailing overhang (where al2 has gaps) */
        for (int i = (int)len - 1; i >= 0; i--) {
            if (al2[i] != '-') break;
            oseq[i] = '-';
        }
    }

    /* Strip all gap characters */
    size_t j = 0;
    for (size_t i = 0; i < len; i++) {
        if (oseq[i] != '-') {
            oseq[j++] = oseq[i];
        }
    }
    oseq[j] = '\0';

    return oseq;
}

/* ----------------------------------------------------------------
 * dada2_rc: Reverse complement an ACGT string.
 * Caller must free result with dada2_free_string().
 * ---------------------------------------------------------------- */
extern "C" char *dada2_rc(const char *seq)
{
    if (!seq) return NULL;

    size_t len = strlen(seq);
    char *out = (char *)malloc(len + 1);
    if (!out) return NULL;

    for (size_t i = 0; i < len; i++) {
        char c = seq[len - 1 - i];
        switch (c) {
            case 'A': out[i] = 'T'; break;
            case 'T': out[i] = 'A'; break;
            case 'C': out[i] = 'G'; break;
            case 'G': out[i] = 'C'; break;
            case 'N': out[i] = 'N'; break;
            default:  out[i] = 'N'; break;
        }
    }
    out[len] = '\0';
    return out;
}

/* ----------------------------------------------------------------
 * dada2_free_string: Free a string returned by the above functions.
 * ---------------------------------------------------------------- */
extern "C" void dada2_free_string(char *s)
{
    free(s);
}

#endif /* NO_RCPP */
