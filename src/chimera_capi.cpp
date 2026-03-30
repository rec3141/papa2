/*
 * chimera_capi.cpp — Standalone C API for dada2 chimera detection.
 *
 * Ports C_is_bimera() and C_table_bimera2() from chimera.cpp,
 * replacing Rcpp/RcppParallel with plain C types and OpenMP.
 */

#include "dada.h"
#include "dada2_capi.h"

#include <string.h>
#include <stdlib.h>
#include <vector>
#include <string>

/* ---- internal helpers (same logic as original chimera.cpp) ---- */

static int get_ham_endsfree(const char *seq1, const char *seq2) {
  size_t len = strlen(seq2);
  bool gap1, gap2;
  // Find start of internal part of alignment
  int i = 0;
  gap1 = (seq1[i] == '-');
  gap2 = (seq2[i] == '-');
  while (gap1 || gap2) {
    i++;
    gap1 = (gap1 && (seq1[i] == '-'));
    gap2 = (gap2 && (seq2[i] == '-'));
  }
  // Find end of internal part of alignment
  int j = (int)len - 1;
  gap1 = (seq1[j] == '-');
  gap2 = (seq2[j] == '-');
  while (gap1 || gap2) {
    j--;
    gap1 = (gap1 && (seq1[j] == '-'));
    gap2 = (gap2 && (seq2[j] == '-'));
  }
  // Calculate hamming distance over internal part
  int ham = 0;
  for (int pos = i; pos <= j; pos++) {
    if (seq1[pos] != seq2[pos]) { ham++; }
  }
  return ham;
}

static void get_lr(char **al, int &left, int &right, int &left_oo, int &right_oo, bool allow_one_off, int max_shift) {
  size_t len = strlen(al[0]);
  int pos = 0;
  left = 0;
  while (al[0][pos] == '-' && pos < (int)len) {
    pos++; // Scan in until query starts
  }
  while (al[1][pos] == '-' && pos < max_shift) {
    pos++; left++; // Credit as ends-free coverage until parent starts
  }
  while (pos < (int)len && al[0][pos] == al[1][pos]) {
    pos++; left++; // Credit as covered until a mismatch
  }
  if (allow_one_off) {
    left_oo = left;
    pos++;
    if (pos < (int)len && al[0][pos] != '-') { left_oo++; }
    while (pos < (int)len && al[0][pos] == al[1][pos]) {
      pos++; left_oo++;
    }
  }

  pos = (int)len - 1;
  right = 0;
  while (al[0][pos] == '-' && pos >= 0) {
    pos--;
  }
  while (al[1][pos] == '-' && pos > (int)(len - max_shift)) {
    pos--; right++;
  }
  while (pos >= 0 && al[0][pos] == al[1][pos]) {
    pos--; right++;
  }
  if (allow_one_off) {
    right_oo = right;
    pos--;
    if (pos >= 0 && al[0][pos] != '-') { right_oo++; }
    while (pos >= 0 && al[0][pos] == al[1][pos]) {
      pos--; right_oo++;
    }
  }
}

/* ---- helper: convert ACGT string to int-encoded (in-place alloc) ---- */

static char *acgt_to_int(const char *seq, size_t len) {
  char *buf = (char *)malloc(len + 1);
  nt2int(buf, seq);
  buf[len] = '\0';
  return buf;
}

/* nwalign_vectorized2 expects int-encoded sequences and returns
   int-encoded alignments.  The chimera helpers (get_lr, get_ham_endsfree)
   work on ACGT+'-' strings.  So we convert the alignment output back. */

static char **align_acgt(const char *seq_acgt, size_t len1,
                          const char *par_acgt, size_t len2,
                          int16_t match, int16_t mismatch, int16_t gap_p,
                          int max_shift) {
  char *iseq1 = acgt_to_int(seq_acgt, len1);
  char *iseq2 = acgt_to_int(par_acgt, len2);

  char **al = nwalign_vectorized2(iseq1, len1, iseq2, len2,
                                   match, mismatch, gap_p, 0, max_shift);
  free(iseq1);
  free(iseq2);

  /* Convert alignment back to ACGT+'-' */
  size_t alen = strlen(al[0]);
  for (size_t p = 0; p < alen; p++) {
    if (al[0][p] != '-') {
      switch ((int)(unsigned char)al[0][p]) {
        case 1: al[0][p] = 'A'; break;
        case 2: al[0][p] = 'C'; break;
        case 3: al[0][p] = 'G'; break;
        case 4: al[0][p] = 'T'; break;
        default: break; // keep as-is (gap)
      }
    }
    if (al[1][p] != '-') {
      switch ((int)(unsigned char)al[1][p]) {
        case 1: al[1][p] = 'A'; break;
        case 2: al[1][p] = 'C'; break;
        case 3: al[1][p] = 'G'; break;
        case 4: al[1][p] = 'T'; break;
        default: break;
      }
    }
  }
  return al;
}

/* ================================================================
   dada2_is_bimera — check if seq is a bimera of parents
   ================================================================ */

extern "C"
int dada2_is_bimera(
    const char *seq,
    const char **parents,
    int n_parents,
    int allow_one_off,
    int min_one_off_par_dist,
    int match, int mismatch, int gap_p, int max_shift)
{
  int left, right, left_oo = 0, right_oo = 0;
  int max_left = 0, max_right = 0;
  int oo_max_left = 0, oo_max_right = 0, oo_max_left_oo = 0, oo_max_right_oo = 0;
  bool rval = false;
  size_t sqlen = strlen(seq);

  for (int i = 0; i < n_parents && !rval; i++) {
    size_t plen = strlen(parents[i]);
    char **al = align_acgt(seq, sqlen, parents[i], plen,
                           (int16_t)match, (int16_t)mismatch, (int16_t)gap_p, max_shift);
    get_lr(al, left, right, left_oo, right_oo, (bool)allow_one_off, max_shift);

    if ((left + right) >= (int)sqlen) {
      // Toss id/pure-shift/internal-indel "parents"
      free(al[0]); free(al[1]); free(al);
      continue;
    }
    if (left > max_left) { max_left = left; }
    if (right > max_right) { max_right = right; }

    if (allow_one_off && get_ham_endsfree(al[0], al[1]) >= min_one_off_par_dist) {
      if (left > oo_max_left) { oo_max_left = left; }
      if (right > oo_max_right) { oo_max_right = right; }
      if (left_oo > oo_max_left_oo) { oo_max_left_oo = left_oo; }
      if (right_oo > oo_max_right_oo) { oo_max_right_oo = right_oo; }
    }

    if ((max_right + max_left) >= (int)sqlen) {
      rval = true;
    }
    if (allow_one_off) {
      if ((oo_max_left + oo_max_right_oo) >= (int)sqlen ||
          (oo_max_left_oo + oo_max_right) >= (int)sqlen) {
        rval = true;
      }
    }
    free(al[0]); free(al[1]); free(al);
  }
  return rval ? 1 : 0;
}

/* ================================================================
   dada2_table_bimera — consensus table-level chimera detection
   ================================================================ */

extern "C"
ChimeraResult* dada2_table_bimera(
    const int *mat,
    int nrow, int ncol,
    const char **seqs,
    double min_fold,
    int min_abund,
    int allow_one_off,
    int min_one_off_par_dist,
    int match, int mismatch, int gap_p, int max_shift)
{
  ChimeraResult *res = (ChimeraResult *)malloc(sizeof(ChimeraResult));
  res->n_seqs = ncol;
  res->nflag = (int *)calloc(ncol, sizeof(int));
  res->nsam  = (int *)calloc(ncol, sizeof(int));

  #pragma omp parallel for schedule(dynamic)
  for (int j = 0; j < ncol; j++) {
    int nsam_local = 0, nflag_local = 0;
    int sqlen = (int)strlen(seqs[j]);

    /* Per-parent cached alignment results */
    std::vector<int> lefts(ncol, -1);
    std::vector<int> rights(ncol, -1);
    std::vector<int> lefts_oo(ncol, -1);
    std::vector<int> rights_oo(ncol, -1);
    std::vector<bool> allowed(ncol, false);

    for (int i = 0; i < nrow; i++) {  // each sample (row)
      if (mat[i + j * nrow] <= 0) { continue; }
      nsam_local++;

      int max_left = 0, max_right = 0;
      int oo_max_left = 0, oo_max_right = 0, oo_max_left_oo = 0, oo_max_right_oo = 0;

      for (int k = 0; k < ncol; k++) {  // each potential parent (column)
        if (mat[i + k * nrow] > (min_fold * mat[i + j * nrow]) &&
            mat[i + k * nrow] >= min_abund) {

          if (lefts[k] < 0) {  // alignment not yet done for this parent
            int left, right, left_oo = 0, right_oo = 0;
            char **al = align_acgt(seqs[j], sqlen, seqs[k], strlen(seqs[k]),
                                    (int16_t)match, (int16_t)mismatch, (int16_t)gap_p, max_shift);
            get_lr(al, left, right, left_oo, right_oo, (bool)allow_one_off, max_shift);

            if (allow_one_off && get_ham_endsfree(al[0], al[1]) >= min_one_off_par_dist) {
              allowed[k] = true;
            }

            if ((left + right) < sqlen) {
              lefts[k] = left;
              rights[k] = right;
              if (allow_one_off) {
                lefts_oo[k] = left_oo;
                rights_oo[k] = right_oo;
              }
            } else {  // Ignore id/pure-shift/internal-indel "parents"
              lefts[k] = 0;
              rights[k] = 0;
              if (allow_one_off) {
                lefts_oo[k] = 0;
                rights_oo[k] = 0;
              }
            }
            free(al[0]); free(al[1]); free(al);
          }

          // Compare to best parents yet found
          if (lefts[k] > max_left) { max_left = lefts[k]; }
          if (rights[k] > max_right) { max_right = rights[k]; }
          if (allow_one_off && allowed[k]) {
            if (lefts[k] > oo_max_left) { oo_max_left = lefts[k]; }
            if (rights[k] > oo_max_right) { oo_max_right = rights[k]; }
            if (lefts_oo[k] > oo_max_left_oo) { oo_max_left_oo = lefts_oo[k]; }
            if (rights_oo[k] > oo_max_right_oo) { oo_max_right_oo = rights_oo[k]; }
          }
        }
      }  // for k

      // Flag if chimeric model exists
      if ((max_right + max_left) >= sqlen) {
        nflag_local++;
      } else if (allow_one_off) {
        if ((oo_max_left + oo_max_right_oo) >= sqlen ||
            (oo_max_left_oo + oo_max_right) >= sqlen) {
          nflag_local++;
        }
      }
    }  // for i (samples)

    res->nflag[j] = nflag_local;
    res->nsam[j]  = nsam_local;
  }  // for j (ASVs)

  return res;
}

extern "C"
void dada2_chimera_result_free(ChimeraResult *res) {
  if (res) {
    free(res->nflag);
    free(res->nsam);
    free(res);
  }
}
