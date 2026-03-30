"""High-level chimera detection functions matching R's dada2 interface."""

import numpy as np
from ._cdada import is_bimera, table_bimera


def is_bimera_denovo(seqtab_row, seqs, allow_one_off=False,
                     min_one_off_par_dist=4, min_fold=1.5, min_abund=2,
                     match=5, mismatch=-4, gap_p=-8, max_shift=16):
    """Check whether a single sequence is a bimera of more-abundant parents.

    This mirrors R's isBimeraDenovo for a single sample (one row of the
    sequence table).

    Args:
        seqtab_row: 1-D array-like of int abundances, one per ASV.
        seqs: list of ASV sequences (str, ACGT), same length as seqtab_row.
        allow_one_off: allow one mismatch in chimera model.
        min_one_off_par_dist: min hamming distance between parents for one-off.
        min_fold: parent must be this-fold more abundant than the query.
        min_abund: parent minimum absolute abundance.
        match, mismatch, gap_p, max_shift: NW alignment parameters.

    Returns:
        numpy bool array (n_seqs,): True where ASV is flagged as bimera.
    """
    row = np.asarray(seqtab_row, dtype=np.int32)
    n = len(seqs)
    if n == 0:
        return np.array([], dtype=bool)

    # Sort by decreasing abundance for parent ordering
    order = np.argsort(-row)
    flags = np.zeros(n, dtype=bool)

    for idx in range(n):
        j = order[idx]
        if row[j] <= 0:
            continue
        # Gather parents: more abundant sequences
        parents = []
        for k in order:
            if row[k] > min_fold * row[j] and row[k] >= min_abund:
                parents.append(seqs[k])
        if len(parents) == 0:
            continue
        flags[j] = is_bimera(seqs[j], parents,
                             allow_one_off=allow_one_off,
                             min_one_off_par_dist=min_one_off_par_dist,
                             match=match, mismatch=mismatch,
                             gap_p=gap_p, max_shift=max_shift)
    return flags


def remove_bimera_denovo(seqtab, method="consensus", min_fold=1.5,
                         min_abund=2, allow_one_off=False,
                         min_one_off_par_dist=4, min_sample_fraction=0.9,
                         ignore_n_negatives=1,
                         match=5, mismatch=-4, gap_p=-8, max_shift=16,
                         verbose=False):
    """Remove bimeric sequences from a sequence table.

    Mirrors R's removeBimeraDenovo.

    Args:
        seqtab: dict with keys:
            "table": numpy int32 array (samples x ASVs), column-major preferred.
            "seqs": list of ASV sequences (str, ACGT), length = ncol.
        method: "consensus" (default), "pooled", or "per-sample".
            - "consensus": flag per-sample, then remove ASVs flagged in
              enough samples (controlled by min_sample_fraction and
              ignore_n_negatives).
            - "pooled": sum across samples, treat as single sample.
            - "per-sample": zero only the sample/ASV cells flagged as
              chimeric, then drop all-zero ASV columns.
        min_fold: parent fold-abundance threshold.
        min_abund: parent minimum absolute abundance.
        allow_one_off: allow one mismatch in chimera model.
        min_one_off_par_dist: min hamming distance for one-off parents.
        min_sample_fraction: fraction of present samples that must flag
            chimeric for consensus removal (default 0.9).
        ignore_n_negatives: ignore this many non-flagging samples
            (default 1). An ASV is chimeric if
            nflag >= nsam - ignore_n_negatives, provided
            nflag/nsam >= min_sample_fraction.
        match, mismatch, gap_p, max_shift: NW alignment parameters.
        verbose: print progress information.

    Returns:
        dict with:
            "table": numpy int32 array with chimeric columns removed.
            "seqs": list of non-chimeric ASV sequences.
            "is_chimera": numpy bool array (ncol,) for "pooled"/"consensus",
              or numpy bool array (nrow, ncol) for "per-sample".
    """
    if isinstance(seqtab, dict):
        mat = seqtab["table"]
        seqs = list(seqtab["seqs"])
    else:
        raise TypeError("seqtab must be a dict with 'table' and 'seqs' keys")

    mat = np.asarray(mat, dtype=np.int32)
    if mat.ndim == 1:
        mat = mat.reshape(1, -1)
    nrow, ncol = mat.shape

    if ncol != len(seqs):
        raise ValueError("Number of sequences ({}) must match number of "
                         "columns ({})".format(len(seqs), ncol))

    if method == "pooled":
        # Sum across samples, treat as single-sample
        pooled = mat.sum(axis=0).reshape(1, -1).astype(np.int32)
        flags = is_bimera_denovo(
            pooled[0], seqs,
            allow_one_off=allow_one_off,
            min_one_off_par_dist=min_one_off_par_dist,
            min_fold=min_fold, min_abund=min_abund,
            match=match, mismatch=mismatch, gap_p=gap_p, max_shift=max_shift
        )
        is_chimera = flags
    elif method == "consensus":
        # Use the table-level C function for consensus detection
        mat_f = np.asfortranarray(mat, dtype=np.int32)
        result = table_bimera(
            mat_f, seqs,
            min_fold=min_fold, min_abund=min_abund,
            allow_one_off=allow_one_off,
            min_one_off_par_dist=min_one_off_par_dist,
            match=match, mismatch=mismatch, gap_p=gap_p, max_shift=max_shift
        )
        nflag = result["nflag"]
        nsam = result["nsam"]

        # Apply consensus logic (matching R's isBimeraDenovoTable)
        # R logic: flag if nflag > 0 AND
        #   (nflag >= nsam OR nflag >= (nsam - ignoreNNegatives) * minSampleFraction)
        is_chimera = np.zeros(ncol, dtype=bool)
        for j in range(ncol):
            if nflag[j] <= 0 or nsam[j] <= 0:
                continue
            if nflag[j] >= nsam[j]:
                is_chimera[j] = True
            elif nflag[j] >= (nsam[j] - ignore_n_negatives) * min_sample_fraction:
                is_chimera[j] = True
    elif method == "per-sample":
        per_sample = np.zeros((nrow, ncol), dtype=bool)
        new_mat = mat.copy()
        for i in range(nrow):
            flags = is_bimera_denovo(
                new_mat[i], seqs,
                allow_one_off=allow_one_off,
                min_one_off_par_dist=min_one_off_par_dist,
                min_fold=min_fold, min_abund=min_abund,
                match=match, mismatch=mismatch, gap_p=gap_p, max_shift=max_shift
            )
            per_sample[i] = flags
            new_mat[i, flags] = 0

        keep = new_mat.sum(axis=0) > 0
        new_mat = new_mat[:, keep]
        new_seqs = [s for s, k in zip(seqs, keep) if k]
        return {
            "table": new_mat,
            "seqs": new_seqs,
            "is_chimera": per_sample,
        }
    else:
        raise ValueError("method must be 'consensus', 'pooled', or 'per-sample'")

    if verbose:
        n_chim = int(is_chimera.sum())
        n_total = ncol
        print("Identified {} bimeras out of {} input sequences.".format(
            n_chim, n_total))

    # Filter out chimeric columns
    keep = ~is_chimera
    new_mat = mat[:, keep].copy()
    new_seqs = [s for s, k in zip(seqs, keep) if k]

    return {
        "table": new_mat,
        "seqs": new_seqs,
        "is_chimera": is_chimera,
    }
