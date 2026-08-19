"""High-level paired-read merging, mirroring R's mergePairs()."""

import os

import numpy as np
from collections import Counter
from ._cdada import nwalign, eval_pair, pair_consensus, rc


def _merge_one(args):
    """Worker for parallel multi-sample merge_pairs (picklable)."""
    dF, drF, dR, drR, kw = args
    return merge_pairs(dF, drF, dR, drR, **kw)


def merge_pairs(dadaF, derepF, dadaR, derepR,
                min_overlap=12, max_mismatch=0,
                return_rejects=False,
                trim_overhang=False, just_concatenate=False,
                verbose=False):
    """Merge denoised forward and reverse reads into full amplicon sequences.

    Mirrors the R dada2::mergePairs() function.

    Args:
        dadaF: dict from dada() with keys 'denoised', 'cluster_seqs',
               'cluster_abunds', 'map' (0-indexed cluster assignment per unique).
        derepF: dict from derep_fastq() with key 'map' (0-indexed unique
                assignment per read).
        dadaR: same structure as dadaF, for reverse reads.
        derepR: same structure as derepF, for reverse reads.
        min_overlap: minimum overlap required for merging (default 12).
        max_mismatch: maximum mismatches allowed in the overlap (default 0).
        return_rejects: if True, rejected pairs are kept in the output with
            'accept': False and an empty sequence (R's returnRejects).
        trim_overhang: if True, trim overhanging ends of the merged sequence
                       (default False).
        just_concatenate: if True, concatenate rather than merge (inserts 10 Ns
                          between forward and reverse; default False).
        verbose: print progress information (default False).

    Returns:
        List of dicts sorted by abundance (descending), each with keys:
            'sequence': merged sequence (str)
            'abundance': number of reads supporting this pair
            'forward': 0-indexed forward cluster index
            'reverse': 0-indexed reverse cluster index
            'nmatch': number of matches in overlap
            'nmismatch': number of mismatches in overlap
            'nindel': number of indels in overlap
            'prefer': which read's base wins mismatch positions (1=forward,
                      2=reverse), chosen by the larger denoised n0
                      (R's mergePairs behaviour); None when concatenating
            'accept': whether the pair passed the overlap/mismatch criteria

    Only accepted pairs are returned unless return_rejects is True.

    Like R's mergePairs, all four leading arguments may be lists (one
    entry per sample); a list of per-sample results is then returned,
    computed in parallel across samples.
    """
    if isinstance(dadaF, list):
        if not (isinstance(dadaR, list) and isinstance(derepF, list)
                and isinstance(derepR, list)
                and len(dadaF) == len(derepF) == len(dadaR) == len(derepR)):
            raise ValueError("dadaF/derepF/dadaR/derepR must be lists of "
                             "the same length.")
        kw = dict(min_overlap=min_overlap, max_mismatch=max_mismatch,
                  return_rejects=return_rejects, trim_overhang=trim_overhang,
                  just_concatenate=just_concatenate, verbose=verbose)
        tasks = list(zip(dadaF, derepF, dadaR, derepR, [kw] * len(dadaF)))
        n_workers = min(len(tasks), os.cpu_count() or 1)
        if n_workers > 1:
            try:
                from loky import get_reusable_executor
                with get_reusable_executor(max_workers=n_workers) as pool:
                    return list(pool.map(_merge_one, tasks))
            except ImportError:
                pass
        return [_merge_one(t) for t in tasks]

    # --- Step 1: Map each read to its denoised cluster ---
    # derepF['map'] maps read -> unique index (0-indexed)
    # dadaF['map'] maps unique index -> cluster index (0-indexed, -1 = unassigned)
    dada_map_F = np.asarray(dadaF["map"], dtype=np.int32)
    derep_map_F = np.asarray(derepF["map"], dtype=np.int32)
    dada_map_R = np.asarray(dadaR["map"], dtype=np.int32)
    derep_map_R = np.asarray(derepR["map"], dtype=np.int32)

    # Per-read cluster assignment: rF[i] = cluster of read i
    rF = dada_map_F[derep_map_F]
    rR = dada_map_R[derep_map_R]

    # --- Step 2: Find unique (forward, reverse) cluster pairs with counts ---
    # Only keep reads that were assigned to a cluster (not -1)
    valid = (rF >= 0) & (rR >= 0)
    pairs = list(zip(rF[valid].tolist(), rR[valid].tolist()))
    pair_counts = Counter(pairs)

    if verbose:
        print(f"Found {len(pair_counts)} unique F/R cluster pairs from "
              f"{sum(pair_counts.values())} reads.")

    # Get cluster sequences
    fwd_seqs = dadaF["cluster_seqs"]
    rev_seqs = dadaR["cluster_seqs"]

    # n0 per cluster decides which read is preferred at mismatch positions
    # (R: prefer = 1 + (dadaR$clustering$n0 > dadaF$clustering$n0))
    n0_F = np.asarray(dadaF.get("cluster_n0", np.zeros(len(fwd_seqs))), dtype=np.int64)
    n0_R = np.asarray(dadaR.get("cluster_n0", np.zeros(len(rev_seqs))), dtype=np.int64)

    # Alignment scoring: stringent defaults for merging
    # match=1, mismatch=-64, gap=-64 makes NW alignment find the overlap
    # with effectively zero tolerance (for max_mismatch=0)
    if max_mismatch == 0:
        nw_match = 1
        nw_mismatch = -64
        nw_gap = -64
    else:
        # More permissive alignment for non-zero max_mismatch
        nw_match = 1
        nw_mismatch = -8
        nw_gap = -8

    # --- Step 3: Process each unique pair ---
    results = []
    for (fi, ri), abundance in pair_counts.items():
        seq_F = fwd_seqs[fi]
        seq_R_rc = rc(rev_seqs[ri])

        if just_concatenate:
            merged = seq_F + "N" * 10 + seq_R_rc
            results.append({
                "sequence": merged,
                "abundance": abundance,
                "forward": fi,
                "reverse": ri,
                "nmatch": 0,
                "nmismatch": 0,
                "nindel": 0,
                "prefer": None,
                "accept": True,
            })
            continue

        # NW ends-free alignment
        al1, al2 = nwalign(seq_F, seq_R_rc,
                            match=nw_match, mismatch=nw_mismatch,
                            gap_p=nw_gap, band=-1)

        # Evaluate the alignment
        nmatch, nmismatch, nindel = eval_pair(al1, al2)

        # Accept/reject based on criteria
        accept = (nmatch >= min_overlap) and ((nmismatch + nindel) <= max_mismatch)

        prefer = 2 if n0_R[ri] > n0_F[fi] else 1
        if accept:
            merged = pair_consensus(al1, al2, prefer=prefer,
                                    trim_overhang=trim_overhang)
        else:
            merged = ""

        if verbose and not accept:
            print(f"Pair F{fi}/R{ri} (abd={abundance}): "
                  f"{nmatch} match, {nmismatch} mismatch, {nindel} indel -- REJECTED")

        results.append({
            "sequence": merged,
            "abundance": abundance,
            "forward": fi,
            "reverse": ri,
            "nmatch": nmatch,
            "nmismatch": nmismatch,
            "nindel": nindel,
            "prefer": prefer,
            "accept": accept,
        })

    # --- Step 4: Sort by abundance descending (stable, like R's order) ---
    results.sort(key=lambda x: x["abundance"], reverse=True)

    if not return_rejects:
        results = [r for r in results if r["accept"]]

    if verbose:
        n_accept = sum(1 for r in results if r["accept"])
        n_total = len(results)
        total_reads = sum(r["abundance"] for r in results if r["accept"])
        print(f"Accepted {n_accept}/{n_total} unique pairs, "
              f"representing {total_reads} merged reads.")

    return results
