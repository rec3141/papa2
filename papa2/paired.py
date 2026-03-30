"""High-level paired-read merging, mirroring R's mergePairs()."""

import numpy as np
from collections import Counter
from ._cdada import nwalign, eval_pair, pair_consensus, rc


def merge_pairs(dadaF, derepF, dadaR, derepR,
                min_overlap=12, max_mismatch=0,
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
            'accept': whether the pair passed the overlap/mismatch criteria
    """
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

        if accept:
            # Build consensus (prefer forward sequence on mismatch)
            merged = pair_consensus(al1, al2, prefer=1,
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
            "accept": accept,
        })

    # --- Step 4: Sort by abundance descending ---
    results.sort(key=lambda x: x["abundance"], reverse=True)

    if verbose:
        n_accept = sum(1 for r in results if r["accept"])
        n_total = len(results)
        total_reads = sum(r["abundance"] for r in results if r["accept"])
        print(f"Accepted {n_accept}/{n_total} unique pairs, "
              f"representing {total_reads} merged reads.")

    return results
