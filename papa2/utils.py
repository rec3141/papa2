"""Utility functions for the dada2 Python package.

Pure Python implementations of common dada2 helper functions including
taxonomy assignment, sequence table operations, quality profiling,
FASTA I/O, PhiX detection, and sequence complexity analysis.
"""

import gzip
import logging
import math
import os
from collections import defaultdict
from typing import Dict, List, Optional, Sequence, Tuple, Union

import numpy as np

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# PhiX reference: first 100 bp of the PhiX174 genome (GenBank: NC_001422.1)
# Used for quick kmer-based matching in is_phix().
# ---------------------------------------------------------------------------
_PHIX_100BP = (
    "GAGTTTTATCGCTTCCATGACGCAGAAGTTAACACTTTCGGATATTTCTGATGAGTCGAAAAATTATCTT"
    "GATAAAGCAGGAATTACTACTGCTTGTTTACGA"
)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _rc(seq: str) -> str:
    """Reverse complement a DNA sequence (ACGT only, pure-Python fallback).

    Tries to use the fast C implementation from _cdada first.
    """
    try:
        from . import rc as _c_rc
        return _c_rc(seq)
    except Exception:
        comp = str.maketrans("ACGTacgt", "TGCAtgca")
        return seq.translate(comp)[::-1]


def _open_maybe_gz(path: str):
    """Open a file, auto-detecting gzip compression."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def _parse_fasta(path: str) -> List[Tuple[str, str]]:
    """Parse a FASTA file into (header, sequence) pairs.

    Supports gzip-compressed files.
    """
    records: List[Tuple[str, str]] = []
    header = None
    seq_parts: List[str] = []
    with _open_maybe_gz(path) as fh:
        for line in fh:
            line = line.rstrip("\n\r")
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_parts).upper()))
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line.strip())
        if header is not None:
            records.append((header, "".join(seq_parts).upper()))
    return records


def _sindex(counts) -> float:
    """Shannon effective number of categories (exp of Shannon entropy)."""
    total = sum(counts)
    if total == 0:
        return 0.0
    freqs = [c / total for c in counts if c > 0]
    H = -sum(f * math.log(f) for f in freqs)
    return math.exp(H)


# ---------------------------------------------------------------------------
# 10. get_sequences -- extract sequences from various input types
# ---------------------------------------------------------------------------

def get_sequences(obj, collapse: bool = False) -> List[str]:
    """Extract a list of sequences from various dada2 object types.

    Supported inputs:
        - list / tuple of strings (returned as-is, upper-cased)
        - dict ``{seq: abundance}``  (keys returned)
        - pandas DataFrame with a ``"sequence"`` column
        - pandas DataFrame where column names are sequences (sequence table)
        - numpy character matrix with row names (taxonomy table -- not common
          in Python, included for completeness)
        - a single file path to a FASTA/FASTQ file

    Args:
        obj: Input object.
        collapse: If True, remove duplicate sequences.

    Returns:
        List of upper-case DNA sequence strings.
    """
    # String -- could be a file path or a single sequence
    if isinstance(obj, str):
        if os.path.isfile(obj):
            recs = _parse_fasta(obj)
            return [seq.upper() for _, seq in recs]
        return [obj.upper()]

    # Dict {seq: abundance}
    if isinstance(obj, dict):
        seqs = [str(k).upper() for k in obj.keys()]
        if collapse:
            seqs = list(dict.fromkeys(seqs))
        return seqs

    # Try pandas objects
    try:
        import pandas as pd
        if isinstance(obj, pd.DataFrame):
            if "sequence" in obj.columns:
                seqs = obj["sequence"].astype(str).str.upper().tolist()
            else:
                # Sequence table: column names are sequences
                seqs = [str(c).upper() for c in obj.columns]
            if collapse:
                seqs = list(dict.fromkeys(seqs))
            return seqs
        if isinstance(obj, pd.Series):
            return [str(v).upper() for v in obj.values]
    except ImportError:
        pass

    # Iterable of strings
    if hasattr(obj, "__iter__"):
        seqs = [str(s).upper() for s in obj]
        if collapse:
            seqs = list(dict.fromkeys(seqs))
        return seqs

    raise TypeError(
        "Cannot extract sequences from object of type "
        f"{type(obj).__name__}. Provide a list of strings, a dict "
        "{{seq: abundance}}, a pandas DataFrame, or a FASTA file path."
    )


# ---------------------------------------------------------------------------
# 11. get_uniques -- extract {seq: abundance} dict from various types
# ---------------------------------------------------------------------------

def get_uniques(obj, collapse: bool = True) -> Dict[str, int]:
    """Extract a ``{sequence: abundance}`` dictionary from various types.

    Supported inputs:
        - dict ``{seq: abundance}`` -- returned directly (optionally collapsed)
        - list / tuple of strings -- each occurrence counted
        - pandas DataFrame with ``"sequence"`` and ``"abundance"`` columns
        - pandas DataFrame where columns are sequences (sequence table) --
          column sums used as abundances
        - a single FASTA file path (each record gets abundance 1)

    Args:
        obj: Input object.
        collapse: If True, merge duplicate sequences by summing abundances.

    Returns:
        Dictionary mapping upper-case sequence strings to integer abundances.
    """
    # Already a dict
    if isinstance(obj, dict):
        out: Dict[str, int] = {}
        for k, v in obj.items():
            k_up = str(k).upper()
            out[k_up] = out.get(k_up, 0) + int(v) if collapse else int(v)
        return out

    # String file path
    if isinstance(obj, str) and os.path.isfile(obj):
        recs = _parse_fasta(obj)
        out = {}
        for _, seq in recs:
            seq = seq.upper()
            out[seq] = out.get(seq, 0) + 1
        return out

    # Try pandas
    try:
        import pandas as pd
        if isinstance(obj, pd.DataFrame):
            if "sequence" in obj.columns and "abundance" in obj.columns:
                out = {}
                for seq, ab in zip(obj["sequence"], obj["abundance"]):
                    seq = str(seq).upper()
                    ab = int(ab)
                    if collapse:
                        out[seq] = out.get(seq, 0) + ab
                    else:
                        out[seq] = ab
                return out
            else:
                # Sequence table: column names are seqs, values are counts
                sums = obj.sum(axis=0)
                return {str(seq).upper(): int(ab) for seq, ab in sums.items()}
    except ImportError:
        pass

    # Iterable of strings -- count occurrences
    if hasattr(obj, "__iter__"):
        out = {}
        for s in obj:
            s = str(s).upper()
            out[s] = out.get(s, 0) + 1
        return out

    raise TypeError(
        "Cannot extract uniques from object of type "
        f"{type(obj).__name__}."
    )


# ---------------------------------------------------------------------------
# 1. assign_species
# ---------------------------------------------------------------------------

def assign_species(
    seqs,
    ref_fasta: str,
    allow_multiple: Union[bool, int] = False,
    try_rc: bool = False,
) -> np.ndarray:
    """Taxonomic assignment to species level by exact matching.

    Each query sequence is searched as a substring against reference sequences.
    Reference FASTA headers must be in the format::

        >SeqID Genus species
        ACGAATGTGAAGTAA...

    Args:
        seqs: Sequences to classify (any type accepted by ``get_sequences``).
        ref_fasta: Path to reference FASTA file (may be gzipped).
        allow_multiple: If False, only unambiguous (single) species matches
            are returned.  If True, all matching species are returned
            (concatenated with ``/``).  If an integer, at most that many
            matching species are returned.
        try_rc: If True, also search the reverse complement of each query.

    Returns:
        A numpy character array of shape ``(N, 2)`` with columns
        ``["Genus", "Species"]``.  Rows correspond to input sequences.
        ``None`` entries indicate no match.
    """
    seqs = get_sequences(seqs)
    logger.info("[INFO] assign_species: classifying %d sequences against %s",
                len(seqs), ref_fasta)

    if isinstance(allow_multiple, bool):
        keep = float("inf") if allow_multiple else 1
    else:
        keep = int(allow_multiple)

    # Parse reference
    ref_records = _parse_fasta(ref_fasta)
    ref_seqs: List[str] = []
    ref_genus: List[str] = []
    ref_species: List[str] = []
    for header, rseq in ref_records:
        parts = header.split()
        if len(parts) >= 3:
            ref_genus.append(parts[1])
            ref_species.append(parts[2])
        elif len(parts) == 2:
            ref_genus.append(parts[1])
            ref_species.append("")
        else:
            ref_genus.append("")
            ref_species.append("")
        ref_seqs.append(rseq)

    n = len(seqs)
    result = np.full((n, 2), None, dtype=object)

    # Build a kmer index over reference sequences for fast candidate lookup.
    # Instead of checking every query against every ref (O(n*m)), we only
    # verify substring containment for refs that share a kmer with the query.
    _KMER_SIZE = 8
    kmer_index: Dict[str, List[int]] = {}
    for j, rseq in enumerate(ref_seqs):
        rseq_upper = rseq.upper()
        for pos in range(0, len(rseq_upper) - _KMER_SIZE + 1, _KMER_SIZE):
            kmer = rseq_upper[pos:pos + _KMER_SIZE]
            if kmer not in kmer_index:
                kmer_index[kmer] = []
            kmer_index[kmer].append(j)

    def _find_hits(query_up):
        """Find ref indices containing query_up as a substring, using kmer pre-filter."""
        candidates: set = set()
        for pos in range(len(query_up) - _KMER_SIZE + 1):
            kmer = query_up[pos:pos + _KMER_SIZE]
            if kmer in kmer_index:
                candidates.update(kmer_index[kmer])
        # Verify actual substring containment on candidates only
        return [j for j in candidates if query_up in ref_seqs[j].upper()]

    for i, query in enumerate(seqs):
        query_up = query.upper()
        hit_indices = _find_hits(query_up)
        # Also check reverse complement
        if try_rc and not hit_indices:
            query_rc = _rc(query_up)
            hit_indices = _find_hits(query_rc)

        if not hit_indices:
            continue

        # Map hits to genus/species, collapsing Escherichia/Shigella
        genera = [ref_genus[j] for j in hit_indices]
        species = [ref_species[j] for j in hit_indices]

        # Genus: only unambiguous (single unique genus)
        genera_clean = list(genera)
        for k in range(len(genera_clean)):
            if "Escherichia" in genera_clean[k] or "Shigella" in genera_clean[k]:
                genera_clean[k] = "Escherichia/Shigella"
        unique_genera = sorted(set(genera_clean))
        if len(unique_genera) == 1:
            result[i, 0] = unique_genera[0]
        else:
            result[i, 0] = None

        # Species: apply keep threshold
        species_clean = list(species)
        for k in range(len(species_clean)):
            g = genera[k]
            if "Escherichia" in g or "Shigella" in g:
                species_clean[k] = species_clean[k]  # keep species as-is
        unique_species = sorted(set(s for s in species_clean if s))
        if len(unique_species) == 0:
            result[i, 1] = None
        elif len(unique_species) <= keep:
            result[i, 1] = "/".join(unique_species)
        else:
            result[i, 1] = None

    logger.info("[INFO] assign_species: %d of %d assigned to species level.",
                sum(result[:, 1] != None), n)  # noqa: E711
    return result


# ---------------------------------------------------------------------------
# 2. add_species
# ---------------------------------------------------------------------------

def add_species(
    tax_table,
    ref_fasta: str,
    allow_multiple: Union[bool, int] = False,
    try_rc: bool = False,
):
    """Add species-level annotation to an existing taxonomy DataFrame.

    Wraps ``assign_species`` and appends a ``"Species"`` column.  Only species
    assignments whose genus is consistent with the genus already present in
    ``tax_table`` are kept.

    Args:
        tax_table: A pandas DataFrame with sequences as the index and
            taxonomic levels as columns.  Must include a ``"Genus"`` column
            (or the last column is used).
        ref_fasta: Path to species-level reference FASTA.
        allow_multiple: Passed to ``assign_species``.
        try_rc: Passed to ``assign_species``.

    Returns:
        The input DataFrame with an added ``"Species"`` column.
    """
    import pandas as pd

    seqs = list(tax_table.index)
    binom = assign_species(seqs, ref_fasta,
                           allow_multiple=allow_multiple, try_rc=try_rc)

    # Determine genus column
    if "Genus" in tax_table.columns:
        gcol = "Genus"
    else:
        gcol = tax_table.columns[-1]

    species_col = []
    for i, seq in enumerate(seqs):
        tax_genus = tax_table.loc[seq, gcol] if seq in tax_table.index else None
        binom_genus = binom[i, 0]
        binom_species = binom[i, 1]

        if _match_genera(tax_genus, binom_genus):
            species_col.append(binom_species)
        else:
            species_col.append(None)

    tax_table = tax_table.copy()
    tax_table["Species"] = species_col
    return tax_table


def _match_genera(gen_tax, gen_binom, split_glyph: str = "/") -> bool:
    """Check if a curated genus name matches a binomial genus name.

    Handles split genera (e.g. ``"Clostridium/Ruminococcus"``).
    """
    if gen_tax is None or gen_binom is None:
        return False
    gen_tax = str(gen_tax).strip()
    gen_binom = str(gen_binom).strip()
    if not gen_tax or not gen_binom:
        return False
    if gen_tax == gen_binom:
        return True
    # gen_binom at the start followed by space, underscore, or split glyph
    if gen_tax.startswith(gen_binom) and len(gen_tax) > len(gen_binom):
        next_char = gen_tax[len(gen_binom)]
        if next_char in (" ", "_", split_glyph):
            return True
    # gen_binom at the end after split glyph
    if gen_tax.endswith(split_glyph + gen_binom):
        return True
    return False


# ---------------------------------------------------------------------------
# 3. collapse_no_mismatch
# ---------------------------------------------------------------------------

def collapse_no_mismatch(seqtab) -> dict:
    """Merge ASVs whose sequences are identical except for length differences.

    If sequence A is a substring of sequence B (allowing only leading/trailing
    gaps, i.e. one is a prefix or suffix or internal subsequence with no
    mismatches), their abundances are merged under the longer sequence.

    Sequences are processed in decreasing order of total abundance.  For each
    query, if it is found as a substring of an already-accepted sequence (or
    vice versa), it is collapsed into that representative.

    Args:
        seqtab: Either a ``dict`` mapping ``{sequence: abundance}`` or a
            pandas DataFrame (samples x sequences, column names are seqs).

    Returns:
        If input was a dict, returns a collapsed ``{sequence: abundance}``
        dict.  If input was a DataFrame, returns a collapsed DataFrame.
    """
    import pandas as pd

    is_df = isinstance(seqtab, pd.DataFrame)

    if is_df:
        # Column names are sequences; collapse column-wise
        seq_abundances = seqtab.sum(axis=0).to_dict()
    else:
        seq_abundances = dict(seqtab)

    logger.info("[INFO] collapse_no_mismatch: %d input sequences.",
                len(seq_abundances))

    # Sort by decreasing total abundance
    sorted_seqs = sorted(seq_abundances.keys(),
                         key=lambda s: seq_abundances[s], reverse=True)

    # Build mapping: query seq -> representative seq
    representatives: List[str] = []
    merge_map: Dict[str, str] = {}  # query -> representative

    for query in sorted_seqs:
        added = False
        q_upper = query.upper()
        for idx, ref in enumerate(representatives):
            r_upper = ref.upper()
            if q_upper in r_upper:
                # query is a subsequence of ref — merge under ref (longer)
                merge_map[query] = ref
                added = True
                break
            elif r_upper in q_upper:
                # ref is a subsequence of query — query is longer, so it
                # becomes the new representative and ref merges into it
                merge_map[query] = query
                # Re-point everything that was merged into ref → query
                for k, v in merge_map.items():
                    if v == ref:
                        merge_map[k] = query
                representatives[idx] = query
                added = True
                break
        if not added:
            representatives.append(query)
            merge_map[query] = query

    if is_df:
        # Build collapsed DataFrame
        collapsed = pd.DataFrame(0, index=seqtab.index,
                                 columns=representatives)
        for seq in seqtab.columns:
            rep = merge_map[seq]
            collapsed[rep] = collapsed[rep] + seqtab[seq].values
        # Order by total abundance descending
        col_order = collapsed.sum(axis=0).sort_values(ascending=False).index
        collapsed = collapsed[col_order]
        logger.info("[INFO] collapse_no_mismatch: %d output sequences.",
                    len(collapsed.columns))
        return collapsed
    else:
        out: Dict[str, int] = {}
        for seq, ab in seq_abundances.items():
            rep = merge_map[seq]
            out[rep] = out.get(rep, 0) + int(ab)
        # Sort by decreasing abundance
        out = dict(sorted(out.items(), key=lambda x: x[1], reverse=True))
        logger.info("[INFO] collapse_no_mismatch: %d output sequences.",
                    len(out))
        return out


# ---------------------------------------------------------------------------
# 4. make_sequence_table
# ---------------------------------------------------------------------------

def make_sequence_table(
    samples_dict: Dict[str, Dict[str, int]],
    order_by: str = "abundance",
):
    """Construct a sample-by-sequence observation matrix.

    Args:
        samples_dict: ``{sample_name: {sequence: abundance}}``.
        order_by: How to order the columns.  ``"abundance"`` (default) sorts
            by total abundance across all samples (descending).
            ``"nsamples"`` sorts by number of samples in which the sequence
            appears.  ``None`` for no ordering.

    Returns:
        A pandas DataFrame with samples as rows and sequences as columns.
    """
    import pandas as pd

    # Accept a list of dada results or merger lists (extract denoised dicts)
    if isinstance(samples_dict, list):
        converted = {}
        for i, item in enumerate(samples_dict):
            name = f"sample_{i+1}"
            if isinstance(item, dict) and "denoised" in item:
                converted[name] = item["denoised"]
            elif isinstance(item, dict):
                converted[name] = item
            elif isinstance(item, list):
                # List of merger dicts (from merge_pairs)
                merged = {}
                for m in item:
                    if m.get("accept", True):
                        seq = m["sequence"]
                        merged[seq] = merged.get(seq, 0) + m["abundance"]
                converted[name] = merged
            else:
                raise TypeError(f"Unexpected item type in samples list: {type(item)}")
        samples_dict = converted

    logger.info("[INFO] make_sequence_table: %d samples.", len(samples_dict))

    # Gather all unique sequences
    all_seqs: Dict[str, int] = {}
    for sample_name, uniques in samples_dict.items():
        for seq in uniques:
            all_seqs[seq] = 0

    seq_list = list(all_seqs.keys())
    sample_names = list(samples_dict.keys())

    mat = np.zeros((len(sample_names), len(seq_list)), dtype=int)
    seq_index = {s: i for i, s in enumerate(seq_list)}

    for i, sample_name in enumerate(sample_names):
        for seq, ab in samples_dict[sample_name].items():
            mat[i, seq_index[seq]] = int(ab)

    df = pd.DataFrame(mat, index=sample_names, columns=seq_list)

    # Order columns
    if order_by == "abundance":
        col_sums = df.sum(axis=0)
        df = df[col_sums.sort_values(ascending=False).index]
    elif order_by == "nsamples":
        n_samples = (df > 0).sum(axis=0)
        df = df[n_samples.sort_values(ascending=False).index]

    return df


# ---------------------------------------------------------------------------
# 5. plot_quality_profile
# ---------------------------------------------------------------------------

def plot_quality_profile(
    fastq_files: Union[str, List[str]],
    n: int = 500000,
    output_pdf: Optional[str] = None,
):
    """Plot mean quality score per cycle position from FASTQ file(s).

    Reads up to *n* records from each file, computes per-position quality
    statistics, and produces a matplotlib figure with mean (green), median
    (orange), and 25th/75th percentile (dashed orange) quality lines,
    similar to the R ``plotQualityProfile``.

    Args:
        fastq_files: One or more FASTQ file paths (may be gzipped).
        n: Maximum number of reads to sample per file.
        output_pdf: If provided, save the figure to this PDF path.

    Returns:
        A numpy array of quality scores with shape ``(n_reads, max_len)``
        (values are NaN where reads are shorter than max_len).
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if isinstance(fastq_files, str):
        fastq_files = [fastq_files]

    logger.info("[INFO] plot_quality_profile: processing %d file(s), n=%d.",
                len(fastq_files), n)

    all_quals: List[List[List[int]]] = []  # per-file list of per-read quals
    file_labels: List[str] = []

    for fpath in fastq_files:
        quals_per_read: List[List[int]] = []
        opener = gzip.open if fpath.endswith(".gz") else open
        count = 0
        with opener(fpath, "rt") as fh:
            while count < n:
                header = fh.readline()
                if not header:
                    break
                seq = fh.readline()
                plus = fh.readline()
                qual_line = fh.readline().rstrip("\n\r")
                if not qual_line:
                    break
                quals = [ord(c) - 33 for c in qual_line]
                quals_per_read.append(quals)
                count += 1
        all_quals.append(quals_per_read)
        file_labels.append(os.path.basename(fpath))
        logger.info("[INFO] plot_quality_profile: read %d records from %s.",
                    count, os.path.basename(fpath))

    # Build quality matrix (all files combined for return value)
    max_len = max(
        max((len(q) for q in qlist), default=0)
        for qlist in all_quals
    ) if all_quals else 0
    total_reads = sum(len(qlist) for qlist in all_quals)
    qual_matrix = np.full((total_reads, max_len), np.nan)
    row = 0
    for qlist in all_quals:
        for q in qlist:
            qual_matrix[row, :len(q)] = q
            row += 1

    # Plot
    n_files = len(fastq_files)
    fig, axes = plt.subplots(1, n_files, figsize=(6 * n_files, 4),
                             squeeze=False)

    row_offset = 0
    for fi, qlist in enumerate(all_quals):
        ax = axes[0, fi]
        n_reads = len(qlist)
        if n_reads == 0:
            ax.set_title(file_labels[fi])
            continue

        local_max = max(len(q) for q in qlist)
        local_mat = np.full((n_reads, local_max), np.nan)
        for ri, q in enumerate(qlist):
            local_mat[ri, :len(q)] = q

        positions = np.arange(1, local_max + 1)
        means = np.nanmean(local_mat, axis=0)
        q25 = np.nanpercentile(local_mat, 25, axis=0)
        q50 = np.nanpercentile(local_mat, 50, axis=0)
        q75 = np.nanpercentile(local_mat, 75, axis=0)

        ax.plot(positions, means, color="#66C2A5", linewidth=1.0,
                label="Mean")
        ax.plot(positions, q50, color="#FC8D62", linewidth=0.8,
                label="Median")
        ax.plot(positions, q25, color="#FC8D62", linewidth=0.5,
                linestyle="--", label="Q25")
        ax.plot(positions, q75, color="#FC8D62", linewidth=0.5,
                linestyle="--", label="Q75")
        ax.set_xlabel("Cycle")
        ax.set_ylabel("Quality Score")
        ax.set_title(file_labels[fi])

        read_label = (f"Reads >= {n}" if n_reads >= n
                      else f"Reads: {n_reads}")
        ax.text(0.02, 0.02, read_label, transform=ax.transAxes,
                color="red", fontsize=8, verticalalignment="bottom")
        ax.set_ylim(bottom=0)

        row_offset += n_reads

    plt.tight_layout()
    if output_pdf:
        fig.savefig(output_pdf, format="pdf")
        logger.info("[INFO] plot_quality_profile: saved to %s.", output_pdf)
    plt.close(fig)

    return qual_matrix


# ---------------------------------------------------------------------------
# 6. uniquesto_fasta
# ---------------------------------------------------------------------------

def uniquesto_fasta(
    uniques,
    fasta_path: str,
    ids: Optional[List[str]] = None,
) -> None:
    """Write a uniques dict or list of sequences to a FASTA file.

    If *uniques* is a dict ``{seq: abundance}``, headers are formatted as
    ``>sq1;size=1234;`` (uchime-compatible) unless custom *ids* are given.

    Args:
        uniques: A dict ``{sequence: abundance}`` or any object accepted by
            ``get_uniques``.
        fasta_path: Output FASTA file path.
        ids: Optional custom sequence identifiers.
    """
    uniqs = get_uniques(uniques)
    seqs = list(uniqs.keys())
    abunds = list(uniqs.values())

    if ids is not None and len(ids) != len(seqs):
        raise ValueError(
            f"Length of ids ({len(ids)}) != number of sequences ({len(seqs)})."
        )

    with open(fasta_path, "w") as fh:
        for i, (seq, ab) in enumerate(zip(seqs, abunds)):
            if ids is not None:
                header = ids[i]
            else:
                header = f"sq{i + 1};size={ab};"
            fh.write(f">{header}\n{seq}\n")


# ---------------------------------------------------------------------------
# 7. write_fasta
# ---------------------------------------------------------------------------

def write_fasta(
    seqs: Sequence[str],
    fasta_path: str,
    ids: Optional[List[str]] = None,
) -> None:
    """Write a list of sequences to a FASTA file.

    Args:
        seqs: Iterable of DNA sequence strings.
        fasta_path: Output file path.
        ids: Optional list of identifiers (one per sequence).
            Defaults to ``1, 2, 3, ...``.
    """
    seqs = list(seqs)
    if ids is not None and len(ids) != len(seqs):
        raise ValueError(
            f"Length of ids ({len(ids)}) != number of sequences ({len(seqs)})."
        )
    with open(fasta_path, "w") as fh:
        for i, seq in enumerate(seqs):
            header = ids[i] if ids is not None else str(i + 1)
            fh.write(f">{header}\n{seq}\n")


# ---------------------------------------------------------------------------
# 8. is_phix
# ---------------------------------------------------------------------------

def is_phix(
    seqs,
    ref_path: Optional[str] = None,
    word_size: int = 16,
    min_matches: int = 2,
) -> np.ndarray:
    """Check sequences against the PhiX genome using kmer matching.

    For each query sequence, kmers of size ``word_size`` are extracted and
    compared to the PhiX reference (forward and reverse complement).  A
    sequence is flagged as PhiX if at least ``min_matches`` kmers hit.

    Args:
        seqs: Sequences to check (any type accepted by ``get_sequences``).
        ref_path: Path to a FASTA file containing the PhiX genome.  If None,
            a built-in 100 bp snippet is used for matching.
        word_size: Kmer size for matching.
        min_matches: Minimum number of kmer hits to call PhiX.

    Returns:
        A boolean numpy array, True where a sequence matches PhiX.
    """
    seqs = get_sequences(seqs)

    if ref_path is not None:
        recs = _parse_fasta(ref_path)
        phix_seq = "".join(seq for _, seq in recs).upper()
    else:
        phix_seq = _PHIX_100BP.upper()

    phix_rc = _rc(phix_seq)

    # Build kmer sets for forward and RC PhiX
    def _kmer_set(s, k):
        return set(s[i:i+k] for i in range(len(s) - k + 1))

    phix_kmers = _kmer_set(phix_seq, word_size)
    phix_rc_kmers = _kmer_set(phix_rc, word_size)
    all_phix_kmers = phix_kmers | phix_rc_kmers

    result = np.zeros(len(seqs), dtype=bool)
    for i, seq in enumerate(seqs):
        seq = seq.upper()
        hits = 0
        for j in range(len(seq) - word_size + 1):
            kmer = seq[j:j + word_size]
            if kmer in all_phix_kmers:
                hits += 1
                if hits >= min_matches:
                    break
        result[i] = hits >= min_matches

    return result


# ---------------------------------------------------------------------------
# 9. seq_complexity
# ---------------------------------------------------------------------------

def seq_complexity(
    seqs,
    kmer_size: int = 2,
    window: Optional[int] = None,
    by: int = 5,
) -> np.ndarray:
    """Calculate sequence complexity as Shannon effective number of kmers.

    Complexity is the exponential of the Shannon entropy of kmer frequencies.
    A perfectly random sequence of sufficient length will approach
    ``4**kmer_size``.  Repetitive / low-complexity sequences will have values
    well below this maximum.

    If a *window* is provided, the minimum complexity observed over a sliding
    window along each sequence is returned.

    Args:
        seqs: Sequences (any type accepted by ``get_sequences``).
        kmer_size: Size of kmers to count.  Default 2 (dinucleotides).
        window: Width of sliding window in nucleotides.  If None, the
            whole sequence is used.
        by: Step size for the sliding window.

    Returns:
        A numpy array of complexity values, one per input sequence.
    """
    seqs = get_sequences(seqs)
    si_max = 4 ** kmer_size

    if window is not None and kmer_size >= window:
        raise ValueError("window must be larger than kmer_size.")

    def _kmer_complexity(s: str) -> float:
        """Shannon effective kmer number for a single string."""
        k = kmer_size
        if len(s) < k:
            return 0.0
        counts: Dict[str, int] = {}
        for i in range(len(s) - k + 1):
            kmer = s[i:i+k]
            # Ignore kmers with non-ACGT
            if all(c in "ACGT" for c in kmer):
                counts[kmer] = counts.get(kmer, 0) + 1
        if not counts:
            return 0.0
        return _sindex(list(counts.values()))

    result = np.zeros(len(seqs), dtype=float)
    for i, seq in enumerate(seqs):
        seq = seq.upper()
        if window is None:
            result[i] = _kmer_complexity(seq)
        else:
            min_si = si_max
            seq_len = len(seq)
            if seq_len < window:
                result[i] = _kmer_complexity(seq)
                continue
            for start in range(0, seq_len - window + 1, by):
                subseq = seq[start:start + window]
                si = _kmer_complexity(subseq)
                if si < min_si:
                    min_si = si
            result[i] = min_si

    return result


def plot_sankey(
    track: dict,
    title: str = "Read tracking through papa2 pipeline",
    output: Optional[str] = None,
    width: int = 900,
    height: int = 500,
):
    """Create a Sankey diagram showing read flow through pipeline stages.

    Args:
        track: Dict mapping stage names to read/sequence counts.
            Typical keys (in order):
            ``{"input": 50000, "filtered": 45000, "denoised": 42000,
               "merged": 40000, "non-chimeric": 38000}``

            Can also be a list of such dicts (one per sample) — values
            will be summed across samples.

            Or a pandas DataFrame with stage columns and sample rows
            (as produced by a read-tracking table).

        title: Plot title.
        output: If given, save to this path (.html for interactive,
            .png/.svg/.pdf require kaleido). If None, returns the
            plotly Figure object.
        width: Figure width in pixels.
        height: Figure height in pixels.

    Returns:
        plotly.graph_objects.Figure if output is None, else None.
    """
    try:
        import plotly.graph_objects as go
    except ImportError:
        raise ImportError("plotly is required for plot_sankey: pip install plotly")

    # Normalize input
    if hasattr(track, "iterrows"):  # DataFrame
        # Sum columns to get totals per stage
        stages = list(track.columns)
        counts = [int(track[col].sum()) for col in stages]
    elif isinstance(track, list):
        # List of dicts — sum values per key
        all_keys = list(track[0].keys())
        counts = [sum(d.get(k, 0) for d in track) for k in all_keys]
        stages = all_keys
    elif isinstance(track, dict):
        stages = list(track.keys())
        counts = [int(track[k]) for k in stages]
    else:
        raise TypeError(f"track must be a dict, list of dicts, or DataFrame, got {type(track)}")

    if len(stages) < 2:
        raise ValueError("Need at least 2 stages for a Sankey diagram")

    # Build Sankey data
    # Nodes: each stage + a "lost" node for each transition
    node_labels = list(stages)
    node_colors = []

    # Color palette: blues for pipeline stages, greens for taxonomy ranks, red for lost
    stage_colors = [
        "#2196F3", "#1E88E5", "#1976D2", "#1565C0",
        "#0D47A1", "#0A3D91",
        # taxonomy ranks (greens)
        "#4CAF50", "#43A047", "#388E3C", "#2E7D32",
        "#1B5E20", "#145214", "#0E3E0E",
    ]
    lost_color = "#EF5350"

    for i in range(len(stages)):
        node_colors.append(stage_colors[i % len(stage_colors)])

    # Add "lost" nodes
    for i in range(len(stages) - 1):
        lost = counts[i] - counts[i + 1]
        if lost > 0:
            node_labels.append(f"lost ({stages[i]}→{stages[i+1]})")
            node_colors.append(lost_color)

    source = []
    target = []
    value = []
    link_colors = []

    lost_idx = len(stages)  # first "lost" node index
    for i in range(len(stages) - 1):
        retained = counts[i + 1]
        lost = counts[i] - counts[i + 1]

        # Retained flow: stage[i] → stage[i+1]
        if retained > 0:
            source.append(i)
            target.append(i + 1)
            value.append(retained)
            link_colors.append("rgba(33, 150, 243, 0.4)")

        # Lost flow: stage[i] → lost node
        if lost > 0:
            source.append(i)
            target.append(lost_idx)
            value.append(lost)
            link_colors.append("rgba(239, 83, 80, 0.4)")
            lost_idx += 1

    # Format counts for labels
    def fmt(n):
        if n >= 1_000_000:
            return f"{n/1_000_000:.1f}M"
        if n >= 1_000:
            return f"{n/1_000:.1f}K"
        return str(n)

    node_labels_fmt = []
    for i, label in enumerate(node_labels):
        if i < len(counts):
            node_labels_fmt.append(f"{label}<br>{fmt(counts[i])}")
        else:
            # Lost node — find its value
            lost_val = 0
            for j, s in enumerate(source):
                if target[j] == i:
                    lost_val = value[j]
                    break
            node_labels_fmt.append(f"{label}<br>{fmt(lost_val)}")

    fig = go.Figure(data=[go.Sankey(
        node=dict(
            pad=20,
            thickness=30,
            line=dict(color="black", width=0.5),
            label=node_labels_fmt,
            color=node_colors,
        ),
        link=dict(
            source=source,
            target=target,
            value=value,
            color=link_colors,
        ),
    )])

    fig.update_layout(
        title_text=title,
        font_size=13,
        width=width,
        height=height,
    )

    if output:
        if output.endswith(".html"):
            fig.write_html(output)
        else:
            fig.write_image(output)
        return None

    return fig


def _seqtab_sum(st):
    """Return total read count from a sequence table (array, dict, or DataFrame)."""
    if hasattr(st, "values"):  # DataFrame
        return int(st.values.sum())
    if isinstance(st, dict) and "table" in st:
        return int(st["table"].sum())
    if hasattr(st, "sum"):  # numpy array
        return int(st.sum())
    return 0


def _seqtab_seqs(st):
    """Return the list of ASV sequences from a sequence table."""
    if isinstance(st, dict) and "seqs" in st:
        return st["seqs"]
    if hasattr(st, "columns"):  # DataFrame with seq columns
        return list(st.columns)
    return None


def track_reads(
    dereps=None,
    dadas=None,
    mergers=None,
    seqtab=None,
    seqtab_nochim=None,
    taxa=None,
):
    """Build a read-tracking dict from pipeline stage outputs.

    Pass whichever stages you have — earlier stages are required for
    later ones to make sense, but all are optional.

    Args:
        dereps: List of derep dicts (from derep_fastq)
        dadas: List of dada result dicts
        mergers: List of merger lists (from merge_pairs)
        seqtab: Sequence table (dict with 'table'/'seqs', DataFrame, or array)
        seqtab_nochim: Chimera-filtered sequence table (same formats as seqtab)
        taxa: Taxonomy array from assign_species / add_species, shape (N, K).
            Used together with seqtab_nochim (or seqtab) to count reads
            assigned to classified ASVs.

    Returns:
        Dict mapping stage names to total read counts.
        Pass directly to plot_sankey().
    """
    track = {}
    if dereps is not None:
        track["input"] = int(sum(d["abundances"].sum() for d in dereps))
    if dadas is not None:
        if isinstance(dadas, dict):
            dadas = [dadas]
        track["denoised"] = int(sum(sum(d["denoised"].values()) for d in dadas))
    if mergers is not None:
        total = 0
        for m in mergers:
            if isinstance(m, list):
                total += sum(x["abundance"] for x in m if x.get("accept", True))
            elif isinstance(m, dict) and "abundance" in m:
                if m.get("accept", True):
                    total += m["abundance"]
        track["merged"] = int(total)
    if seqtab is not None:
        track["tabled"] = _seqtab_sum(seqtab)
    if seqtab_nochim is not None:
        track["non-chimeric"] = _seqtab_sum(seqtab_nochim)
    if taxa is not None:
        # Break down reads by taxonomic rank.
        # taxa can be:
        #   - pandas DataFrame with rank columns (Kingdom, Phylum, ..., Species)
        #   - numpy array of shape (N_asvs, K_ranks)
        #   - numpy array of shape (N_asvs,) for single-rank
        import numpy as np

        # Get per-ASV read totals from the best available table
        ref_st = seqtab_nochim if seqtab_nochim is not None else seqtab
        if ref_st is not None:
            if isinstance(ref_st, dict) and "table" in ref_st:
                tbl = ref_st["table"]
            elif hasattr(ref_st, "values"):
                tbl = ref_st.values
            else:
                tbl = np.asarray(ref_st)
            per_asv = tbl.sum(axis=0) if tbl.ndim == 2 else tbl

            # Determine rank names and values
            if hasattr(taxa, "columns"):
                # pandas DataFrame — use column names
                rank_names = list(taxa.columns)
                taxa_vals = taxa.values
            else:
                taxa_vals = np.asarray(taxa)
                if taxa_vals.ndim == 1:
                    taxa_vals = taxa_vals.reshape(-1, 1)
                rank_names = [
                    "Kingdom", "Phylum", "Class", "Order",
                    "Family", "Genus", "Species",
                ][:taxa_vals.shape[1]]

            if len(per_asv) == taxa_vals.shape[0]:
                for k, rank in enumerate(rank_names):
                    mask = np.array(
                        [taxa_vals[i, k] is not None
                         and str(taxa_vals[i, k]).strip() != ""
                         for i in range(taxa_vals.shape[0])]
                    )
                    track[rank.lower()] = int(per_asv[mask].sum())
    return track


# ---------------------------------------------------------------------------
# get_errors -- extract error information from various dada2 objects
# ---------------------------------------------------------------------------

def get_errors(
    obj,
    detailed: bool = False,
    enforce: bool = True,
):
    """Extract error rate information from various dada2 object types.

    Mirrors R's ``getErrors``.

    Supported inputs:
        - numpy array: used directly as ``err_out``.
        - dada result dict: extracts ``err_out``, ``err_in``, ``trans``.
        - list of dada result dicts: verifies all share the same ``err_out``,
          accumulates ``trans`` across samples.

    Args:
        obj: Input object containing error information.
        detailed: If True, return a dict with ``err_out``, ``err_in``, and
            ``trans`` keys.  If False (default), return just the ``err_out``
            matrix.
        enforce: If True, validate that the error matrix has 16 rows,
            is numeric, and all values are in [0, 1].

    Returns:
        numpy array (16, ncol) if ``detailed=False``, or a dict with
        ``err_out``, ``err_in``, ``trans`` keys if ``detailed=True``.
    """
    err_out = None
    err_in = None
    trans = None

    if isinstance(obj, np.ndarray):
        err_out = obj

    elif isinstance(obj, dict):
        err_out = obj.get("err_out", None)
        err_in = obj.get("err_in", None)
        trans = obj.get("trans", None)
        # If the dict doesn't have err_out but has an error matrix shape
        if err_out is None and "err" in obj:
            err_out = obj["err"]

    elif isinstance(obj, list):
        # List of dada result dicts
        if len(obj) == 0:
            raise ValueError("Empty list provided to get_errors.")

        # Use the first element's err_out as reference
        first = obj[0]
        if isinstance(first, dict):
            err_out = first.get("err_out", None)
            err_in = first.get("err_in", None)
            trans = first.get("trans", None)
            if trans is not None:
                trans = trans.copy().astype(np.float64)

            for item in obj[1:]:
                if not isinstance(item, dict):
                    raise TypeError(
                        f"Expected dict in list, got {type(item).__name__}."
                    )
                item_err = item.get("err_out", None)
                if item_err is not None and err_out is not None:
                    if not np.array_equal(item_err, err_out):
                        raise ValueError(
                            "Not all dada results share the same err_out."
                        )
                item_trans = item.get("trans", None)
                if item_trans is not None:
                    if trans is None:
                        trans = item_trans.copy().astype(np.float64)
                    else:
                        trans = trans + item_trans.astype(np.float64)
        else:
            raise TypeError(
                f"Expected list of dicts, got list of {type(first).__name__}."
            )
    else:
        raise TypeError(
            f"Cannot extract errors from object of type "
            f"{type(obj).__name__}."
        )

    if err_out is None:
        raise ValueError("Could not extract error matrix from input.")

    err_out = np.asarray(err_out, dtype=np.float64)

    if enforce:
        if err_out.ndim != 2 or err_out.shape[0] != 16:
            raise ValueError(
                f"Error matrix must have 16 rows, got shape {err_out.shape}."
            )
        if not np.issubdtype(err_out.dtype, np.number):
            raise ValueError("Error matrix must be numeric.")
        if np.any(err_out < 0.0) or np.any(err_out > 1.0):
            raise ValueError(
                "Error matrix values must be in [0, 1]."
            )

    if detailed:
        return {
            "err_out": err_out,
            "err_in": np.asarray(err_in, dtype=np.float64) if err_in is not None else None,
            "trans": np.asarray(trans, dtype=np.float64) if trans is not None else None,
        }
    return err_out


# ---------------------------------------------------------------------------
# merge_sequence_tables -- merge multiple sequence tables
# ---------------------------------------------------------------------------

def merge_sequence_tables(
    *tables,
    order_by: str = "abundance",
    try_rc: bool = False,
):
    """Merge multiple sequence tables into one.

    Mirrors R's ``mergeSequenceTables``.

    Each table is a dict with:
        - ``'table'``: numpy array (nsamples x nseqs)
        - ``'seqs'``: list of sequence strings

    Sequences present in multiple tables have their counts summed.
    If ``try_rc`` is True, reverse-complement sequences are detected
    and re-oriented to match the majority orientation.

    Args:
        *tables: Variable number of sequence table dicts.
        order_by: Column ordering: ``"abundance"`` (total counts, descending)
            or ``"nsamples"`` (number of samples present, descending).
        try_rc: If True, check for reverse-complement sequences across
            tables and re-orient them before merging.

    Returns:
        A dict with ``'table'`` (numpy array) and ``'seqs'`` (list of str).
    """
    if len(tables) == 0:
        raise ValueError("No tables provided to merge.")

    if len(tables) == 1:
        return tables[0]

    # Collect all unique sequences and their sample counts
    # Each table contributes rows (samples) with counts for its sequences
    all_seqs_set: Dict[str, int] = {}  # seq -> index in merged
    sample_rows: List[np.ndarray] = []  # list of arrays, one per sample
    seq_lists: List[List[str]] = []

    # If try_rc, build a mapping of RC sequences
    rc_map: Dict[str, str] = {}  # maps RC form -> canonical form

    # First pass: gather all sequences
    all_seqs_ordered: List[str] = []
    for tbl in tables:
        tseqs = tbl["seqs"]
        for seq in tseqs:
            canon = seq
            if try_rc and seq not in all_seqs_set:
                seq_rc = _rc(seq)
                if seq_rc in all_seqs_set:
                    # This sequence's RC is already known; map it
                    rc_map[seq] = seq_rc
                    canon = seq_rc
            if canon not in all_seqs_set:
                all_seqs_set[canon] = len(all_seqs_ordered)
                all_seqs_ordered.append(canon)

    n_seqs = len(all_seqs_ordered)
    seq_index = {s: i for i, s in enumerate(all_seqs_ordered)}

    # Second pass: build merged count matrix
    all_sample_counts: List[np.ndarray] = []

    for tbl in tables:
        tseqs = tbl["seqs"]
        tmat = np.asarray(tbl["table"])
        if tmat.ndim == 1:
            tmat = tmat.reshape(1, -1)
        nsamples = tmat.shape[0]

        for si in range(nsamples):
            row = np.zeros(n_seqs, dtype=np.int64)
            for ji, seq in enumerate(tseqs):
                canon = rc_map.get(seq, seq)
                idx = seq_index[canon]
                row[idx] += tmat[si, ji]
            all_sample_counts.append(row)

    merged_mat = np.vstack(all_sample_counts)

    # Order columns
    if order_by == "abundance":
        col_totals = merged_mat.sum(axis=0)
        order = np.argsort(-col_totals)
    elif order_by == "nsamples":
        n_present = (merged_mat > 0).sum(axis=0)
        order = np.argsort(-n_present)
    else:
        order = np.arange(n_seqs)

    merged_mat = merged_mat[:, order]
    merged_seqs = [all_seqs_ordered[i] for i in order]

    return {"table": merged_mat, "seqs": merged_seqs}


# ---------------------------------------------------------------------------
# nwhamming -- Hamming distance via NW alignment
# ---------------------------------------------------------------------------

def nwhamming(s1, s2, **kwargs):
    """Compute Hamming distance between sequences via NW alignment.

    Aligns the two sequences using ``nwalign`` and then counts mismatches
    and indels using ``eval_pair``.

    Vectorized: if *s1* and *s2* are both lists (of the same length),
    returns a list of distances.

    Args:
        s1: DNA sequence string, or list of strings.
        s2: DNA sequence string, or list of strings.
        **kwargs: Additional keyword arguments passed to ``nwalign``
            (e.g. ``match``, ``mismatch``, ``gap_p``, ``band``).

    Returns:
        int (scalar inputs) or list of int (list inputs): mismatch + indel
        count for each pair.
    """
    from ._cdada import nwalign as _nwalign, eval_pair as _eval_pair

    is_list = isinstance(s1, (list, tuple))
    if is_list:
        if not isinstance(s2, (list, tuple)) or len(s1) != len(s2):
            raise ValueError(
                "s1 and s2 must be lists of the same length for vectorized nwhamming."
            )
        results = []
        for a, b in zip(s1, s2):
            al1, al2 = _nwalign(a, b, **kwargs)
            nmatch, nmismatch, nindel = _eval_pair(al1, al2)
            results.append(nmismatch + nindel)
        return results

    al1, al2 = _nwalign(s1, s2, **kwargs)
    nmatch, nmismatch, nindel = _eval_pair(al1, al2)
    return nmismatch + nindel


# ---------------------------------------------------------------------------
# is_shift_denovo -- detect shifted sequences
# ---------------------------------------------------------------------------

def is_shift_denovo(
    unqs,
    min_overlap: int = 20,
    verbose: bool = False,
) -> np.ndarray:
    """Check if sequences are shifted versions of more-abundant sequences.

    For each sequence (sorted by decreasing abundance), check whether it
    is a "shift" of any more-abundant sequence.  A shifted pair has:

    - ``match < len(sq1)`` AND ``match < len(sq2)``
    - ``match >= min_overlap``
    - ``mismatch == 0`` AND ``indel == 0``

    This identifies reads that are identical subsequences but offset
    (shifted) relative to each other.

    Args:
        unqs: dict ``{sequence: abundance}`` (as from ``get_uniques``),
            or a list of sequences sorted by decreasing abundance.
        min_overlap: Minimum overlap (match) length to call a shift.
        verbose: If True, log details about detected shifts.

    Returns:
        Boolean numpy array, True where a sequence is a shifted duplicate
        of a more-abundant sequence.
    """
    from ._cdada import nwalign as _nwalign, eval_pair as _eval_pair

    if isinstance(unqs, dict):
        # Sort by decreasing abundance
        sorted_items = sorted(unqs.items(), key=lambda x: x[1], reverse=True)
        seqs = [item[0] for item in sorted_items]
    elif isinstance(unqs, (list, tuple)):
        seqs = list(unqs)
    else:
        seqs = list(unqs)

    n = len(seqs)
    is_shift = np.zeros(n, dtype=bool)

    for i in range(1, n):
        sq_i = seqs[i]
        len_i = len(sq_i)
        for j in range(i):
            if is_shift[j]:
                continue
            sq_j = seqs[j]
            len_j = len(sq_j)

            al1, al2 = _nwalign(sq_j, sq_i)
            nmatch, nmismatch, nindel = _eval_pair(al1, al2)

            if (nmatch < len_j and nmatch < len_i
                    and nmatch >= min_overlap
                    and nmismatch == 0 and nindel == 0):
                is_shift[i] = True
                if verbose:
                    logger.info(
                        "[INFO] is_shift_denovo: seq %d is a shift of seq %d "
                        "(overlap=%d)", i, j, nmatch
                    )
                break

    return is_shift


# ---------------------------------------------------------------------------
# plot_errors -- error rate diagnostic plot
# ---------------------------------------------------------------------------

def plot_errors(
    dq,
    nti: Tuple[str, ...] = ("A", "C", "G", "T"),
    ntj: Tuple[str, ...] = ("A", "C", "G", "T"),
    obs: bool = True,
    err_out: bool = True,
    err_in: bool = False,
    nominal_q: bool = False,
    output: Optional[str] = None,
):
    """Error rate diagnostic plot using plotly.

    Mirrors R's ``plotErrors``.

    Creates a faceted plot showing error rates for each nucleotide
    transition (from_nuc -> to_nuc).  Self-transitions are blanked out.

    Args:
        dq: A dada result dict (with ``err_out``, optionally ``err_in``
            and ``trans``), or a numpy error matrix (16, ncol).
        nti: Tuple of source nucleotides to include.
        ntj: Tuple of target nucleotides to include.
        obs: If True, show observed error rates as scatter points.
        err_out: If True, show the estimated (output) error rates as a line.
        err_in: If True, show the input error rates as a dashed line.
        nominal_q: If True, show nominal Q-score error rates as a red line.
        output: If given, save to this path (.html for interactive,
            .png/.svg/.pdf require kaleido).  If None, returns the
            plotly Figure object.

    Returns:
        plotly.graph_objects.Figure if output is None, else None.
    """
    try:
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
    except ImportError:
        raise ImportError("plotly is required for plot_errors: pip install plotly")

    nt_map = {"A": 0, "C": 1, "G": 2, "T": 3}

    # Extract error info
    err_info = get_errors(dq, detailed=True, enforce=False)
    err_out_mat = err_info["err_out"]
    err_in_mat = err_info.get("err_in", None) if err_in else None
    trans_mat = err_info.get("trans", None) if obs else None

    ncol = err_out_mat.shape[1]
    qq = np.arange(ncol)

    nti_indices = [nt_map[nt] for nt in nti]
    ntj_indices = [nt_map[nt] for nt in ntj]

    n_rows = len(nti_indices)
    n_cols = len(ntj_indices)

    subplot_titles = []
    for ni in nti:
        for nj in ntj:
            subplot_titles.append(f"{ni} -> {nj}" if ni != nj else "")

    fig = make_subplots(
        rows=n_rows, cols=n_cols,
        subplot_titles=subplot_titles,
        shared_xaxes=True, shared_yaxes=True,
        horizontal_spacing=0.03, vertical_spacing=0.05,
    )

    for ri, ni_idx in enumerate(nti_indices):
        for ci, nj_idx in enumerate(ntj_indices):
            row_plot = ri + 1
            col_plot = ci + 1
            row_idx = ni_idx * 4 + nj_idx

            if ni_idx == nj_idx:
                # Self-transition: blank
                continue

            show_legend = (ri == 0 and ci == 1)

            # Observed error rates (scatter)
            if obs and trans_mat is not None:
                from .error import _BASE_ROWS
                base_rows = _BASE_ROWS[ni_idx]
                tot = trans_mat[base_rows].sum(axis=0).astype(np.float64)
                errs = trans_mat[row_idx].astype(np.float64)
                with np.errstate(divide="ignore", invalid="ignore"):
                    obs_rate = errs / tot
                obs_rate[~np.isfinite(obs_rate)] = np.nan

                # Point sizes proportional to sqrt(total counts)
                sizes = np.sqrt(np.maximum(tot, 0)) / np.sqrt(np.nanmax(tot) + 1) * 12 + 2

                fig.add_trace(
                    go.Scatter(
                        x=qq, y=obs_rate,
                        mode="markers",
                        marker=dict(
                            size=sizes,
                            color="rgba(0, 0, 0, 0.4)",
                        ),
                        name="Observed",
                        showlegend=show_legend,
                        legendgroup="obs",
                    ),
                    row=row_plot, col=col_plot,
                )

            # Estimated error rates (line)
            if err_out:
                fig.add_trace(
                    go.Scatter(
                        x=qq, y=err_out_mat[row_idx],
                        mode="lines",
                        line=dict(color="#1f77b4", width=2),
                        name="Estimated",
                        showlegend=show_legend,
                        legendgroup="est",
                    ),
                    row=row_plot, col=col_plot,
                )

            # Input error rates (dashed line)
            if err_in and err_in_mat is not None:
                fig.add_trace(
                    go.Scatter(
                        x=qq, y=err_in_mat[row_idx],
                        mode="lines",
                        line=dict(color="#ff7f0e", width=1.5, dash="dash"),
                        name="Input",
                        showlegend=show_legend,
                        legendgroup="in",
                    ),
                    row=row_plot, col=col_plot,
                )

            # Nominal Q-score error rates (red line)
            if nominal_q:
                nominal_rates = 10.0 ** (-qq / 10.0)
                fig.add_trace(
                    go.Scatter(
                        x=qq, y=nominal_rates,
                        mode="lines",
                        line=dict(color="red", width=1, dash="dot"),
                        name="Nominal Q",
                        showlegend=show_legend,
                        legendgroup="nom",
                    ),
                    row=row_plot, col=col_plot,
                )

            # Log scale for y-axis
            fig.update_yaxes(type="log", row=row_plot, col=col_plot)

    fig.update_layout(
        title="Error rates by nucleotide transition",
        height=200 * n_rows + 100,
        width=200 * n_cols + 100,
    )
    fig.update_xaxes(title_text="Quality Score")
    fig.update_yaxes(title_text="Error Rate")

    if output:
        if output.endswith(".html"):
            fig.write_html(output)
        else:
            fig.write_image(output)
        return None

    return fig


# ---------------------------------------------------------------------------
# plot_complexity -- sequence complexity histogram
# ---------------------------------------------------------------------------

def plot_complexity(
    files: Union[str, List[str]],
    kmer_size: int = 2,
    n: int = 100000,
    bins: int = 100,
    output: Optional[str] = None,
):
    """Sequence complexity histogram using plotly.

    Samples *n* reads from each FASTQ file, computes sequence complexity
    (Shannon effective number of kmers), and plots a faceted histogram.

    Args:
        files: One or more FASTQ file paths (may be gzipped).
        kmer_size: Kmer size for complexity calculation (default 2).
        n: Maximum number of reads to sample per file.
        bins: Number of histogram bins.
        output: If given, save to this path (.html for interactive,
            .png/.svg/.pdf require kaleido).  If None, returns the
            plotly Figure object.

    Returns:
        plotly.graph_objects.Figure if output is None, else None.
    """
    try:
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
    except ImportError:
        raise ImportError("plotly is required for plot_complexity: pip install plotly")

    if isinstance(files, str):
        files = [files]

    n_files = len(files)
    fig = make_subplots(
        rows=n_files, cols=1,
        subplot_titles=[os.path.basename(f) for f in files],
        shared_xaxes=True,
        vertical_spacing=0.05 if n_files > 1 else 0.0,
    )

    for fi, fpath in enumerate(files):
        # Read sequences from FASTQ
        seqs_read: List[str] = []
        opener = gzip.open if fpath.endswith(".gz") else open
        count = 0
        with opener(fpath, "rt") as fh:
            while count < n:
                header = fh.readline()
                if not header:
                    break
                seq = fh.readline().rstrip("\n\r")
                plus = fh.readline()
                qual = fh.readline()
                if not seq:
                    break
                seqs_read.append(seq)
                count += 1

        # Compute complexity
        complexities = seq_complexity(seqs_read, kmer_size=kmer_size)

        fig.add_trace(
            go.Histogram(
                x=complexities,
                nbinsx=bins,
                name=os.path.basename(fpath),
                showlegend=True,
            ),
            row=fi + 1, col=1,
        )

    fig.update_layout(
        title="Sequence Complexity",
        height=300 * n_files + 100,
        width=700,
    )
    fig.update_xaxes(title_text="Complexity (Shannon effective kmers)")
    fig.update_yaxes(title_text="Count")

    if output:
        if output.endswith(".html"):
            fig.write_html(output)
        else:
            fig.write_image(output)
        return None

    return fig


# ---------------------------------------------------------------------------
# remove_primers -- match and trim primer sequences from reads
# ---------------------------------------------------------------------------

def _hamming_match(seq_region: str, primer: str, max_mismatch: int) -> bool:
    """Check if seq_region matches primer with at most max_mismatch mismatches.

    Both strings must be the same length.
    """
    if len(seq_region) != len(primer):
        return False
    mismatches = 0
    for a, b in zip(seq_region, primer):
        if a != b:
            mismatches += 1
            if mismatches > max_mismatch:
                return False
    return True


def remove_primers(
    fn: str,
    fout: str,
    primer_fwd: str,
    primer_rev: Optional[str] = None,
    max_mismatch: int = 2,
    trim_fwd: bool = True,
    trim_rev: bool = True,
    orient: bool = True,
    compress: bool = True,
    verbose: bool = False,
) -> Tuple[int, int]:
    """Match and trim primer sequences from FASTQ reads.

    Reads the input FASTQ file, matches the forward primer at the start
    of each read and the reverse-complement of the reverse primer at
    the end.  Reads that match are trimmed (if ``trim_fwd`` / ``trim_rev``)
    and written to the output file.

    If ``orient=True``, reads that don't match in forward orientation
    are checked in reverse-complement and flipped if they match.

    Args:
        fn: Input FASTQ file path (may be gzipped).
        fout: Output FASTQ file path.
        primer_fwd: Forward primer sequence (5' to 3').
        primer_rev: Reverse primer sequence (5' to 3'), optional.
            Its reverse complement is matched at the 3' end of reads.
        max_mismatch: Maximum allowed mismatches per primer match.
        trim_fwd: If True, trim the matched forward primer region.
        trim_rev: If True, trim the matched reverse primer region.
        orient: If True, check reverse complement of reads and flip
            if primers match in that orientation.
        compress: If True and output path ends with ``.gz``, gzip-compress
            the output.
        verbose: If True, log progress and match statistics.

    Returns:
        Tuple ``(reads_in, reads_out)`` with the number of reads read
        and the number of reads written.
    """
    primer_fwd = primer_fwd.upper()
    len_fwd = len(primer_fwd)

    if primer_rev is not None:
        primer_rev = primer_rev.upper()
        primer_rev_rc = _rc(primer_rev)
        len_rev_rc = len(primer_rev_rc)
    else:
        primer_rev_rc = None
        len_rev_rc = 0

    reads_in = 0
    reads_out = 0

    opener_in = gzip.open if fn.endswith(".gz") else open
    if compress and fout.endswith(".gz"):
        opener_out = gzip.open
        mode_out = "wt"
    else:
        opener_out = open
        mode_out = "w"

    with opener_in(fn, "rt") as fh_in, opener_out(fout, mode_out) as fh_out:
        while True:
            header = fh_in.readline()
            if not header:
                break
            seq = fh_in.readline().rstrip("\n\r")
            plus = fh_in.readline()
            qual = fh_in.readline().rstrip("\n\r")
            if not seq:
                break

            reads_in += 1
            seq_upper = seq.upper()
            matched = False
            is_rc = False

            # Try forward orientation
            fwd_ok = False
            rev_ok = False

            if len(seq_upper) >= len_fwd:
                fwd_ok = _hamming_match(seq_upper[:len_fwd], primer_fwd, max_mismatch)

            if primer_rev_rc is not None and len(seq_upper) >= len_rev_rc:
                rev_ok = _hamming_match(
                    seq_upper[len(seq_upper) - len_rev_rc:],
                    primer_rev_rc, max_mismatch,
                )
            elif primer_rev_rc is None:
                rev_ok = True  # No reverse primer to check

            if fwd_ok and rev_ok:
                matched = True
            elif orient and not matched:
                # Try reverse complement
                seq_rc = _rc(seq_upper)
                qual_rc = qual[::-1]

                fwd_ok_rc = False
                rev_ok_rc = False

                if len(seq_rc) >= len_fwd:
                    fwd_ok_rc = _hamming_match(seq_rc[:len_fwd], primer_fwd, max_mismatch)

                if primer_rev_rc is not None and len(seq_rc) >= len_rev_rc:
                    rev_ok_rc = _hamming_match(
                        seq_rc[len(seq_rc) - len_rev_rc:],
                        primer_rev_rc, max_mismatch,
                    )
                elif primer_rev_rc is None:
                    rev_ok_rc = True

                if fwd_ok_rc and rev_ok_rc:
                    matched = True
                    is_rc = True
                    seq = seq_rc
                    qual = qual_rc

            if not matched:
                continue

            # Trim primers
            start = len_fwd if (trim_fwd and fwd_ok) or (trim_fwd and is_rc and fwd_ok_rc) else 0
            end = len(seq) - len_rev_rc if (trim_rev and primer_rev_rc is not None) else len(seq)
            if end <= start:
                continue

            seq_out = seq[start:end]
            qual_out = qual[start:end]

            fh_out.write(header)
            fh_out.write(seq_out + "\n")
            fh_out.write("+\n")
            fh_out.write(qual_out + "\n")
            reads_out += 1

    if verbose:
        logger.info(
            "[INFO] remove_primers: %d reads in, %d reads out (%.1f%% matched).",
            reads_in, reads_out,
            100.0 * reads_out / reads_in if reads_in > 0 else 0.0,
        )

    return reads_in, reads_out
