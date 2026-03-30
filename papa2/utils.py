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

    for i, query in enumerate(seqs):
        query_up = query.upper()
        # Find all reference sequences that contain the query as a substring
        hit_indices = []
        for j, rseq in enumerate(ref_seqs):
            if query_up in rseq:
                hit_indices.append(j)
        # Also check reverse complement
        if try_rc and not hit_indices:
            query_rc = _rc(query_up)
            for j, rseq in enumerate(ref_seqs):
                if query_rc in rseq:
                    hit_indices.append(j)

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
        for ref in representatives:
            r_upper = ref.upper()
            # Check if one is a substring of the other (no-mismatch collapse)
            if q_upper in r_upper or r_upper in q_upper:
                merge_map[query] = ref
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
