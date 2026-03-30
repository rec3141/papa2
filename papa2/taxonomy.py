"""Taxonomy assignment for amplicon sequences using the RDP Naive Bayesian Classifier.

This module wraps the C ``run_taxonomy`` function from ``papa2._cdada`` to
classify amplicon sequence variants (ASVs) against a reference taxonomy
database.  It faithfully ports the logic of dada2's ``assignTaxonomy()``
R function, including UNITE fungal database support.
"""

import logging
import re
from typing import Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd

from .utils import _parse_fasta, _rc, get_sequences

logger = logging.getLogger(__name__)

# Sentinel used internally to mark levels with no taxonomic assignment.
_UNSPECIFIED = "_DADA2_UNSPECIFIED"


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _parse_ref_fasta(path: str) -> List[Tuple[str, str]]:
    """Read a (possibly gzipped) reference FASTA file.

    Returns a list of ``(header, sequence)`` tuples.  This is a thin wrapper
    around :func:`papa2.utils._parse_fasta` kept here for clarity and to
    provide a taxonomy-module-specific entry point.
    """
    return _parse_fasta(path)


def _is_unite(headers: List[str]) -> bool:
    """Detect whether reference headers follow the UNITE format.

    UNITE headers contain ``FU|refs`` or ``FU|reps`` as a pipe-delimited
    field, e.g.::

        >SeqID|FU|refs|...|k__Fungi;p__Basidiomycota;...
    """
    # Check first few headers for the UNITE signature.
    for hdr in headers[:20]:
        if "FU|refs" in hdr or "FU|reps" in hdr:
            return True
    return False


def _parse_taxonomy_unite(header: str) -> str:
    """Extract semicolon-delimited taxonomy from a UNITE-style header.

    UNITE headers are pipe-delimited; the taxonomy string lives in the 5th
    field (0-indexed field 4).  Additionally:
      - Rank prefixes of the form ``[kpcofg]__unidentified`` are replaced
        with ``_DADA2_UNSPECIFIED``.
      - The species-level epithet (``s__``) is handled specially: if
        ``s__unidentified`` it becomes ``_DADA2_UNSPECIFIED``, otherwise the
        ``s__`` prefix is stripped but the species epithet is kept.
    """
    fields = header.split("|")
    if len(fields) < 5:
        # Fallback: treat whole header as taxonomy.
        tax_str = header.split(None, 1)[-1] if " " in header or "\t" in header else header
    else:
        tax_str = fields[4].strip()

    # Replace [kpcofg]__unidentified with _DADA2_UNSPECIFIED
    tax_str = re.sub(r"[kpcofgs]__unidentified", _UNSPECIFIED, tax_str)

    # Strip rank prefixes like k__, p__, etc.
    tax_str = re.sub(r"[kpcofgs]__", "", tax_str)

    return tax_str


def _parse_taxonomy_standard(header: str) -> str:
    """Extract semicolon-delimited taxonomy from a standard dada2 header.

    Standard format::

        >SeqID Kingdom;Phylum;Class;Order;Family;Genus;Species

    Everything after the first whitespace is treated as the taxonomy string.
    """
    parts = header.split(None, 1)
    if len(parts) < 2:
        return ""
    return parts[1].strip()


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def assign_taxonomy(
    seqs,
    ref_fasta: str,
    min_boot: int = 50,
    try_rc: bool = False,
    output_bootstraps: bool = False,
    tax_levels: Sequence[str] = (
        "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species",
    ),
    verbose: bool = False,
) -> Union[pd.DataFrame, Dict[str, pd.DataFrame]]:
    """Classify sequences against a reference taxonomy using the RDP
    Naive Bayesian Classifier method.

    This is a Python port of dada2's ``assignTaxonomy()`` R function.

    Args:
        seqs: Query sequences -- any type accepted by
            :func:`papa2.utils.get_sequences` (list of strings, FASTA path,
            pandas DataFrame, dict, etc.).
        ref_fasta: Path to a reference FASTA file (may be gzipped).  Headers
            must contain semicolon-delimited taxonomy.  UNITE-formatted
            databases are detected and parsed automatically.
        min_boot: Minimum bootstrap confidence (0--100) for retaining a
            taxonomic assignment at each rank.  Ranks below this threshold
            are set to ``None``.
        try_rc: If ``True``, sequences that fail to classify on the forward
            strand are re-tried on the reverse complement.
        output_bootstraps: If ``True``, return a dict with keys ``"tax"``
            (the taxonomy DataFrame) and ``"boot"`` (a DataFrame of bootstrap
            values).
        tax_levels: Column names for the taxonomic ranks in the output
            DataFrame.  Must match the number of semicolon-delimited levels
            in the reference database headers (shorter reference taxonomies
            are padded with ``None``).
        verbose: Print progress messages.

    Returns:
        A :class:`pandas.DataFrame` with query sequences as the index and
        ``tax_levels`` as columns.  If *output_bootstraps* is ``True``,
        returns a dict ``{"tax": <DataFrame>, "boot": <DataFrame>}``.
    """
    from ._cdada import run_taxonomy

    # ------------------------------------------------------------------
    # 1. Get query sequences
    # ------------------------------------------------------------------
    query_seqs = get_sequences(seqs)
    n_seqs = len(query_seqs)
    n_levels = len(tax_levels)

    if verbose:
        logger.info("assign_taxonomy: classifying %d sequences against %s",
                     n_seqs, ref_fasta)

    # ------------------------------------------------------------------
    # 2. Read and parse the reference database
    # ------------------------------------------------------------------
    ref_records = _parse_ref_fasta(ref_fasta)
    headers = [hdr for hdr, _ in ref_records]
    ref_sequences = [seq for _, seq in ref_records]

    unite = _is_unite(headers)
    if verbose:
        logger.info("assign_taxonomy: UNITE format detected: %s", unite)

    # Extract raw taxonomy strings
    if unite:
        raw_taxa = [_parse_taxonomy_unite(h) for h in headers]
    else:
        raw_taxa = [_parse_taxonomy_standard(h) for h in headers]

    # ------------------------------------------------------------------
    # 3. Split into levels, pad to n_levels
    # ------------------------------------------------------------------
    # Each tax string is "K;P;C;O;F;G;S" -- split on ';', stripping empty
    # trailing entries.
    split_taxa: List[List[str]] = []
    for ts in raw_taxa:
        parts = [p.strip() for p in ts.split(";") if p.strip()]
        # Pad with sentinel if shorter than n_levels
        while len(parts) < n_levels:
            parts.append(_UNSPECIFIED)
        # Truncate if longer
        parts = parts[:n_levels]
        split_taxa.append(parts)

    # ------------------------------------------------------------------
    # 4. Build integer maps for the C function
    # ------------------------------------------------------------------
    # A "genus" in dada2 terminology is a unique full taxonomy string
    # (all levels concatenated).  Build a mapping from each reference
    # sequence to a genus index, and a matrix mapping each genus to
    # integer-encoded level values.

    # First, collect unique level values per rank.
    level_to_int: List[Dict[str, int]] = [dict() for _ in range(n_levels)]
    for taxa in split_taxa:
        for lvl_idx, val in enumerate(taxa):
            if val not in level_to_int[lvl_idx]:
                level_to_int[lvl_idx][val] = len(level_to_int[lvl_idx])

    # Build the "genus" key as a tuple of all level values.
    genus_map: Dict[Tuple[str, ...], int] = {}
    ref_to_genus = np.empty(len(ref_sequences), dtype=np.int32)

    for ref_idx, taxa in enumerate(split_taxa):
        key = tuple(taxa)
        if key not in genus_map:
            genus_map[key] = len(genus_map)
        ref_to_genus[ref_idx] = genus_map[key]

    n_genus = len(genus_map)

    # genusmat: (n_genus, n_levels), integer-encoded level values.
    genusmat = np.zeros((n_genus, n_levels), dtype=np.int32)
    for key, gidx in genus_map.items():
        for lvl_idx, val in enumerate(key):
            genusmat[gidx, lvl_idx] = level_to_int[lvl_idx][val]

    # Reverse lookup: int -> string for each level.
    int_to_level: List[Dict[int, str]] = []
    for lvl_idx in range(n_levels):
        int_to_level.append({v: k for k, v in level_to_int[lvl_idx].items()})

    # Reverse lookup: genus index -> taxonomy tuple.
    idx_to_genus: Dict[int, Tuple[str, ...]] = {v: k for k, v in genus_map.items()}

    if verbose:
        logger.info(
            "assign_taxonomy: %d references, %d unique genera, %d levels",
            len(ref_sequences), n_genus, n_levels,
        )

    # ------------------------------------------------------------------
    # 5. Call the C classifier
    # ------------------------------------------------------------------
    result = run_taxonomy(
        query_seqs, ref_sequences,
        ref_to_genus, genusmat,
        n_genus, n_levels,
        verbose=verbose,
    )
    rval = result["rval"]    # (n_seqs,) int32, 1-indexed best genus (0=NA)
    rboot = result["rboot"]  # (n_seqs, n_levels) int32, bootstrap counts

    # ------------------------------------------------------------------
    # 5b. Optionally try reverse complement for unassigned sequences
    # ------------------------------------------------------------------
    if try_rc:
        unassigned = np.where(rval == 0)[0]
        if len(unassigned) > 0:
            rc_seqs = [_rc(query_seqs[i]) for i in unassigned]
            if verbose:
                logger.info(
                    "assign_taxonomy: trying reverse complement for %d "
                    "unassigned sequences", len(unassigned),
                )
            rc_result = run_taxonomy(
                rc_seqs, ref_sequences,
                ref_to_genus, genusmat,
                n_genus, n_levels,
                verbose=verbose,
            )
            rc_rval = rc_result["rval"]
            rc_rboot = rc_result["rboot"]
            # Fill in previously unassigned positions
            for j, orig_idx in enumerate(unassigned):
                if rc_rval[j] != 0:
                    rval[orig_idx] = rc_rval[j]
                    rboot[orig_idx, :] = rc_rboot[j, :]

    # ------------------------------------------------------------------
    # 6. Parse results into taxonomy strings, masking low-bootstrap ranks
    # ------------------------------------------------------------------
    tax_out = np.full((n_seqs, n_levels), None, dtype=object)
    boot_out = np.zeros((n_seqs, n_levels), dtype=np.int64)

    for i in range(n_seqs):
        gidx = rval[i]
        if gidx == 0:
            # No assignment
            continue
        # rval is 1-indexed; convert to 0-indexed.
        gidx_0 = gidx - 1
        taxa_tuple = idx_to_genus.get(gidx_0)
        if taxa_tuple is None:
            continue

        for lvl in range(n_levels):
            boot_out[i, lvl] = rboot[i, lvl]
            val = taxa_tuple[lvl]
            # Mask ranks below min_boot or with the unspecified sentinel.
            if val == _UNSPECIFIED:
                tax_out[i, lvl] = None
            elif rboot[i, lvl] < min_boot:
                # Once a rank is masked, all deeper ranks should also be
                # masked (dada2 behaviour).
                tax_out[i, lvl:] = None
                break
            else:
                tax_out[i, lvl] = val

    # ------------------------------------------------------------------
    # 7. Build output DataFrame(s)
    # ------------------------------------------------------------------
    tax_levels_list = list(tax_levels)

    tax_df = pd.DataFrame(tax_out, columns=tax_levels_list, index=query_seqs)
    tax_df.index.name = "sequence"

    if verbose:
        assigned_counts = [
            int((tax_df[col].notna()).sum()) for col in tax_levels_list
        ]
        for col, cnt in zip(tax_levels_list, assigned_counts):
            logger.info("assign_taxonomy: %s assigned: %d / %d",
                         col, cnt, n_seqs)

    if output_bootstraps:
        boot_df = pd.DataFrame(
            boot_out, columns=tax_levels_list, index=query_seqs,
        )
        boot_df.index.name = "sequence"
        return {"tax": tax_df, "boot": boot_df}

    return tax_df
