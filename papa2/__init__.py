"""papa2 Python bindings."""

from .dada import dada, learn_errors, DADA_OPTS, set_dada_opt, get_dada_opt
from .io import derep_fastq
from .error import loess_errfun, noqual_errfun, inflate_err
from ._cdada import nwalign, eval_pair, pair_consensus, rc
from .paired import merge_pairs
from .chimera import remove_bimera_denovo
from .utils import (
    assign_species,
    add_species,
    collapse_no_mismatch,
    make_sequence_table,
    plot_quality_profile,
    plot_sankey,
    track_reads,
    uniquesto_fasta,
    write_fasta,
    is_phix,
    seq_complexity,
    get_sequences,
    get_uniques,
)

__version__ = "0.1.0"
