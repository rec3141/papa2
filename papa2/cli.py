"""Command-line interface for papa2.

Usage:
    papa2 filter-trim <fwd> <filt_fwd> <rev> <filt_rev> [options]
    papa2 assign-taxonomy <seqs> <ref_fasta> <output_prefix> [options]
    papa2 --version
"""

import argparse
import os
import pickle
import sys

import pandas as pd

from . import __version__


def _cmd_filter_trim(args):
    """Run paired-end quality filtering."""
    from .filter import fastq_paired_filter

    reads_in, reads_out = fastq_paired_filter(
        args.fwd, args.filt_fwd,
        args.rev, args.filt_rev,
        max_ee=(args.max_ee, args.max_ee),
        trunc_q=(args.trunc_q, args.trunc_q) if args.trunc_q is not None else (None, None),
        max_n=(args.max_n, args.max_n),
        trunc_len=(args.trunc_len_fwd, args.trunc_len_rev),
        min_len=(args.min_len, args.min_len),
        rm_phix=not args.no_rm_phix,
        compress=not args.no_compress,
        verbose=args.verbose,
    )

    # Write stats TSV if requested
    if args.stats:
        pct = round(reads_out / max(reads_in, 1) * 100, 1)
        with open(args.stats, "w") as f:
            f.write("sample\treads_in\treads_out\tpct_retained\n")
            f.write(f"{args.sample_id or os.path.basename(args.fwd)}\t"
                    f"{reads_in}\t{reads_out}\t{pct}\n")

    # Handle zero-output: create empty files so downstream can detect
    if reads_out == 0:
        for p in (args.filt_fwd, args.filt_rev):
            if os.path.exists(p):
                os.remove(p)
            open(p, "w").close()

    pct = round(reads_out / max(reads_in, 1) * 100, 1)
    print(f"{reads_out}/{reads_in} ({pct}%) passed filter")


def _cmd_auto_trim(args):
    """Analyze quality profiles and recommend truncation lengths."""
    import gzip
    import glob

    import numpy as np

    fwd_files = sorted(glob.glob(args.input_dir + "/*_1.fastq.gz") +
                       glob.glob(args.input_dir + "/*_R1*.fastq.gz"))
    rev_files = sorted(glob.glob(args.input_dir + "/*_2.fastq.gz") +
                       glob.glob(args.input_dir + "/*_R2*.fastq.gz"))

    if not fwd_files or not rev_files:
        print("ERROR: no paired FASTQ files found in " + args.input_dir,
              file=sys.stderr)
        sys.exit(1)

    # Sample from a subset of files
    n_files = min(args.n_files, len(fwd_files))
    step = max(1, len(fwd_files) // n_files)
    fwd_sample = fwd_files[::step][:n_files]
    rev_sample = rev_files[::step][:n_files]

    def _read_quals(files, n_reads):
        """Read quality scores from FASTQ files, return per-position arrays."""
        all_quals = []
        count = 0
        for fpath in files:
            opener = gzip.open if fpath.endswith(".gz") else open
            with opener(fpath, "rt") as fh:
                while count < n_reads:
                    header = fh.readline()
                    if not header:
                        break
                    fh.readline()  # seq
                    fh.readline()  # +
                    qual_str = fh.readline().rstrip()
                    if not qual_str:
                        break
                    all_quals.append([ord(c) - 33 for c in qual_str])
                    count += 1
            if count >= n_reads:
                break
        return all_quals

    if args.verbose:
        print(f"Sampling {args.n_reads} reads from {n_files} file pairs...",
              file=sys.stderr)

    fwd_quals = _read_quals(fwd_sample, args.n_reads)
    rev_quals = _read_quals(rev_sample, args.n_reads)

    def _find_trunc_pos(quals_list, min_q, window):
        """Find position where rolling median quality drops below threshold."""
        if not quals_list:
            return 0
        max_len = max(len(q) for q in quals_list)
        # Build quality matrix (NaN-padded)
        mat = np.full((len(quals_list), max_len), np.nan)
        for i, q in enumerate(quals_list):
            mat[i, :len(q)] = q
        medians = np.nanmedian(mat, axis=0)
        # Rolling mean of medians over window
        if len(medians) < window:
            return len(medians)
        rolling = np.convolve(medians, np.ones(window) / window, mode="valid")
        for i, val in enumerate(rolling):
            if val < min_q:
                return i + window // 2
        return len(medians)

    trunc_fwd = _find_trunc_pos(fwd_quals, args.min_quality, args.window)
    trunc_rev = _find_trunc_pos(rev_quals, args.min_quality, args.window)

    # Compute summary stats for reporting
    def _median_at(quals_list, pos):
        vals = [q[pos] for q in quals_list if len(q) > pos]
        return float(np.median(vals)) if vals else 0.0

    fwd_len = max(len(q) for q in fwd_quals) if fwd_quals else 0
    rev_len = max(len(q) for q in rev_quals) if rev_quals else 0

    if args.verbose:
        print(f"Forward: max_len={fwd_len}, median_Q@end={_median_at(fwd_quals, fwd_len-1):.0f}",
              file=sys.stderr)
        print(f"Reverse: max_len={rev_len}, median_Q@end={_median_at(rev_quals, rev_len-1):.0f}",
              file=sys.stderr)

    # Output as key=value for easy parsing
    print(f"trunc_len_fwd={trunc_fwd}")
    print(f"trunc_len_rev={trunc_rev}")
    print(f"fwd_read_len={fwd_len}")
    print(f"rev_read_len={rev_len}")
    print(f"n_reads_sampled={len(fwd_quals)}")

    # Write to file if requested
    if args.output:
        with open(args.output, "w") as f:
            f.write(f"trunc_len_fwd\t{trunc_fwd}\n")
            f.write(f"trunc_len_rev\t{trunc_rev}\n")
            f.write(f"fwd_read_len\t{fwd_len}\n")
            f.write(f"rev_read_len\t{rev_len}\n")
            f.write(f"min_quality\t{args.min_quality}\n")
            f.write(f"n_reads_sampled\t{len(fwd_quals)}\n")


def _cmd_assign_taxonomy(args):
    """Run taxonomy assignment."""
    from .taxonomy import assign_taxonomy

    # Load query sequences
    data = pickle.load(open(args.seqs, "rb"))
    if isinstance(data, pd.DataFrame):
        seqs = sorted(data["sequence"].unique())
    elif isinstance(data, list):
        seqs = data
    else:
        seqs = list(data)

    if args.verbose:
        print(f"{len(seqs)} query sequences", file=sys.stderr)

    os.environ["OMP_NUM_THREADS"] = str(args.threads)

    kw = dict(
        seqs=seqs,
        ref_fasta=args.ref_fasta,
        min_boot=args.min_boot,
        try_rc=args.try_rc,
        output_bootstraps=True,
        verbose=args.verbose,
    )
    if args.tax_levels:
        kw["tax_levels"] = args.tax_levels

    result = assign_taxonomy(**kw)
    tax_df = result["tax"]
    boot_df = result["boot"]

    prefix = args.output_prefix

    # Save pickle outputs
    pickle.dump(result, open(f"{prefix}_taxonomy.pkl", "wb"))
    pickle.dump(boot_df, open(f"{prefix}_bootstrap.pkl", "wb"))

    # Save TSV
    tax_out = tax_df.copy()
    tax_out.insert(0, "sequence", tax_out.index)
    tax_out.to_csv(f"{prefix}_taxonomy.tsv", sep="\t", index=False)

    # Summary
    for col in tax_df.columns:
        n = int(tax_df[col].notna().sum())
        pct = round(n / max(len(tax_df), 1) * 100, 1)
        print(f"{col}: {n}/{len(tax_df)} ({pct}%)")


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="papa2",
        description="Amplicon denoising and analysis toolkit",
    )
    parser.add_argument("--version", action="version", version=f"papa2 {__version__}")
    sub = parser.add_subparsers(dest="command")

    # -- filter-trim --
    ft = sub.add_parser("filter-trim", help="Paired-end quality filtering")
    ft.add_argument("fwd", help="Forward input FASTQ")
    ft.add_argument("filt_fwd", help="Filtered forward output FASTQ")
    ft.add_argument("rev", help="Reverse input FASTQ")
    ft.add_argument("filt_rev", help="Filtered reverse output FASTQ")
    ft.add_argument("--max-ee", type=float, default=float("inf"), help="Max expected errors [default: inf]")
    ft.add_argument("--trunc-q", type=int, default=None, help="Truncate at first base with Q <= this")
    ft.add_argument("--max-n", type=int, default=0, help="Max ambiguous bases [default: 0]")
    ft.add_argument("--trunc-len-fwd", type=int, default=0, help="Truncate forward reads to length [default: 0 = off]")
    ft.add_argument("--trunc-len-rev", type=int, default=0, help="Truncate reverse reads to length [default: 0 = off]")
    ft.add_argument("--min-len", type=int, default=20, help="Min read length [default: 20]")
    ft.add_argument("--no-rm-phix", action="store_true", help="Skip PhiX removal")
    ft.add_argument("--no-compress", action="store_true", help="Don't gzip output")
    ft.add_argument("--stats", help="Write filtering stats TSV to this path")
    ft.add_argument("--sample-id", help="Sample ID for stats TSV [default: input filename]")
    ft.add_argument("-v", "--verbose", action="store_true")

    # -- auto-trim --
    aq = sub.add_parser("auto-trim", help="Analyze quality profiles and recommend truncation lengths")
    aq.add_argument("input_dir", help="Directory containing paired FASTQ files")
    aq.add_argument("--min-quality", type=float, default=25, help="Min rolling median Q to keep [default: 25]")
    aq.add_argument("--window", type=int, default=10, help="Rolling window size for quality smoothing [default: 10]")
    aq.add_argument("--n-reads", type=int, default=10000, help="Number of reads to sample [default: 10000]")
    aq.add_argument("--n-files", type=int, default=20, help="Number of files to sample from [default: 20]")
    aq.add_argument("-o", "--output", help="Write results to TSV file")
    aq.add_argument("-v", "--verbose", action="store_true")

    # -- assign-taxonomy --
    at = sub.add_parser("assign-taxonomy", help="Naive Bayesian taxonomy assignment")
    at.add_argument("seqs", help="Query sequences (pickle file)")
    at.add_argument("ref_fasta", help="Reference FASTA (may be gzipped)")
    at.add_argument("output_prefix", help="Output prefix for .pkl and .tsv files")
    at.add_argument("--tax-levels", nargs="+", help="Taxonomy level names")
    at.add_argument("--min-boot", type=int, default=50, help="Min bootstrap confidence [default: 50]")
    at.add_argument("--try-rc", action="store_true", help="Try reverse complement for unassigned")
    at.add_argument("--threads", type=int, default=1, help="Number of threads [default: 1]")
    at.add_argument("-v", "--verbose", action="store_true")

    args = parser.parse_args(argv)

    if args.command is None:
        parser.print_help()
        sys.exit(1)

    if args.command == "filter-trim":
        _cmd_filter_trim(args)
    elif args.command == "auto-trim":
        _cmd_auto_trim(args)
    elif args.command == "assign-taxonomy":
        _cmd_assign_taxonomy(args)
