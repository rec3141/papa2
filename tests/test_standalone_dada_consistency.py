import importlib
from pathlib import Path

import numpy as np

py_dada = importlib.import_module("papa2.dada")
from papa2.io import derep_fastq


def test_standalone_abundances_match_assigned_reads():
    repo_root = Path(__file__).resolve().parents[1]
    files = [
        str(repo_root / "tests" / "data" / "sam1F.fastq.gz"),
        str(repo_root / "tests" / "data" / "sam2F.fastq.gz"),
    ]

    err = py_dada.learn_errors(files, verbose=False)
    derep = derep_fastq(files[0], verbose=False)
    res = py_dada.dada(derep, err=err, verbose=False)

    assigned = np.asarray(res["map"]) >= 0
    assigned_reads = int(derep["abundances"][assigned].sum())
    denoised_reads = int(sum(res["denoised"].values()))
    cluster_reads = int(np.asarray(res["cluster_abunds"]).sum())

    assert denoised_reads == assigned_reads
    assert cluster_reads == assigned_reads
