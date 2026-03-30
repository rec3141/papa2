import importlib

import numpy as np

dada_mod = importlib.import_module("papa2.dada")


def test_file_inputs_stream_derep_and_denoise(monkeypatch):
    calls = []

    def fake_derep_fastq(filepath, verbose=False):
        calls.append(("derep", filepath))
        return {
            "seqs": ["ACGT"],
            "abundances": np.array([1], dtype=np.int32),
            "quals": np.array([[30.0, 30.0, 30.0, 30.0]], dtype=np.float64),
        }

    def fake_run_one_sample(args):
        drp, erri, opts, max_clust, verbose = args
        calls.append(("run", tuple(drp["seqs"])))
        return {
            "cluster_seqs": ["ACGT"],
            "cluster_abunds": np.array([1], dtype=np.int32),
            "trans": np.zeros((16, erri.shape[1]), dtype=np.int32),
            "map": np.array([0], dtype=np.int32),
            "pval": np.array([1.0], dtype=np.float64),
            "denoised": {"ACGT": 1},
        }

    monkeypatch.setattr(dada_mod, "derep_fastq", fake_derep_fastq)
    monkeypatch.setattr(dada_mod, "_run_one_sample", fake_run_one_sample)
    monkeypatch.setenv("DADA2_WORKERS", "1")

    err = np.ones((16, 41), dtype=np.float64)
    files = ["sample_a.fastq.gz", "sample_b.fastq.gz"]
    res = dada_mod.dada(files, err=err, verbose=False)

    assert len(res) == 2
    assert calls == [
        ("derep", "sample_a.fastq.gz"),
        ("run", ("ACGT",)),
        ("derep", "sample_b.fastq.gz"),
        ("run", ("ACGT",)),
    ]
