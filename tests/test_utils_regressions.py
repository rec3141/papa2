from pathlib import Path

import numpy as np

import papa2


def test_assign_species_finds_non_binned_substring(tmp_path):
    ref = tmp_path / "ref.fa"
    ref.write_text(">id1 Genus species\nAACCGGTTAACCGGTT\n")

    out = papa2.assign_species(["ACCGGTTA"], str(ref))

    assert out[0, 0] == "Genus"
    assert out[0, 1] == "species"


def test_assign_species_handles_short_query(tmp_path):
    ref = tmp_path / "ref.fa"
    ref.write_text(">id1 Genus species\nAACCGGTTAACCGGTT\n")

    out = papa2.assign_species(["CCGGT"], str(ref))

    assert out[0, 0] == "Genus"
    assert out[0, 1] == "species"


def test_merge_sequence_tables_try_rc_uses_majority_orientation():
    t1 = {
        "table": np.array([[100]], dtype=np.int64),
        "seqs": ["GTTTT"],
        "sample_names": ["s1"],
    }
    t2 = {
        "table": np.array([[1]], dtype=np.int64),
        "seqs": ["AAAAC"],
        "sample_names": ["s2"],
    }

    out1 = papa2.merge_sequence_tables(t1, t2, try_rc=True)
    out2 = papa2.merge_sequence_tables(t2, t1, try_rc=True)

    assert out1["seqs"] == ["GTTTT"]
    assert out2["seqs"] == ["GTTTT"]
    assert np.array_equal(out1["table"], np.array([[100], [1]], dtype=np.int64))
    assert np.array_equal(out2["table"], np.array([[1], [100]], dtype=np.int64))


def test_merge_sequence_tables_preserves_sample_names():
    t1 = {
        "table": np.array([[1, 2]], dtype=np.int64),
        "seqs": ["AAA", "CCC"],
        "sample_names": ["sample_a"],
    }
    t2 = {
        "table": np.array([[3]], dtype=np.int64),
        "seqs": ["GGG"],
        "sample_names": ["sample_b"],
    }

    out = papa2.merge_sequence_tables(t1, t2)

    assert out["sample_names"] == ["sample_a", "sample_b"]
    assert out["table"].shape == (2, 3)
