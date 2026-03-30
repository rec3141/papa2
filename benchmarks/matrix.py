from __future__ import annotations

from dataclasses import dataclass, asdict


@dataclass(frozen=True)
class BenchmarkCase:
    case_id: str
    sample_count: int
    source_name: str
    source_dir: str

    def to_dict(self) -> dict:
        return asdict(self)


DATA_SOURCES = {
    "raw_big": "/tmp/dada2_gpu_test/raw_big",
    "raw_40m": "/tmp/dada2_gpu_test/raw_40m",
}


def _cases_for_counts(source_name: str, counts: list[int]) -> list[BenchmarkCase]:
    source_dir = DATA_SOURCES[source_name]
    return [
        BenchmarkCase(
            case_id=f"{source_name}_s{sample_count}",
            sample_count=sample_count,
            source_name=source_name,
            source_dir=source_dir,
        )
        for sample_count in counts
    ]


def build_matrix(name: str) -> list[BenchmarkCase]:
    if name == "smoke":
        return _cases_for_counts("raw_big", [1, 10])
    if name == "small":
        return (
            _cases_for_counts("raw_big", [1, 10, 100])
            + _cases_for_counts("raw_40m", [1000])
        )
    if name == "full":
        return (
            _cases_for_counts("raw_big", [1, 10, 100])
            + _cases_for_counts("raw_40m", [1000, 3000])
        )
    raise ValueError(f"Unknown matrix: {name}")
