from __future__ import annotations

from dataclasses import dataclass, asdict


@dataclass(frozen=True)
class BenchmarkCase:
    case_id: str
    total_bases: int
    sample_count: int
    per_sample_bases: int
    learn_nbases: int

    def to_dict(self) -> dict:
        return asdict(self)


def build_matrix(name: str) -> list[BenchmarkCase]:
    if name == "smoke":
        totals = [10_000, 100_000]
        samples = [1, 10]
    elif name == "small":
        totals = [10_000, 100_000, 1_000_000]
        samples = [1, 10, 100]
    elif name == "full":
        totals = [10_000, 100_000, 1_000_000, 10_000_000, 100_000_000, 1_000_000_000]
        samples = [1, 10, 100, 1000]
    else:
        raise ValueError(f"Unknown matrix: {name}")

    cases: list[BenchmarkCase] = []
    for total in totals:
        for sample_count in samples:
            per_sample = max(1, (total + sample_count - 1) // sample_count)
            learn_nbases = min(total, 100_000_000)
            case_id = f"b{total}_s{sample_count}"
            cases.append(
                BenchmarkCase(
                    case_id=case_id,
                    total_bases=total,
                    sample_count=sample_count,
                    per_sample_bases=per_sample,
                    learn_nbases=learn_nbases,
                )
            )
    return cases
