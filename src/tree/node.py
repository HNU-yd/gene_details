from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

@dataclass
class RegionNode:
    node_type: str
    node_id: str
    chrom: Optional[str] = None
    start: Optional[int] = None
    end: Optional[int] = None

    variant_count: int = 0
    missing_N_count: int = 0
    missing_N_runs: int = 0
    length: Optional[int] = None

    meta: Dict[str, Any] = field(default_factory=dict)
    children: List["RegionNode"] = field(default_factory=list)

    def add(self, child: "RegionNode") -> None:
        self.children.append(child)

    def compute_lengths(self) -> None:
        for c in self.children:
            c.compute_lengths()
        if self.start is not None and self.end is not None:
            self.length = self.end - self.start + 1
        elif self.children:
            self.length = sum((c.length or 0) for c in self.children)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "type": self.node_type,
            "id": self.node_id,
            "chrom": self.chrom,
            "start": self.start,
            "end": self.end,
            "length": self.length,
            "variant_count": int(self.variant_count),
            "missing_N_count": int(self.missing_N_count),
            "missing_N_runs": int(self.missing_N_runs),
            "meta": self.meta,
            "children": [c.to_dict() for c in self.children],
        }
