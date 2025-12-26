from pathlib import Path
from typing import Iterator, List, Tuple
import gzip
from utils.naming import normalize_chrom

def _open_text(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt", encoding="utf-8", errors="replace")

def iter_vcf_variants(paths: List[Path]) -> Iterator[Tuple[str,int]]:
    for p in paths:
        with _open_text(p) as f:
            for line in f:
                if not line or line.startswith("#"):
                    continue
                cols = line.split("\\t")
                if len(cols) < 2:
                    continue
                chrom = normalize_chrom(cols[0])
                try:
                    pos = int(cols[1])
                except ValueError:
                    continue
                yield chrom, pos
