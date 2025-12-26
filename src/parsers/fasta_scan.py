from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import gzip
from utils.naming import normalize_chrom

@dataclass
class FastaChromInfo:
    name: str
    length: int
    missing_N_count: int
    missing_N_runs: int
    n_runs: List[Tuple[int,int]]  # 1-based inclusive

def _open_text(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt", encoding="utf-8", errors="replace")

def scan_fasta_files(paths: List[Path]) -> Dict[str, FastaChromInfo]:
    out: Dict[str, FastaChromInfo] = {}
    for p in paths:
        with _open_text(p) as f:
            cur = None
            pos = 0
            n_count = 0
            n_runs = 0
            runs: List[Tuple[int,int]] = []
            run_start: Optional[int] = None

            def flush():
                nonlocal cur, pos, n_count, n_runs, runs, run_start
                if cur is None:
                    return
                if run_start is not None:
                    runs.append((run_start, pos))
                    n_runs += 1
                    run_start = None
                out[cur] = FastaChromInfo(cur, pos, n_count, n_runs, runs)

            for line in f:
                if line.startswith(">"):
                    flush()
                    header = line[1:].strip().split()[0]
                    cur = normalize_chrom(header)
                    pos = 0; n_count = 0; n_runs = 0; runs = []; run_start = None
                    continue
                if cur is None:
                    continue
                seq = line.strip().upper()
                if not seq:
                    continue
                for ch in seq:
                    pos += 1
                    if ch == "N":
                        n_count += 1
                        if run_start is None:
                            run_start = pos
                    else:
                        if run_start is not None:
                            runs.append((run_start, pos-1))
                            n_runs += 1
                            run_start = None
            flush()
    return out
