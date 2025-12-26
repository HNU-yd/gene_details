from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import gzip
from collections import defaultdict
from utils.naming import normalize_chrom
from utils.intervals import merge_intervals

def _open_text(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt", encoding="utf-8", errors="replace")

def _parse_attrs(attr: str) -> Dict[str,str]:
    out = {}
    for part in attr.split(";"):
        if not part:
            continue
        if "=" in part:
            k,v = part.split("=",1)
            out[k] = v
    return out

@dataclass
class GeneModel:
    gene_id: str
    start: int
    end: int
    biotype: Optional[str]
    exon: List[Tuple[int,int]]
    cds: List[Tuple[int,int]]
    utr: List[Tuple[int,int]]

def parse_gff_files(paths: List[Path]) -> Dict[str, Dict[str, GeneModel]]:
    genes: Dict[str, Tuple[str,int,int,Optional[str]]] = {}
    tx2gene: Dict[str, str] = {}
    exon_by_gene: Dict[str, List[Tuple[int,int]]] = defaultdict(list)
    cds_by_gene: Dict[str, List[Tuple[int,int]]] = defaultdict(list)
    utr_by_gene: Dict[str, List[Tuple[int,int]]] = defaultdict(list)

    for p in paths:
        with _open_text(p) as f:
            for line in f:
                if not line or line.startswith("#"):
                    continue
                cols = line.rstrip("\\n").split("\\t")
                if len(cols) < 9:
                    continue
                seqid, source, ftype, s, e, score, strand, phase, attr = cols[:9]
                chrom = normalize_chrom(seqid)
                try:
                    start = int(s); end = int(e)
                except ValueError:
                    continue
                attrs = _parse_attrs(attr)

                if ftype == "gene":
                    gid = attrs.get("ID") or attrs.get("gene_id") or attrs.get("Name")
                    if not gid:
                        continue
                    biotype = attrs.get("gene_biotype") or attrs.get("biotype") or attrs.get("gene_type")
                    genes[gid] = (chrom, start, end, biotype)
                elif ftype in ("transcript","mRNA"):
                    tid = attrs.get("ID")
                    parent = attrs.get("Parent")
                    if tid and parent:
                        tx2gene[tid] = parent.split(",")[0]
                elif ftype == "exon":
                    parent = attrs.get("Parent")
                    if not parent:
                        continue
                    tid = parent.split(",")[0]
                    gid = tx2gene.get(tid)
                    if gid:
                        exon_by_gene[gid].append((start,end))
                elif ftype == "CDS":
                    parent = attrs.get("Parent")
                    if not parent:
                        continue
                    tid = parent.split(",")[0]
                    gid = tx2gene.get(tid)
                    if gid:
                        cds_by_gene[gid].append((start,end))
                elif ftype in ("five_prime_UTR","three_prime_UTR","UTR"):
                    parent = attrs.get("Parent")
                    if not parent:
                        continue
                    tid = parent.split(",")[0]
                    gid = tx2gene.get(tid)
                    if gid:
                        utr_by_gene[gid].append((start,end))

    by_chrom: Dict[str, Dict[str, GeneModel]] = {}
    for gid,(chrom,start,end,biotype) in genes.items():
        gm = GeneModel(
            gene_id=gid,
            start=start,
            end=end,
            biotype=biotype,
            exon=merge_intervals(exon_by_gene.get(gid, [])),
            cds=merge_intervals(cds_by_gene.get(gid, [])),
            utr=merge_intervals(utr_by_gene.get(gid, [])),
        )
        by_chrom.setdefault(chrom, {})[gid] = gm
    return by_chrom
