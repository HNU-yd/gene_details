#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Genome stats tree builder (per-chrom outputs + final species summary).

Input layout (typical):
  data/<species>/
    *.fa / *.fa.gz                  (one or many)
    *.gff / *.gff3 / *.gff.gz ...   (one or per-chrom many)
    *.vcf / *.vcf.gz                (one or per-chrom many)

Outputs:
  ./stats/<species>/
    chrom_<chrom>.tree.json.gz
    chrom_<chrom>.nodes.tsv.gz
    species.total.summary.json.gz
    species.total.summary.nodes.tsv.gz

Progress:
  --progress
  --progress-every N  (GFF lines / VCF variants)

Multithreading:
  --threads N  (parallel per-gene subtree build within each chromosome)

Notes:
  - gene-level uses UNION across transcripts.
  - leaf segments: cds_segment, utr_segment, exon_other_segment, intron_segment, non_gene_segment
  - per-node: coverage intervals (merged), covered_length, missing_N_ratio over coverage, variant counts via overlap.
  - gene_union_group and non_gene_group variant_count are UNIQUE counts (by definition of union / complement).
  - chromosome.variant_count is UNIQUE per VCF record in that chrom.
  - For final species summary: only keep chromosome-level stats (drop children) to avoid huge memory.
"""

from __future__ import annotations

import argparse
import bisect
import dataclasses
import gzip
import io
import json
import os
import re
import time
import glob
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Any, Dict, Iterable, List, Optional, Tuple


# =========================
# Logging / IO
# =========================

def log(msg: str) -> None:
    print(msg, flush=True)


def open_text_auto(path: str) -> io.TextIOBase:
    if path.endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def normalize_chrom(name: str) -> str:
    s = name.strip()
    # correct inline flag usage
    s = re.sub(r"(?i)^chr", "", s)
    return s


def is_dir(p: str) -> bool:
    return os.path.isdir(p)


def expand_inputs(path_or_glob: Optional[str], base_dir: str, exts: Tuple[str, ...]) -> List[str]:
    """
    If None: scan base_dir for matching extensions.
    If directory: scan that directory.
    If file: return [file]
    If glob: expand.
    """
    def scan_dir(d: str) -> List[str]:
        out = []
        for fn in os.listdir(d):
            low = fn.lower()
            if any(low.endswith(ext) for ext in exts):
                out.append(os.path.join(d, fn))
        return sorted(out)

    if path_or_glob is None:
        return scan_dir(base_dir)

    p = path_or_glob
    if any(ch in p for ch in ["*", "?", "[", "]"]):
        return sorted(glob.glob(p))
    if is_dir(p):
        return scan_dir(p)
    return [p]


# =========================
# Interval ops (1-based inclusive)
# =========================

Interval = Tuple[int, int]


def interval_len(iv: Interval) -> int:
    return iv[1] - iv[0] + 1


def merge_intervals(intervals: List[Interval]) -> List[Interval]:
    if not intervals:
        return []
    intervals = sorted(intervals, key=lambda x: (x[0], x[1]))
    out: List[Interval] = []
    cs, ce = intervals[0]
    for s, e in intervals[1:]:
        if s <= ce + 1:
            ce = max(ce, e)
        else:
            out.append((cs, ce))
            cs, ce = s, e
    out.append((cs, ce))
    return out


def subtract_intervals(span: Interval, cut: List[Interval]) -> List[Interval]:
    s0, e0 = span
    if s0 > e0:
        return []
    cut = merge_intervals([(max(s0, s), min(e0, e)) for s, e in cut if not (e < s0 or s > e0)])
    if not cut:
        return [(s0, e0)]
    out: List[Interval] = []
    cur = s0
    for s, e in cut:
        if cur < s:
            out.append((cur, s - 1))
        cur = max(cur, e + 1)
        if cur > e0:
            break
    if cur <= e0:
        out.append((cur, e0))
    return out


def intersect_with_interval(intervals: List[Interval], target: Interval) -> List[Interval]:
    if not intervals:
        return []
    ts, te = target
    ends = [e for _, e in intervals]
    i = bisect.bisect_left(ends, ts)
    out: List[Interval] = []
    while i < len(intervals):
        s, e = intervals[i]
        if s > te:
            break
        os_ = max(s, ts)
        oe_ = min(e, te)
        if os_ <= oe_:
            out.append((os_, oe_))
        i += 1
    return out


def overlaps_any(merged: List[Interval], qs: int, qe: int) -> bool:
    if not merged:
        return False
    ends = [e for _, e in merged]
    i = bisect.bisect_left(ends, qs)
    if i >= len(merged):
        return False
    s, e = merged[i]
    return not (e < qs or s > qe)


# =========================
# Stats accumulators
# =========================

@dataclasses.dataclass
class Welford:
    n: int = 0
    mean: float = 0.0
    m2: float = 0.0
    min_v: float = float("inf")
    max_v: float = float("-inf")

    def add(self, x: float) -> None:
        self.n += 1
        self.min_v = min(self.min_v, x)
        self.max_v = max(self.max_v, x)
        delta = x - self.mean
        self.mean += delta / self.n
        delta2 = x - self.mean
        self.m2 += delta * delta2

    def finalize(self) -> Dict[str, float]:
        if self.n == 0:
            return {"min": 0.0, "max": 0.0, "mean": 0.0, "var": 0.0, "n": 0}
        # population variance
        return {"min": float(self.min_v), "max": float(self.max_v), "mean": float(self.mean), "var": float(self.m2 / self.n), "n": int(self.n)}


# =========================
# Node model
# =========================

LEAF_TYPES = {"cds_segment", "utr_segment", "exon_other_segment", "intron_segment", "non_gene_segment"}


@dataclasses.dataclass
class Node:
    type: str
    id: str
    chrom: Optional[str] = None
    start: Optional[int] = None
    end: Optional[int] = None

    coverage: List[Interval] = dataclasses.field(default_factory=list)

    covered_length: int = 0
    length_span: int = 0

    variant_count: int = 0
    variant_count_fixed: bool = False

    missing_N_count: int = 0
    missing_N_ratio: float = 0.0
    missing_N_runs: int = 0

    child_quant: Dict[str, Any] = dataclasses.field(default_factory=dict)
    desc_leaf_quant: Dict[str, Any] = dataclasses.field(default_factory=dict)

    children: List["Node"] = dataclasses.field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "type": self.type,
            "id": self.id,
            "chrom": self.chrom,
            "start": self.start,
            "end": self.end,
            "length_span": self.length_span,
            "coverage": self.coverage,
            "covered_length": self.covered_length,
            "variant_count": self.variant_count,
            "missing_N_count": self.missing_N_count,
            "missing_N_ratio": self.missing_N_ratio,
            "missing_N_runs": self.missing_N_runs,
            "child_quant": self.child_quant,
            "desc_leaf_quant": self.desc_leaf_quant,
            "children": [c.to_dict() for c in self.children],
        }


# =========================
# FASTA N-runs
# =========================

@dataclasses.dataclass
class ChromSeqIndex:
    chrom: str
    length: int
    n_runs: List[Interval]


def build_fasta_n_runs_one(fasta_path: str, progress: bool) -> Dict[str, ChromSeqIndex]:
    """
    Parse one FASTA file (may contain multiple records), return chrom -> ChromSeqIndex.
    """
    idx: Dict[str, ChromSeqIndex] = {}

    def commit_record(chrom: Optional[str], length: int, runs: List[Interval]) -> None:
        if chrom is None:
            return
        chrom_n = normalize_chrom(chrom)
        merged = merge_intervals(runs)
        idx[chrom_n] = ChromSeqIndex(chrom=chrom_n, length=length, n_runs=merged)
        if progress:
            log(f"[fasta] done chrom={chrom_n} length={length} N_runs={len(merged)} file={os.path.basename(fasta_path)}")

    chrom: Optional[str] = None
    pos = 0
    runs: List[Interval] = []
    in_run = False
    run_start = 0

    with open_text_auto(fasta_path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if in_run:
                    runs.append((run_start, pos))
                    in_run = False
                commit_record(chrom, pos, runs)
                header = line[1:].strip()
                chrom = header.split()[0]
                pos = 0
                runs = []
                in_run = False
                run_start = 0
                if progress:
                    log(f"[fasta] start raw_id={chrom} file={os.path.basename(fasta_path)}")
                continue

            seq = line.strip()
            i = 0
            L = len(seq)
            while i < L:
                j = seq.find("N", i)
                if j == -1:
                    if in_run:
                        runs.append((run_start, pos))
                        in_run = False
                    pos += (L - i)
                    break

                if j > i:
                    if in_run:
                        runs.append((run_start, pos))
                        in_run = False
                    pos += (j - i)

                k = j
                while k < L and seq[k] == "N":
                    k += 1

                if not in_run:
                    in_run = True
                    run_start = pos + 1
                pos += (k - j)

                if k < L:
                    runs.append((run_start, pos))
                    in_run = False

                i = k

    if in_run:
        runs.append((run_start, pos))
    commit_record(chrom, pos, runs)
    return idx


def merge_fasta_indices(idxs: List[Dict[str, ChromSeqIndex]], progress: bool) -> Dict[str, ChromSeqIndex]:
    """
    Merge multiple fasta indices.
    If same chrom appears multiple times:
      - keep the one with larger length
      - if length equal: merge N-runs
    """
    out: Dict[str, ChromSeqIndex] = {}
    for idx in idxs:
        for chrom, ci in idx.items():
            if chrom not in out:
                out[chrom] = ci
                continue
            prev = out[chrom]
            if ci.length > prev.length:
                if progress:
                    log(f"[fasta] warn chrom={chrom} length {prev.length}->{ci.length} (keep longer)")
                out[chrom] = ci
            elif ci.length == prev.length:
                merged_runs = merge_intervals(prev.n_runs + ci.n_runs)
                out[chrom] = ChromSeqIndex(chrom=chrom, length=prev.length, n_runs=merged_runs)
    return out


def n_overlap_stats(n_runs: List[Interval], qs: int, qe: int) -> Tuple[int, int]:
    if not n_runs or qs > qe:
        return 0, 0
    ends = [e for _, e in n_runs]
    i = bisect.bisect_left(ends, qs)
    n_count = 0
    n_hit = 0
    while i < len(n_runs):
        s, e = n_runs[i]
        if s > qe:
            break
        os_ = max(s, qs)
        oe_ = min(e, qe)
        if os_ <= oe_:
            n_count += (oe_ - os_ + 1)
            n_hit += 1
        i += 1
    return n_count, n_hit


def n_overlap_over_coverage(n_runs: List[Interval], coverage: List[Interval]) -> Tuple[int, int]:
    total_n = 0
    total_runs = 0
    for s, e in coverage:
        n_cnt, n_hit = n_overlap_stats(n_runs, s, e)
        total_n += n_cnt
        total_runs += n_hit
    return total_n, total_runs


# =========================
# Infer chrom from GFF/VCF file
# =========================

def infer_chrom_from_gff(path: str, max_lines: int = 5000) -> Optional[str]:
    """
    Try:
      - ##sequence-region <seqid> ...
      - first non-comment feature line: seqid in col0
    """
    try:
        with open_text_auto(path) as f:
            for i, line in enumerate(f):
                if i >= max_lines:
                    break
                line = line.strip()
                if not line:
                    continue
                if line.startswith("##sequence-region"):
                    parts = line.split()
                    if len(parts) >= 2:
                        return normalize_chrom(parts[1])
                    continue
                if line.startswith("#"):
                    continue
                cols = line.split("\t")
                if len(cols) >= 1:
                    return normalize_chrom(cols[0])
    except Exception:
        return None
    return None


def infer_chrom_from_vcf(path: str, max_lines: int = 20000) -> Optional[str]:
    try:
        with open_text_auto(path) as f:
            for i, line in enumerate(f):
                if i >= max_lines:
                    break
                if not line:
                    continue
                if line.startswith("#"):
                    continue
                cols = line.rstrip("\n").split("\t")
                if cols and cols[0]:
                    return normalize_chrom(cols[0])
    except Exception:
        return None
    return None


# =========================
# GFF parsing (two-pass for robustness)
# =========================

@dataclasses.dataclass
class TranscriptModel:
    tid: str
    exons: List[Interval] = dataclasses.field(default_factory=list)
    cds: List[Interval] = dataclasses.field(default_factory=list)
    utr: List[Interval] = dataclasses.field(default_factory=list)


@dataclasses.dataclass
class GeneModel:
    gid: str
    chrom: str
    start: int
    end: int
    strand: str
    transcripts: Dict[str, TranscriptModel] = dataclasses.field(default_factory=dict)


def parse_gff_attributes(attr: str) -> Dict[str, str]:
    d: Dict[str, str] = {}
    for part in attr.split(";"):
        part = part.strip()
        if not part:
            continue
        if "=" in part:
            k, v = part.split("=", 1)
            d[k] = v
        else:
            d[part] = ""
    return d


def strip_prefix_id(raw: str) -> str:
    return raw.split(":", 1)[-1] if ":" in raw else raw


def parse_gff_two_pass_onefile(
    gff_path: str,
    progress: bool,
    progress_every: int,
) -> Dict[str, Dict[str, GeneModel]]:
    """
    Return genes_by_chrom for this file.
    Two-pass:
      pass1: collect genes, transcript->gene mapping
      pass2: fill exons/cds/utr
    """
    genes_by_chrom: Dict[str, Dict[str, GeneModel]] = {}
    transcript_to_gene: Dict[str, Tuple[str, str]] = {}  # tid -> (chrom,gid)

    t0 = time.time()

    # pass 1
    line_n = 0
    gene_n = 0
    tx_n = 0
    with open_text_auto(gff_path) as f:
        for line in f:
            line_n += 1
            if progress and (line_n % progress_every == 0):
                log(f"[gff-1] file={os.path.basename(gff_path)} lines={line_n:,} genes={gene_n:,} tx={tx_n:,}")

            line = line.strip()
            if not line or line.startswith("#") or line == "###":
                continue
            cols = line.split("\t")
            if len(cols) < 9:
                continue
            seqid, source, ftype, s, e, score, strand, phase, attr = cols
            chrom = normalize_chrom(seqid)
            start = int(s)
            end = int(e)
            attrs = parse_gff_attributes(attr)

            if ftype == "gene":
                gid_raw = attrs.get("ID", "") or attrs.get("gene_id", "")
                if not gid_raw:
                    continue
                gid = strip_prefix_id(gid_raw)
                genes_by_chrom.setdefault(chrom, {})
                genes_by_chrom[chrom][gid] = GeneModel(gid=gid, chrom=chrom, start=start, end=end, strand=strand, transcripts={})
                gene_n += 1

            elif ftype in ("mRNA", "transcript"):
                tid_raw = attrs.get("ID", "")
                parent_raw = attrs.get("Parent", "")
                if not tid_raw or not parent_raw:
                    continue
                tid = strip_prefix_id(tid_raw)
                parent_gid = strip_prefix_id(parent_raw).split(",")[0]
                transcript_to_gene[tid] = (chrom, parent_gid)
                gm = genes_by_chrom.get(chrom, {}).get(parent_gid)
                if gm is not None:
                    gm.transcripts.setdefault(tid, TranscriptModel(tid=tid))
                tx_n += 1

    # pass 2
    feat_n = 0
    line2 = 0
    with open_text_auto(gff_path) as f:
        for line in f:
            line2 += 1
            if progress and (line2 % progress_every == 0):
                log(f"[gff-2] file={os.path.basename(gff_path)} lines={line2:,} feats={feat_n:,}")

            line = line.strip()
            if not line or line.startswith("#") or line == "###":
                continue
            cols = line.split("\t")
            if len(cols) < 9:
                continue
            seqid, source, ftype, s, e, score, strand, phase, attr = cols
            if ftype not in ("exon", "CDS", "five_prime_UTR", "three_prime_UTR", "UTR"):
                continue

            chrom = normalize_chrom(seqid)
            start = int(s)
            end = int(e)
            attrs = parse_gff_attributes(attr)
            parent_raw = attrs.get("Parent", "")
            if not parent_raw:
                continue
            parents = [strip_prefix_id(x) for x in parent_raw.split(",")]

            for tid in parents:
                cg = transcript_to_gene.get(tid)
                if cg is None:
                    continue
                chrom2, gid = cg
                if chrom2 != chrom:
                    continue
                gm = genes_by_chrom.get(chrom, {}).get(gid)
                if gm is None:
                    continue
                tm = gm.transcripts.setdefault(tid, TranscriptModel(tid=tid))
                if ftype == "exon":
                    tm.exons.append((start, end))
                elif ftype == "CDS":
                    tm.cds.append((start, end))
                else:
                    tm.utr.append((start, end))
                feat_n += 1

    # merge transcript intervals
    for chrom, gmap in genes_by_chrom.items():
        for gm in gmap.values():
            for tm in gm.transcripts.values():
                tm.exons = merge_intervals(tm.exons)
                tm.cds = merge_intervals(tm.cds)
                tm.utr = merge_intervals(tm.utr)

    if progress:
        log(f"[gff] done file={os.path.basename(gff_path)} genes={gene_n:,} tx={tx_n:,} feats={feat_n:,} elapsed={time.time()-t0:.1f}s")

    return genes_by_chrom


# =========================
# Build gene subtree (union across transcripts)
# =========================

@dataclasses.dataclass
class LeafSeg:
    start: int
    end: int
    seg_type: str
    node: Node


def build_exon_subsegments(exon_iv: Interval, cds_cov: List[Interval], utr_cov: List[Interval]) -> Tuple[List[Interval], List[Interval], List[Interval]]:
    cds_in = merge_intervals(intersect_with_interval(cds_cov, exon_iv))
    utr_in = merge_intervals(intersect_with_interval(utr_cov, exon_iv))

    estart, eend = exon_iv
    cut_points = {estart, eend + 1}
    for s, e in cds_in:
        cut_points.add(s); cut_points.add(e + 1)
    for s, e in utr_in:
        cut_points.add(s); cut_points.add(e + 1)

    pts = sorted(p for p in cut_points if estart <= p <= eend + 1)
    segs: List[Interval] = []
    for i in range(len(pts) - 1):
        s = pts[i]
        e = pts[i + 1] - 1
        if s <= e:
            segs.append((s, e))

    def covered_by(merged_cov: List[Interval], seg: Interval) -> bool:
        if not merged_cov:
            return False
        s, e = seg
        starts = [a for a, _ in merged_cov]
        j = bisect.bisect_right(starts, s) - 1
        if j < 0:
            return False
        a, b = merged_cov[j]
        return a <= s and b >= e

    cds_segs: List[Interval] = []
    utr_segs: List[Interval] = []
    other_segs: List[Interval] = []

    for seg in segs:
        if covered_by(cds_in, seg):
            cds_segs.append(seg)
        elif covered_by(utr_in, seg):
            utr_segs.append(seg)
        else:
            other_segs.append(seg)

    return merge_intervals(cds_segs), merge_intervals(utr_segs), merge_intervals(other_segs)


def compute_lengths_and_missing(node: Node, n_runs: List[Interval]) -> None:
    for c in node.children:
        compute_lengths_and_missing(c, n_runs)

    if not node.coverage and node.children:
        cov: List[Interval] = []
        for c in node.children:
            cov.extend(c.coverage)
        node.coverage = merge_intervals(cov)

    if node.start is None or node.end is None:
        if node.coverage:
            node.start = min(s for s, _ in node.coverage)
            node.end = max(e for _, e in node.coverage)

    node.length_span = (node.end - node.start + 1) if (node.start is not None and node.end is not None) else 0
    node.covered_length = sum(interval_len(iv) for iv in node.coverage) if node.coverage else 0

    if node.coverage:
        n_cnt, n_hit = n_overlap_over_coverage(n_runs, node.coverage)
        node.missing_N_count = int(n_cnt)
        node.missing_N_runs = int(n_hit)
        node.missing_N_ratio = (n_cnt / node.covered_length) if node.covered_length > 0 else 0.0


def build_gene_node(gm: GeneModel, n_runs: List[Interval]) -> Tuple[Node, List[LeafSeg]]:
    chrom = gm.chrom
    gene_id = gm.gid

    gene_node = Node(type="gene", id=f"gene:{gene_id}", chrom=chrom, start=gm.start, end=gm.end, coverage=[(gm.start, gm.end)])

    all_exons: List[Interval] = []
    all_cds: List[Interval] = []
    all_utr: List[Interval] = []
    for tm in gm.transcripts.values():
        all_exons.extend(tm.exons)
        all_cds.extend(tm.cds)
        all_utr.extend(tm.utr)

    exon_union = merge_intervals(all_exons)
    cds_union = merge_intervals(all_cds)
    utr_union = merge_intervals(all_utr)

    exon_group = Node(type="exon_group", id=f"{gene_node.id}:exon_group", chrom=chrom, coverage=exon_union)
    non_exon_group = Node(type="non_exon_group", id=f"{gene_node.id}:non_exon_group", chrom=chrom)

    leafs: List[LeafSeg] = []

    for k, exon_iv in enumerate(exon_union, start=1):
        es_node = Node(type="exon_segment", id=f"{gene_node.id}:exon_seg:{k}", chrom=chrom, start=exon_iv[0], end=exon_iv[1], coverage=[exon_iv])

        cds_segs, utr_segs, other_segs = build_exon_subsegments(exon_iv, cds_union, utr_union)

        cds_group = Node(type="cds_group", id=f"{es_node.id}:cds_group", chrom=chrom, coverage=cds_segs)
        utr_group = Node(type="utr_group", id=f"{es_node.id}:utr_group", chrom=chrom, coverage=utr_segs)
        other_group = Node(type="exon_other_group", id=f"{es_node.id}:exon_other_group", chrom=chrom, coverage=other_segs)

        for i, iv in enumerate(cds_segs, start=1):
            nd = Node(type="cds_segment", id=f"{cds_group.id}:seg:{i}", chrom=chrom, start=iv[0], end=iv[1], coverage=[iv])
            cds_group.children.append(nd)
            leafs.append(LeafSeg(iv[0], iv[1], "cds_segment", nd))

        for i, iv in enumerate(utr_segs, start=1):
            nd = Node(type="utr_segment", id=f"{utr_group.id}:seg:{i}", chrom=chrom, start=iv[0], end=iv[1], coverage=[iv])
            utr_group.children.append(nd)
            leafs.append(LeafSeg(iv[0], iv[1], "utr_segment", nd))

        for i, iv in enumerate(other_segs, start=1):
            nd = Node(type="exon_other_segment", id=f"{other_group.id}:seg:{i}", chrom=chrom, start=iv[0], end=iv[1], coverage=[iv])
            other_group.children.append(nd)
            leafs.append(LeafSeg(iv[0], iv[1], "exon_other_segment", nd))

        es_node.children.extend([cds_group, utr_group, other_group])
        exon_group.children.append(es_node)

    # intron = gene_span - exon_union
    introns = subtract_intervals((gm.start, gm.end), exon_union)
    non_exon_group.coverage = introns
    for i, iv in enumerate(introns, start=1):
        nd = Node(type="intron_segment", id=f"{non_exon_group.id}:seg:{i}", chrom=chrom, start=iv[0], end=iv[1], coverage=[iv])
        non_exon_group.children.append(nd)
        leafs.append(LeafSeg(iv[0], iv[1], "intron_segment", nd))

    gene_node.children.extend([exon_group, non_exon_group])
    compute_lengths_and_missing(gene_node, n_runs)
    return gene_node, leafs


# =========================
# Quantification (child_quant + desc_leaf_quant)
# =========================

def compute_child_quant(node: Node) -> None:
    for c in node.children:
        compute_child_quant(c)

    by_type: Dict[str, Dict[str, Any]] = {}
    acc_len: Dict[str, Welford] = {}
    acc_miss: Dict[str, Welford] = {}
    var_sum: Dict[str, int] = {}
    count: Dict[str, int] = {}

    for c in node.children:
        t = c.type
        count[t] = count.get(t, 0) + 1
        var_sum[t] = var_sum.get(t, 0) + int(c.variant_count)
        acc_len.setdefault(t, Welford()).add(float(c.covered_length))
        acc_miss.setdefault(t, Welford()).add(float(c.missing_N_ratio))

    for t in count.keys():
        ls = acc_len[t].finalize()
        ms = acc_miss[t].finalize()
        by_type[t] = {
            "count": int(count[t]),
            "len_min": ls["min"],
            "len_max": ls["max"],
            "len_mean": ls["mean"],
            "len_var": ls["var"],
            "missing_ratio_min": ms["min"],
            "missing_ratio_max": ms["max"],
            "missing_ratio_mean": ms["mean"],
            "missing_ratio_var": ms["var"],
            "variant_sum": int(var_sum[t]),
        }

    node.child_quant = {"by_type": by_type}


def compute_desc_leaf_quant(node: Node) -> Dict[str, Any]:
    leaf_acc: Dict[str, Dict[str, int]] = {}
    if node.type in LEAF_TYPES:
        leaf_acc[node.type] = {"count": 1, "variant_sum": int(node.variant_count)}
    else:
        for c in node.children:
            child_leaf = compute_desc_leaf_quant(c)
            for lt, vv in child_leaf.items():
                cur = leaf_acc.setdefault(lt, {"count": 0, "variant_sum": 0})
                cur["count"] += int(vv["count"])
                cur["variant_sum"] += int(vv["variant_sum"])
    node.desc_leaf_quant = {"by_leaf_type": leaf_acc}
    return leaf_acc


# =========================
# Variant counting (per-chrom VCF)
# =========================

def parse_vcf_variant_interval(pos_str: str, ref: str, info: str) -> Interval:
    pos = int(pos_str)
    end = pos + len(ref) - 1
    m = re.search(r"(?:^|;)END=(\d+)(?:;|$)", info)
    if m:
        try:
            end = int(m.group(1))
        except Exception:
            pass
    if end < pos:
        end = pos
    return pos, end


def build_gene_bin_index(genes: List[Tuple[int, int, str, Node, List[LeafSeg]]], bin_size: int) -> Dict[int, List[int]]:
    bins: Dict[int, List[int]] = {}
    for idx, (s, e, gid, node, leafs) in enumerate(genes):
        b0 = (s - 1) // bin_size
        b1 = (e - 1) // bin_size
        for b in range(b0, b1 + 1):
            bins.setdefault(b, []).append(idx)
    return bins


def add_variant_to_gene_leafs(var_s: int, var_e: int, leafs: List[LeafSeg]) -> None:
    if not leafs:
        return
    starts = [x.start for x in leafs]
    i = bisect.bisect_right(starts, var_s) - 1
    if i < 0:
        i = 0
    while i > 0 and leafs[i - 1].end >= var_s:
        i -= 1
    while i < len(leafs):
        seg = leafs[i]
        if seg.start > var_e:
            break
        if not (seg.end < var_s or seg.start > var_e):
            seg.node.variant_count += 1
        i += 1


def count_variants_for_chrom(
    vcf_paths: List[str],
    chrom: str,
    chrom_node: Node,
    gene_union: List[Interval],
    gene_union_group: Node,
    non_gene_group: Node,
    non_gene_leafs: List[LeafSeg],
    genes_list: List[Tuple[int, int, str, Node, List[LeafSeg]]],
    gene_bins: Dict[int, List[int]],
    bin_size: int,
    progress: bool,
    progress_every: int,
) -> int:
    """
    Return number of VCF records processed for this chrom.
    """
    t0 = time.time()
    n = 0
    for vp in vcf_paths:
        if progress:
            log(f"[vcf] start chrom={chrom} file={os.path.basename(vp)}")
        with open_text_auto(vp) as f:
            for line in f:
                if not line or line.startswith("#"):
                    continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 8:
                    continue
                chrom_raw, pos, vid, ref, alt, qual, flt, info = cols[:8]
                c = normalize_chrom(chrom_raw)
                if c != chrom:
                    # for safety: if file contains multiple chroms, ignore others
                    continue

                var_s, var_e = parse_vcf_variant_interval(pos, ref, info)

                chrom_node.variant_count += 1
                if overlaps_any(gene_union, var_s, var_e):
                    gene_union_group.variant_count += 1
                else:
                    non_gene_group.variant_count += 1
                    # fast: assign by var_s containment to one non_gene segment
                    starts = [x.start for x in non_gene_leafs]
                    j = bisect.bisect_right(starts, var_s) - 1
                    if 0 <= j < len(non_gene_leafs) and non_gene_leafs[j].start <= var_s <= non_gene_leafs[j].end:
                        non_gene_leafs[j].node.variant_count += 1

                # gene overlaps via bins
                b0 = (var_s - 1) // bin_size
                b1 = (var_e - 1) // bin_size
                cand = []
                for b in range(b0, b1 + 1):
                    cand.extend(gene_bins.get(b, []))
                if cand:
                    seen = set()
                    for gi in cand:
                        if gi in seen:
                            continue
                        seen.add(gi)
                        gs, ge, gid, gnode, leafs = genes_list[gi]
                        if ge < var_s or gs > var_e:
                            continue
                        add_variant_to_gene_leafs(var_s, var_e, leafs)

                n += 1
                if progress and (n % progress_every == 0):
                    log(f"[vcf] chrom={chrom} variants={n:,} elapsed={time.time()-t0:.1f}s")

    if progress:
        log(f"[vcf] done chrom={chrom} variants={n:,} elapsed={time.time()-t0:.1f}s")
    return n


def aggregate_variant_counts(node: Node) -> int:
    if not node.children:
        return node.variant_count
    s = 0
    for c in node.children:
        s += aggregate_variant_counts(c)
    if not node.variant_count_fixed:
        node.variant_count = int(s)
    return node.variant_count


# =========================
# Flat table export
# =========================

def iter_nodes_preorder(node: Node, depth: int = 0, parent_id: str = "") -> Iterable[Tuple[Node, int, str]]:
    yield node, depth, parent_id
    for c in node.children:
        yield from iter_nodes_preorder(c, depth + 1, node.id)


def write_nodes_tsv_gz(path: str, root: Node) -> None:
    header = [
        "depth", "parent_id", "id", "type", "chrom",
        "start", "end", "length_span", "covered_length",
        "variant_count",
        "missing_N_count", "missing_N_ratio", "missing_N_runs",
    ]
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8") as w:
        w.write("\t".join(header) + "\n")
        for nd, depth, pid in iter_nodes_preorder(root):
            row = [
                str(depth),
                pid,
                nd.id,
                nd.type,
                nd.chrom or "",
                str(nd.start or ""),
                str(nd.end or ""),
                str(nd.length_span),
                str(nd.covered_length),
                str(nd.variant_count),
                str(nd.missing_N_count),
                f"{nd.missing_N_ratio:.8f}",
                str(nd.missing_N_runs),
            ]
            w.write("\t".join(row) + "\n")


def write_json_gz(path: str, obj: Dict[str, Any]) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8") as w:
        json.dump(obj, w, ensure_ascii=False)


# =========================
# Per-chrom build
# =========================

def build_chromosome_tree(
    species: str,
    chrom: str,
    seq: ChromSeqIndex,
    genes_by_chrom: Dict[str, Dict[str, GeneModel]],
    vcf_paths: List[str],
    bin_size: int,
    threads: int,
    progress: bool,
    progress_every: int,
) -> Node:
    """
    Build one chromosome tree, count variants, compute quants, return chrom_node.
    """
    chrom_len = seq.length
    chrom_node = Node(
        type="chromosome",
        id=f"chrom:{chrom}",
        chrom=chrom,
        start=1,
        end=chrom_len,
        coverage=[(1, chrom_len)],
        variant_count=0,
        variant_count_fixed=True,
    )

    # build genes
    gene_models = genes_by_chrom.get(chrom, {})
    gm_items = sorted(gene_models.items(), key=lambda kv: (kv[1].start, kv[1].end, kv[0]))
    if progress:
        log(f"[build] chrom={chrom} genes={len(gm_items):,} threads={threads}")

    t0 = time.time()
    gene_nodes: List[Node] = []
    genes_list: List[Tuple[int, int, str, Node, List[LeafSeg]]] = []
    gene_spans: List[Interval] = []

    if threads > 1 and len(gm_items) > 0:
        results: Dict[str, Tuple[Node, List[LeafSeg]]] = {}
        with ThreadPoolExecutor(max_workers=threads) as ex:
            futs = {ex.submit(build_gene_node, gm, seq.n_runs): gid for gid, gm in gm_items}
            done = 0
            for fut in as_completed(futs):
                gid = futs[fut]
                gnode, leafs = fut.result()
                results[gid] = (gnode, leafs)
                done += 1
                if progress and (done % 5000 == 0):
                    log(f"[build] chrom={chrom} built_genes={done:,}/{len(gm_items):,} elapsed={time.time()-t0:.1f}s")

        for gid, gm in gm_items:
            gnode, leafs = results[gid]
            gene_nodes.append(gnode)
            genes_list.append((gm.start, gm.end, gid, gnode, leafs))
            gene_spans.append((gm.start, gm.end))
    else:
        for i, (gid, gm) in enumerate(gm_items, start=1):
            gnode, leafs = build_gene_node(gm, seq.n_runs)
            gene_nodes.append(gnode)
            genes_list.append((gm.start, gm.end, gid, gnode, leafs))
            gene_spans.append((gm.start, gm.end))
            if progress and (i % 5000 == 0):
                log(f"[build] chrom={chrom} built_genes={i:,}/{len(gm_items):,} elapsed={time.time()-t0:.1f}s")

    # gene union and non-gene
    gene_union = merge_intervals(gene_spans)
    gene_union_group = Node(
        type="gene_union_group",
        id=f"{chrom_node.id}:gene_union_group",
        chrom=chrom,
        coverage=gene_union,
        variant_count=0,
        variant_count_fixed=True,
    )
    gene_union_group.children.extend(gene_nodes)

    non_gene_intervals = subtract_intervals((1, chrom_len), gene_union)
    non_gene_group = Node(
        type="non_gene_group",
        id=f"{chrom_node.id}:non_gene_group",
        chrom=chrom,
        coverage=non_gene_intervals,
        variant_count=0,
        variant_count_fixed=True,
    )

    non_gene_leafs: List[LeafSeg] = []
    for i, iv in enumerate(non_gene_intervals, start=1):
        nd = Node(type="non_gene_segment", id=f"{non_gene_group.id}:seg:{i}", chrom=chrom, start=iv[0], end=iv[1], coverage=[iv])
        non_gene_group.children.append(nd)
        non_gene_leafs.append(LeafSeg(iv[0], iv[1], "non_gene_segment", nd))

    chrom_node.children.extend([gene_union_group, non_gene_group])

    # missing stats
    compute_lengths_and_missing(chrom_node, seq.n_runs)

    # bins
    gene_bins = build_gene_bin_index(genes_list, bin_size=bin_size)

    # count variants
    count_variants_for_chrom(
        vcf_paths=vcf_paths,
        chrom=chrom,
        chrom_node=chrom_node,
        gene_union=gene_union,
        gene_union_group=gene_union_group,
        non_gene_group=non_gene_group,
        non_gene_leafs=non_gene_leafs,
        genes_list=genes_list,
        gene_bins=gene_bins,
        bin_size=bin_size,
        progress=progress,
        progress_every=progress_every,
    )

    # aggregate + quant
    aggregate_variant_counts(chrom_node)
    compute_child_quant(chrom_node)
    compute_desc_leaf_quant(chrom_node)

    if progress:
        log(f"[build] chrom={chrom} done elapsed={time.time()-t0:.1f}s union_blocks={len(gene_union)} non_gene_blocks={len(non_gene_intervals)}")
    return chrom_node


def make_shallow_chrom_summary(chrom_node: Node) -> Node:
    """
    Create a lightweight copy of chromosome node:
      - keep chromosome-level stats and quant fields
      - drop children to keep memory small for species summary
    """
    nd = Node(
        type=chrom_node.type,
        id=chrom_node.id,
        chrom=chrom_node.chrom,
        start=chrom_node.start,
        end=chrom_node.end,
        coverage=chrom_node.coverage,
        covered_length=chrom_node.covered_length,
        length_span=chrom_node.length_span,
        variant_count=chrom_node.variant_count,
        variant_count_fixed=True,
        missing_N_count=chrom_node.missing_N_count,
        missing_N_ratio=chrom_node.missing_N_ratio,
        missing_N_runs=chrom_node.missing_N_runs,
        child_quant=chrom_node.child_quant,
        desc_leaf_quant=chrom_node.desc_leaf_quant,
        children=[],
    )
    return nd


# =========================
# Main
# =========================

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--species-dir", required=True, help="Path to data/<species>/ directory")
    ap.add_argument("--species-name", default=None, help="Species name (default: basename of species-dir)")
    ap.add_argument("--fasta", default=None, help="FASTA file/dir/glob (default: scan species-dir)")
    ap.add_argument("--gff", default=None, help="GFF file/dir/glob (default: scan species-dir)")
    ap.add_argument("--vcf", default=None, help="VCF file/dir/glob (default: scan species-dir)")
    ap.add_argument("--bin-size", type=int, default=100000, help="Gene bin size for overlap candidate lookup")
    ap.add_argument("--threads", type=int, default=1, help="Worker threads for building per-gene subtrees")
    ap.add_argument("--progress", action="store_true", help="Enable progress logs")
    ap.add_argument("--progress-every", type=int, default=1_000_000, help="Progress interval for GFF lines / VCF variants")

    ap.add_argument("--stats-root", default="./stats", help="Output root directory (default: ./stats)")
    args = ap.parse_args()

    species_dir = args.species_dir
    species = args.species_name or os.path.basename(os.path.abspath(species_dir))

    # collect inputs
    fasta_files = expand_inputs(args.fasta, species_dir, (".fa", ".fasta", ".fa.gz", ".fasta.gz"))
    gff_files = expand_inputs(args.gff, species_dir, (".gff3", ".gff", ".gff3.gz", ".gff.gz"))
    vcf_files = expand_inputs(args.vcf, species_dir, (".vcf", ".vcf.gz"))

    if not fasta_files:
        raise SystemExit(f"No FASTA found (species-dir={species_dir})")
    if not gff_files:
        log("[warn] No GFF found; will treat whole chromosome as non-gene.")
    if not vcf_files:
        log("[warn] No VCF found; variant counts will be 0.")

    out_dir = os.path.join(args.stats_root, species)
    os.makedirs(out_dir, exist_ok=True)

    log(f"[info] species      : {species}")
    log(f"[info] species_dir  : {species_dir}")
    log(f"[info] fasta_files  : {len(fasta_files)}")
    log(f"[info] gff_files    : {len(gff_files)}")
    log(f"[info] vcf_files    : {len(vcf_files)}")
    log(f"[info] out_dir      : {out_dir}")
    log(f"[info] threads      : {args.threads}")
    log(f"[info] progress     : {args.progress}")
    log(f"[info] progress_every: {args.progress_every}")
    log(f"[info] bin_size     : {args.bin_size}")

    t_all = time.time()

    # build FASTA index (possibly many files)
    fasta_idxs = []
    for fp in fasta_files:
        if args.progress:
            log(f"[fasta] file={os.path.basename(fp)}")
        fasta_idxs.append(build_fasta_n_runs_one(fp, progress=args.progress))
    fasta_idx = merge_fasta_indices(fasta_idxs, progress=args.progress)
    if args.progress:
        log(f"[fasta] total chromosomes indexed: {len(fasta_idx)}")

    # map gff/vcf files -> chrom by peeking
    gff_by_chrom: Dict[str, List[str]] = {}
    for gp in gff_files:
        c = infer_chrom_from_gff(gp)
        if c is None:
            log(f"[warn] cannot infer chrom for gff: {gp}")
            continue
        gff_by_chrom.setdefault(c, []).append(gp)

    vcf_by_chrom: Dict[str, List[str]] = {}
    for vp in vcf_files:
        c = infer_chrom_from_vcf(vp)
        if c is None:
            log(f"[warn] cannot infer chrom for vcf: {vp}")
            continue
        vcf_by_chrom.setdefault(c, []).append(vp)

    # process chromosomes present in FASTA
    chroms = sorted(fasta_idx.keys(), key=lambda x: (len(x), x))
    log(f"[plan] chromosomes to process (from FASTA): {len(chroms)}")

    # final species summary will keep shallow chrom nodes only
    species_summary = Node(type="species", id=f"species:{species}")

    for chrom in chroms:
        seq = fasta_idx[chrom]

        # parse gff only for this chrom (from mapped files); if missing => empty genes
        genes_by_chrom: Dict[str, Dict[str, GeneModel]] = {}
        gffs = gff_by_chrom.get(chrom, [])
        if gffs:
            # merge multiple files if provided
            merged: Dict[str, Dict[str, GeneModel]] = {}
            for gp in gffs:
                one = parse_gff_two_pass_onefile(gp, progress=args.progress, progress_every=args.progress_every)
                for c, gmap in one.items():
                    merged.setdefault(c, {})
                    merged[c].update(gmap)
            genes_by_chrom = merged
        else:
            genes_by_chrom = {chrom: {}}
            if args.progress:
                log(f"[gff] chrom={chrom} no gff files, genes=0")

        vcfs = vcf_by_chrom.get(chrom, [])
        if not vcfs and args.progress:
            log(f"[vcf] chrom={chrom} no vcf files, variants=0")

        log(f"[chrom] start chrom={chrom} gff_files={len(gffs)} vcf_files={len(vcfs)} len={seq.length}")

        chrom_node = build_chromosome_tree(
            species=species,
            chrom=chrom,
            seq=seq,
            genes_by_chrom=genes_by_chrom,
            vcf_paths=vcfs,
            bin_size=args.bin_size,
            threads=max(1, int(args.threads)),
            progress=args.progress,
            progress_every=max(1, int(args.progress_every)),
        )

        # write per-chrom outputs
        out_tree = os.path.join(out_dir, f"chrom_{chrom}.tree.json.gz")
        out_tsv = os.path.join(out_dir, f"chrom_{chrom}.nodes.tsv.gz")
        write_json_gz(out_tree, chrom_node.to_dict())
        write_nodes_tsv_gz(out_tsv, chrom_node)

        # add shallow chrom summary to species summary
        species_summary.children.append(make_shallow_chrom_summary(chrom_node))

        log(f"[chrom] done chrom={chrom} -> {os.path.basename(out_tree)} / {os.path.basename(out_tsv)}")

        # release memory for this chrom tree explicitly
        del chrom_node

    # compute species summary stats from children (shallow)
    # coverage/span for species node: not meaningful; we keep only aggregate variant_count and child_quant/leaf_quant.
    # Aggregate variant_count as sum of chromosome unique variant_count.
    species_summary.variant_count_fixed = True
    species_summary.variant_count = sum(c.variant_count for c in species_summary.children)

    compute_child_quant(species_summary)
    compute_desc_leaf_quant(species_summary)

    out_species_tree = os.path.join(out_dir, "species.total.summary.json.gz")
    out_species_tsv = os.path.join(out_dir, "species.total.summary.nodes.tsv.gz")
    write_json_gz(out_species_tree, species_summary.to_dict())
    write_nodes_tsv_gz(out_species_tsv, species_summary)

    log(f"[done] total_elapsed={time.time()-t_all:.1f}s")
    log(f"[done] species summary -> {out_species_tree} / {out_species_tsv}")


if __name__ == "__main__":
    main()
'''
python stats_ensembl_tree.py `
  --species-dir ./data/species/日本晴 `
  --threads 1 `
  --progress `
  --progress-every 1000000`
  --stats-root ./stats
'''