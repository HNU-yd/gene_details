from typing import List, Tuple
import bisect

Interval = Tuple[int, int]

def merge_intervals(intervals: List[Interval]) -> List[Interval]:
    if not intervals:
        return []
    intervals = sorted(intervals, key=lambda x: (x[0], x[1]))
    out = []
    s0, e0 = intervals[0]
    for s, e in intervals[1:]:
        if s <= e0 + 1:
            e0 = max(e0, e)
        else:
            out.append((s0, e0))
            s0, e0 = s, e
    out.append((s0, e0))
    return out

def complement_intervals(span: Interval, covered: List[Interval]) -> List[Interval]:
    s0, e0 = span
    if s0 > e0:
        return []
    covered = merge_intervals([(max(s0,s), min(e0,e)) for s,e in covered if not (e < s0 or s > e0)])
    out = []
    cur = s0
    for s, e in covered:
        if cur < s:
            out.append((cur, s-1))
        cur = max(cur, e+1)
    if cur <= e0:
        out.append((cur, e0))
    return out

def intersect_intervals(a: List[Interval], b: List[Interval]) -> List[Interval]:
    a = merge_intervals(a)
    b = merge_intervals(b)
    i = j = 0
    out = []
    while i < len(a) and j < len(b):
        s1,e1 = a[i]; s2,e2 = b[j]
        s = max(s1,s2); e = min(e1,e2)
        if s <= e:
            out.append((s,e))
        if e1 < e2:
            i += 1
        else:
            j += 1
    return out

def subtract_intervals(base: List[Interval], sub: List[Interval]) -> List[Interval]:
    base = merge_intervals(base)
    sub = merge_intervals(sub)
    out = []
    j = 0
    for s, e in base:
        cur = s
        while j < len(sub) and sub[j][1] < cur:
            j += 1
        k = j
        while k < len(sub) and sub[k][0] <= e:
            ss, ee = sub[k]
            if ss > cur:
                out.append((cur, min(e, ss-1)))
            cur = max(cur, ee+1)
            if cur > e:
                break
            k += 1
        if cur <= e:
            out.append((cur, e))
    return merge_intervals(out)

def interval_total_len(iv: List[Interval]) -> int:
    return sum(e-s+1 for s,e in iv)

def count_hit(disjoint_segs, pos: int):
    # disjoint_segs: List[(s,e,node)]
    if not disjoint_segs:
        return None
    starts = [s for s,_,_ in disjoint_segs]
    i = bisect.bisect_right(starts, pos) - 1
    if i >= 0:
        s,e,node = disjoint_segs[i]
        if s <= pos <= e:
            return node
    return None
