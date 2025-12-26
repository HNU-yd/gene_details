import re

def normalize_chrom(name: str) -> str:
    n = name.strip()
    if n.lower().startswith("chr"):
        return "chr" + n[3:]
    if re.fullmatch(r"\d+", n):
        return "chr" + n
    u = n.upper()
    if u in ("X","Y","M","MT"):
        return "chr" + ("M" if u in ("M","MT") else u)
    return n
