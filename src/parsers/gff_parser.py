
import gzip

def scan_gff(gff_dir):
    gene_count = 0
    for gff in gff_dir.glob("*.gff*"):
        opener = gzip.open if gff.suffix.endswith("gz") else open
        with opener(gff, "rt") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                if line.split("\t")[2] == "gene":
                    gene_count += 1
    return {"total": gene_count}
