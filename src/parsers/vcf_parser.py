
import gzip

def scan_vcf(vcf_dir):
    total = 0
    for vcf in vcf_dir.glob("*.vcf*"):
        opener = gzip.open if vcf.suffix.endswith("gz") else open
        with opener(vcf, "rt") as f:
            for line in f:
                if not line.startswith("#"):
                    total += 1
    return {"total": total}
