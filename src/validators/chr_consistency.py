
def check_chr_consistency(fasta, gff, vcf, species):
    if not fasta.get("chr_set"):
        raise ValueError(f"[{species}] No chromosomes detected")
