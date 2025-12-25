
import gzip

def scan_fasta(fasta_dir):
    chromosomes = set()
    genome_len = 0
    for fa in fasta_dir.glob("*.fa*"):
        opener = gzip.open if fa.suffix.endswith("gz") else open
        with opener(fa, "rt") as f:
            for line in f:
                if line.startswith(">"):
                    chromosomes.add(line[1:].split()[0])
                else:
                    genome_len += len(line.strip())
    return {
        "chromosomes": len(chromosomes),
        "genome_length": genome_len,
        "chr_set": sorted(chromosomes)
    }
