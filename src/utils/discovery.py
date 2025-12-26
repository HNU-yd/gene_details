from dataclasses import dataclass
from pathlib import Path
from typing import List

FA_EXT = (".fa", ".fasta", ".fa.gz", ".fasta.gz")
GFF_EXT = (".gff", ".gff3", ".gff.gz", ".gff3.gz")
VCF_EXT = (".vcf", ".vcf.gz")

@dataclass
class AssemblyFiles:
    fasta: List[Path]
    gff: List[Path]
    vcf: List[Path]

def discover_assembly_files(assembly_dir: Path) -> AssemblyFiles:
    fasta, gff, vcf = [], [], []
    for p in sorted(assembly_dir.iterdir()):
        if not p.is_file():
            continue
        name = p.name.lower()
        if any(name.endswith(ext) for ext in FA_EXT):
            fasta.append(p)
        elif any(name.endswith(ext) for ext in GFF_EXT):
            gff.append(p)
        elif any(name.endswith(ext) for ext in VCF_EXT):
            vcf.append(p)
    return AssemblyFiles(fasta=fasta, gff=gff, vcf=vcf)
