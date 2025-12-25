
from pathlib import Path
import json
from parsers.fasta_parser import scan_fasta
from parsers.gff_parser import scan_gff
from parsers.vcf_parser import scan_vcf
from validators.chr_consistency import check_chr_consistency

def run_species(species_dir, out_dir, log_dir):
    species = species_dir.name
    fasta_info = scan_fasta(species_dir / "fasta")
    gff_info = scan_gff(species_dir / "gff")
    vcf_info = scan_vcf(species_dir / "vcf")
    check_chr_consistency(fasta_info, gff_info, vcf_info, species)
    Path(out_dir).mkdir(exist_ok=True)
    with open(Path(out_dir) / f"{species}.json", "w") as f:
        json.dump({
            "species": species,
            "genome": fasta_info,
            "genes": gff_info,
            "variants": vcf_info
        }, f, indent=2)
