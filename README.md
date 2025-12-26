# Genome Statistics – Fully Recursive (Flat Layout, Strict Validation)

## Data Layout (no fasta/gff/vcf subfolders)
以下是 genome_stats_full 项目的完整目录结构及文件说明：
```plaintext
genome_stats_recursive/
├── README.md
├── src/
│  ├── main.py
│   ├── build_tree.py
│   ├── validate_output.py
│   ├── tree/
│   │   ├── __init__.py
│   │   └── node.py
│   ├── parsers/
│   │   ├── __init__.py
│   │   ├── fasta_scan.py
│   │   ├── gff_ensembl.py
│   │   └── vcf_stream.py
│   ├── validators/
│   │   ├── __init__.py
│   │   ├── structural.py
│   │   └── consistency.py
│   └── utils/
│       ├── __init__.py
│       ├── discovery.py
│       ├── intervals.py
│       └── naming.py
│
├──data/
    └── <species>/
        └── <assembly>/
            ├── *.fa / *.fa.gz
            ├── *.gff / *.gff3 / *.gff.gz
            └── *.vcf / *.vcf.gz

```

## Output Tree (per species JSON)
fasta+gff+vcf统计信息结果
```plaintext
Species
└── Assembly
    └── Chromosome
        ├── gene_union_group
        │   └── gene
        │       ├── exon_group
        │       │   └── exon_segment
        │       │       ├── cds_group -> cds_segment*
        │       │       ├── utr_group -> utr_segment*
        │       │       └── exon_other_group -> exon_other_segment*
        │       └── non_exon_group -> intron_segment*
        └── non_gene_group -> non_gene_segment*

Each genomic node includes:
- start/end/length
- variant_count
- missing_N_count + missing_N_runs (from FASTA)

Two validators are included:
- structural (schema-like recursive checks)
- consistency (bounds, containment, partitioning)

## Run

python src/main.py --data-root /path/to/data --out stats
python src/validate_output.py --json stats/<species>.json
