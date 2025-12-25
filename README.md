
# Genome Statistics Project

This project performs hierarchical genome statistics per species using
FASTA / GFF / VCF files, following Ensembl-style organization.

Run:
python src/main.py --data-root data/species --out stats --log-dir logs



'''text
genome_stats_full/
├── README.md
│
├── data/
│   └── species/
│       └── <species_name>/
│           ├── fasta/
│           ├── gff/
│           ├── vcf/
│           └── meta.yaml
│
├── src/
│   ├── main.py
│   │
│   ├── parsers/
│   │   ├── fasta_parser.py     # FASTA（支持按染色体拆分）
│   │   ├── gff_parser.py       # GFF gene/exon/CDS/UTR 结构解析
│   │   └── vcf_parser.py       # VCF 变异位置解析
│   │
│   ├── stats/
│   │   └── aggregator.py       # 层级统计树构建（核心）
│   │
│   ├── validators/
│   │   ├── schema_validate.py  # JSON Schema 校验（结构级）
│   │   └── consistency_validate.py # 数值/逻辑一致性校验
│   │
│   └── schema/
│       └── genome.schema.json  # 统计结果标准定义
│
├── stats/   # 自动生成统计结果
└── logs/
