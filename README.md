# SRR25297534 — Monarch Butterfly (Danaus plexippus) WGS: Coverage & SNVs

**Group 3 (Week 2)** · Melika · Kiran · Farnaz

## Goal (plain English)
1) **Coverage:** How completely and how deeply we read the genome.  
2) **Variants:** Where this butterfly’s DNA differs by a single letter (SNVs) from the reference.

## Dataset
- SRA run: **SRR25297534** (Illumina paired-end WGS)

## Workflow (done)
Download → QC → Trim → Re-QC → Reference → **Align** → **Coverage** → Variants

## Headline results
- **Alignment:** 97.96% mapped; 94.72% properly paired (see `results/flagstat_head.txt`)
- **Coverage:** mean depth **17.092×**; breadth ≥1× **0.9498**; breadth ≥10× **0.7836**
- **Variants (filtered):** total **6175516**; SNPs **5635808**; indels **539708**; Ti/Tv **1.00**

## Reproduce (conda env)
```bash
conda env create -f env_dplex.yml
conda activate dplex
```

## Small artifacts for slides
- `results/coverage_summary.txt`, `results/coverage_table_head.tsv`
- `results/variants_summary.csv`, `results/variants_summary.md`
- `results/fastqc_pretrim_summary.csv`, `results/fastqc_posttrim_summary.csv` (if present)
- alignment heads: `results/flagstat_head.txt`, `results/idxstats_head.tsv`

## Notes
- BAM/FASTQ/VCF kept out of Git (see `.gitignore`); this repo documents the workflow and results.
