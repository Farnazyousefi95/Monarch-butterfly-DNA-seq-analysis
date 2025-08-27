# SRR25297534 — Monarch WGS (Coverage & SNVs)

**Question:** What is the average genomic coverage and what SNVs are present compared to the Dplex_v4 reference?  
**Mean depth:** 17.092× · **Breadth ≥1×:** 0.9498 · **Breadth ≥10×:** 0.7836  
**Variants (filtered):** total 6175516; SNPs 5635808; indels 539708; Ti/Tv 1.00

## Reproduce
```bash
conda env create -f env_dplex.yml
conda activate dplex
```
See **results/** for small artifacts and **logs/** for run logs (index in logs/LOG_INDEX.md).  
Large files (BAM/FASTQ/VCF) intentionally excluded from Git.
