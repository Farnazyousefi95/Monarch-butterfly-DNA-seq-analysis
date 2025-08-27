# Monarch-butterfly-DNA-seq-analysis.
Group3_Week2_Assignment


To submit jobs: qsub download.pbs


# SRR25297534 — Monarch Butterfly (Danaus plexippus) WGS: Coverage & SNVs

**Group 3 (Week 2, Computational Biology Colloquium)**  
**Team:** Farnaz · Kiran ·   
**Focus:** End-to-end, single-sample DNA-seq analysis on **SRR25297534** (Illumina PE WGS) — from download to QC → trimming → alignment → coverage (variants next).

---

## 🔎 Project Goal

**Question:** *What is the average genomic coverage and what notable single nucleotide variants (SNVs) are present in **SRR25297534** compared to a monarch butterfly reference genome?*

Why this is sensible: the dataset is **DNA‑seq** (not RNA‑seq), so coverage and variant discovery are the most direct, low‑assumption analyses.

---

## 📦 Dataset

| Field | Value |
|---|---|
| **SRA Run** | **SRR25297534** |
| **Organism** | *Danaus plexippus* (monarch butterfly) |
| **Assay / Platform** | Whole‑genome sequencing (Illumina HiSeq X Ten) |
| **Layout** | Paired-end |
| **Size (approx.)** | ~16.5 M read pairs; ~4.9 Gbases |

> We analyze **SRR25297534** throughout this repository.

---

## 🧭 Workflow Overview (Week‑2 Scope)

Download → **QC** → **Trim** → **Re‑QC** → **Reference** → **Alignment** → **Coverage** → (Variants & visualization next)

```
raw FASTQ ──FastQC──► trim (Trimmomatic) ──FastQC──►
                        ▼
                  ref (Dplex_v4)  ◄── bwa index + faidx
                        ▼
                bwa mem → samtools sort/index
                        ▼
          coverage (mean depth, breadth, per‑contig)
                        ▼
                 [ next: bcftools variants + Ti/Tv ]
```

---

## ✅ Status & Headline Results (so far)

- **Pre‑trim QC**: typical WGS—quality PASS; composition/duplication flags common pre‑alignment.
- **Trimming (Trimmomatic)** + **post‑trim QC**: quality remains strong; composition/duplication flags persist (expected before mapping).
- **Reference**: *Danaus plexippus* **Dplex_v4** (“toplevel”) downloaded and indexed (BWA + samtools).
- **Alignment (BWA‑MEM)**: **31,818,650** reads total; **97.96%** mapped; **94.72%** properly paired.
- **Coverage (samtools)**: **mean depth = 17.092×**, breadth ≥1× **94.98%**, breadth ≥10× **78.36%**.

> These values come from the files listed below; see **Reproduce the Results** for exact commands.

---

## 📁 Repository Layout (suggested)

```
repo/
├── README.md                      # this file
├── scripts/                       # PBS/SLURM/CLI scripts (portable where possible)
│   ├── 01_download_sra.pbs
│   ├── 02_fastqc_pretrim.pbs
│   ├── 03_trim_trimmomatic_postqc.pbs
│   ├── 04_get_reference_dplex_v4.pbs
│   ├── 05_align_bwa_mem.pbs
│   ├── 06_coverage.pbs
│   └── 07_variants_bcftools.pbs   # (next step)
├── env/                           # environment files (conda YAMLs)
│   └── env_dplex.yml
├── results/                       # small text outputs for documentation
│   ├── flagstat_head.txt
│   ├── mean_depth.txt
│   ├── breadth.txt
│   └── coverage_table_head.tsv
└── docs/                          # figures/tables for slides (optional)
    ├── fastqc_before_after_R1.png
    ├── fastqc_before_after_R2.png
    └── coverage_summary.md
```

> On HPC, we executed under `/scratch/<user>/week2_assignment/`. In this Git repo, we include the **scripts** and **small text outputs** only (avoid uploading large BAM/FASTQ).

---

## 🔧 Environment

We used a small **conda** environment (works on Linux/macOS):

```bash
conda create -n dplex -c conda-forge -c bioconda \
  bwa=0.7.17 samtools=1.17 bcftools=1.17 trimmomatic=0.39 fastqc openjdk
conda activate dplex
```

Notes:
- **FastQC** requires Java; the `openjdk` package in the env provides it.
- If your cluster provides modules, you can mix-and-match (but keep versions similar).

---

## 🚀 Reproduce the Results

All commands assume the working directory contains **SRR25297534** FASTQs and the **Dplex_v4** reference after download. On HPC, adapt headers to your scheduler; below is CLI‑style for clarity.

### 1) Download SRA → FASTQ
```bash
# Using sratoolkit
prefetch SRR25297534
fasterq-dump SRR25297534 --split-files -O .
```

### 2) FastQC (pre‑trim)
```bash
fastqc -t 4 SRR25297534_1.fastq SRR25297534_2.fastq -o qc/fastqc
```

### 3) Trimming (Trimmomatic) + FastQC (post‑trim)
```bash
trimmomatic PE -threads 4 -phred33 \
  SRR25297534_1.fastq SRR25297534_2.fastq \
  trim/SRR25297534_1P.trim.fastq.gz trim/SRR25297534_1U.trim.fastq.gz \
  trim/SRR25297534_2P.trim.fastq.gz trim/SRR25297534_2U.trim.fastq.gz \
  ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50

fastqc -t 4 trim/SRR25297534_1P.trim.fastq.gz trim/SRR25297534_2P.trim.fastq.gz -o qc/fastqc_posttrim
```

### 4) Reference (Dplex_v4) download + indexing
```bash
# Download "toplevel" FASTA from Ensembl Metazoa
curl -L -o ref/Danaus_plexippus.Dplex_v4.dna.toplevel.fa.gz "https://<ensembl_metazoa_url>/Danaus_plexippus.Dplex_v4.dna.toplevel.fa.gz"
gunzip ref/Danaus_plexippus.Dplex_v4.dna.toplevel.fa.gz

# Index
samtools faidx ref/Danaus_plexippus.Dplex_v4.dna.toplevel.fa
bwa index ref/Danaus_plexippus.Dplex_v4.dna.toplevel.fa
```

### 5) Alignment (BWA‑MEM) → sorted BAM + QC
```bash
RG='@RG\tID:SRR25297534\tSM:SRR25297534\tPL:ILLUMINA\tLB:lib1\tPU:unit1'
bwa mem -t 8 -R "$RG" ref/Danaus_plexippus.Dplex_v4.dna.toplevel.fa \
  trim/SRR25297534_1P.trim.fastq.gz trim/SRR25297534_2P.trim.fastq.gz \
| samtools sort -@8 -o SRR25297534.dplexv4.sorted.bam -

samtools index SRR25297534.dplexv4.sorted.bam
samtools flagstat SRR25297534.dplexv4.sorted.bam > SRR25297534.flagstat.txt
samtools idxstats SRR25297534.dplexv4.sorted.bam > SRR25297534.idxstats.tsv
```

**Observed from our run** (for your reference):
- Total reads: **31,818,650**
- Mapped: **97.96%** (31,169,855)
- Properly paired: **94.72%** (29,564,814)

### 6) Coverage (mean depth & breadth)
```bash
samtools coverage -o coverage_table.tsv SRR25297534.dplexv4.sorted.bam
samtools depth -a SRR25297534.dplexv4.sorted.bam | \
  awk '{s+=$3;n++}END{print s/n}' > mean_depth.txt

# breadth at ≥1x and ≥10x (optional)
samtools depth -a SRR25297534.dplexv4.sorted.bam | \
  awk 'BEGIN{c1=0;c10=0;N=0}{N++; if($3>0)c1++; if($3>=10)c10++} \
       END{printf("breadth>=1x\t%.4f\nbreadth>=10x\t%.4f\n", c1/N, c10/N)}' > breadth.txt
```

**Observed from our run:**
- **Mean depth:** **17.092×**
- **Breadth ≥1×:** **0.9498** (94.98%)
- **Breadth ≥10×:** **0.7836** (78.36%)

### 7) Variants (next step)
```bash
# Install once in the env if needed
# conda install -y -c bioconda bcftools=1.17

bcftools mpileup -Ou -f ref/Danaus_plexippus.Dplex_v4.dna.toplevel.fa SRR25297534.dplexv4.sorted.bam \
 | bcftools call -mv -Ou \
 | bcftools view -Oz -o variants.raw.vcf.gz
tabix -p vcf variants.raw.vcf.gz

bcftools filter -e 'QUAL<30 || DP<5' -S . variants.raw.vcf.gz -Oz -o variants.filt.vcf.gz
tabix -p vcf variants.filt.vcf.gz

bcftools stats -F ref/Danaus_plexippus.Dplex_v4.dna.toplevel.fa -s - variants.filt.vcf.gz > variants.stats.txt
```

**Report next:** total SNPs/indels, Ti/Tv, SNV class counts.

---

## 📊 What to Put in a Slide Deck

- **QC (pre vs post):** PASS/WARN/FAIL counts; note that composition/duplication flags pre‑alignment are expected.
- **Alignment metrics:** mapping %, properly paired %, small note on insert size (~200–230 bp median FR).
- **Coverage:** mean depth; breadth at ≥1× and ≥10×; 1–2 lines from `coverage_table.tsv`.
- **Variants (when done):** SNP/indel counts and Ti/Tv from `variants.stats.txt`.

---

## 🛠️ Troubleshooting (things we hit)

- **`LC_ALL: unbound variable`** during PBS prolog: avoid `set -u` while sourcing `/etc/profile`; set `LC_ALL/LANG` to `C` up front.
- **FastQC needs Java:** ensure `openjdk` is in the conda env or the system `java` is in `PATH`.
- **Queue minima:** our `medium` queue required `walltime ≥ 01:00:00`.
- **FastQC HTML inside ZIP:** extract `fastqc_report.html` from `*_fastqc.zip` if no `*_fastqc.html` was written.
- **No cluster modules for bwa/samtools:** use a conda env (`dplex`) with those tools installed.
- **SCP from the **right** machine:** run `scp` from your **local terminal**, not from inside the HPC shell.

---

## 🧪 Reproducibility

- Scripts checked into `scripts/` and referenced in this README.
- Deterministic toolchain; no random seeding involved.
- Keep large binaries (BAM/VCF/FASTQ) **out of Git**; archive small text results in `results/` for documentation.

---

## 📌 Roadmap

- [x] Download → QC → Trim → Re‑QC  
- [x] Reference + indexing  
- [x] Alignment + mapping QC  
- [x] Coverage (mean depth, breadth)  
- [x] Variants (SNVs/indels) + Ti/Tv  
- [ ] Optional: mark duplicates, dedup coverage, contamination screen, functional annotation

---

## 🙌 Acknowledgements

- Alabama Supercomputer Center for compute resources.  
- Open‑source tools: sratoolkit, FastQC, Trimmomatic, BWA, samtools, bcftools.


---

## ✉️ Contact

For questions about this analysis: open an issue or contact the Group 3 maintainers.

