import os, re, pathlib
import pandas as pd
import matplotlib.pyplot as plt

USER = os.environ.get("USER","")
PROJ = f"/scratch/{USER}/week2_assignment"
FIGDIR = pathlib.Path("docs/figs")
FIGDIR.mkdir(parents=True, exist_ok=True)

# ---------- Coverage: per-contig (top 20 by mean depth) ----------
cov_path = f"{PROJ}/coverage_table.tsv"
if os.path.exists(cov_path):
    cov = pd.read_csv(cov_path, sep="\t", header=0)
    if '#rname' in cov.columns:
        cov = cov.rename(columns={'#rname':'rname'})
    cols_needed = {'rname','meandepth','coverage','meanbaseq','meanmapq'}
    miss = cols_needed - set(cov.columns)
    if miss:
        print("[WARN] coverage_table.tsv missing columns:", miss)
    cov_top = cov.sort_values('meandepth', ascending=False).head(20)
    plt.figure(figsize=(10,4.5))
    plt.bar(range(len(cov_top)), cov_top['meandepth'])
    plt.xticks(range(len(cov_top)), cov_top['rname'], rotation=90)
    plt.ylabel("Mean depth (×)")
    plt.title("Top 20 contigs by mean depth")
    plt.tight_layout()
    plt.savefig(FIGDIR/"coverage_top20_meandepth.png", dpi=200); plt.close()

# ---------- Depth histogram (0..100, with 101 = 100+) ----------
dh_path = f"{PROJ}/stats/depth_hist.tsv"
if os.path.exists(dh_path):
    dh = pd.read_csv(dh_path, sep="\t", header=None, names=["depth","count"])
    plt.figure(figsize=(8,4))
    plt.plot(dh['depth'].to_numpy(), dh['count'].to_numpy())
    plt.xlabel("Depth (×)  [101 = 100+]")
    plt.ylabel("# genomic positions")
    plt.title("Genome-wide depth histogram")
    plt.tight_layout()
    plt.savefig(FIGDIR/"depth_histogram.png", dpi=200); plt.close()

# ---------- Insert size distribution (samtools stats) ----------
ss_path = f"{PROJ}/stats/samstats.txt"
if os.path.exists(ss_path):
    is_rows = []
    with open(ss_path) as fh:
        for line in fh:
            if line.startswith("IS\t"):
                # Format: IS\t<insert>\t<count>\t...
                parts = line.strip().split("\t")
                if len(parts) >= 3:
                    try:
                        ins = int(parts[1]); cnt = int(parts[2])
                        is_rows.append((ins, cnt))
                    except: pass
    if is_rows:
        import pandas as pd
        is_df = pd.DataFrame(is_rows, columns=["insert","count"]).sort_values("insert")
        # Limit very large inserts for visualization (e.g., <= 1000)
        is_df = is_df[is_df["insert"]<=1000]
        plt.figure(figsize=(8,4))
        plt.plot(is_df["insert"].to_numpy(), is_df["count"].to_numpy())
        plt.xlabel("Insert size (bp)")
        plt.ylabel("Count")
        plt.title("Insert size distribution (truncated to 1000 bp)")
        plt.tight_layout()
        plt.savefig(FIGDIR/"insert_size_distribution.png", dpi=200); plt.close()

# ---------- Variants: counts & Ti/Tv from bcftools stats ----------
vs_path = f"{PROJ}/variants.stats.txt"
if os.path.exists(vs_path):
    # Pull counts
    sn = {}
    with open(vs_path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if parts[0] == "SN":
                sn[parts[2]] = parts[3]
    snp = int(sn.get("number of SNPs:", 0))
    indel = int(sn.get("number of indels:", 0))
    # TSTV last line
    titv = None
    with open(vs_path) as fh:
        for line in fh:
            if line.startswith("TSTV"):
                titv = line.rstrip("\n").split("\t")
    if titv and len(titv)>=6:
        ti = int(titv[2]); tv = int(titv[3]); ratio = float(titv[4])
        # Variant type bar
        plt.figure(figsize=(5,3.5))
        plt.bar(["SNP","INDEL"], [snp, indel])
        plt.ylabel("Count")
        plt.title("Variant counts (filtered)")
        plt.tight_layout()
        plt.savefig(FIGDIR/"variant_counts.png", dpi=200); plt.close()

        # Ti/Tv bar
        plt.figure(figsize=(5,3.5))
        plt.bar(["Ti","Tv"], [ti, tv])
        plt.title(f"Ti/Tv ratio ≈ {ratio:.2f}")
        plt.tight_layout()
        plt.savefig(FIGDIR/"titv_bar.png", dpi=200); plt.close()

# ---------- FastQC before vs after (PASS/WARN/FAIL counts) ----------
def status_counts(df, colname):
    c = df[colname].value_counts()
    return [c.get("PASS",0), c.get("WARN",0), c.get("FAIL",0)]

pre_csv = f"{PROJ}/qc/fastqc/summary.csv"
post_csv = f"{PROJ}/qc/fastqc_posttrim/summary_posttrim.csv"
if os.path.exists(pre_csv) and os.path.exists(post_csv):
    pre = pd.read_csv(pre_csv)
    post = pd.read_csv(post_csv)
    # Try to infer columns
    # Pre has columns: Module, SRR25297534_1, SRR25297534_2 (as created earlier)
    pre_cols = [c for c in pre.columns if c != "Module"]
    post_cols = [c for c in post.columns if c != "Module"]
    if len(pre_cols)>=2 and len(post_cols)>=2:
        pre_r1, pre_r2 = pre_cols[:2]
        post_r1, post_r2 = post_cols[:2]
        labels = ["PASS","WARN","FAIL"]

        for tag, pcol, qcol, fname in [
            ("R1", pre_r1, post_r1, "fastqc_before_after_R1.png"),
            ("R2", pre_r2, post_r2, "fastqc_before_after_R2.png")
        ]:
            p = status_counts(pre, pcol)
            q = status_counts(post, qcol)
            import numpy as np
            x = np.arange(len(labels)); w = 0.38
            plt.figure(figsize=(6,3.8))
            plt.bar(x - w/2, p, width=w, label="Before")
            plt.bar(x + w/2, q, width=w, label="After")
            plt.xticks(x, labels)
            plt.ylabel("Modules")
            plt.title(f"FastQC {tag}: Before vs After Trimming")
            plt.legend()
            plt.tight_layout()
            plt.savefig(FIGDIR/fname, dpi=200); plt.close()

print("[OK] Figures written to", FIGDIR)
