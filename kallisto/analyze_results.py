#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, math, statistics as st
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def ensure_dir(p):
    if not os.path.exists(p):
        os.makedirs(p, exist_ok=True)

def load_abundance(path):
    df = pd.read_csv(path, sep="\t")
    # انتظار ستون‌ها: target_id, length, eff_length, est_counts, tpm
    need = {"target_id","tpm"}
    if not need.issubset(set(df.columns)):
        raise ValueError(f"file {path} missing required columns: {need}")
    return df[["target_id","length","eff_length","est_counts","tpm"]].copy()

def align_by_target(df_a, df_b, suffix_a="my", suffix_b="kal"):
    a = df_a.sort_values("target_id").reset_index(drop=True)
    b = df_b.sort_values("target_id").reset_index(drop=True)
    # inner join تا فقط مشترک‌ها مقایسه شوند
    m = pd.merge(a, b, on="target_id", suffixes=(f"_{suffix_a}", f"_{suffix_b}"))
    return m

def pearson(xs, ys):
    mx, my = st.mean(xs), st.mean(ys)
    num = sum((x-mx)*(y-my) for x,y in zip(xs,ys))
    den = (sum((x-mx)**2 for x in xs) * sum((y-my)**2 for y in ys))**0.5
    return (num/den) if den else 0.0

def spearman(xs, ys):
    # رتبه‌ها
    def rank(vs):
        order = np.argsort(vs)
        ranks = np.empty_like(order, dtype=float)
        ranks[order] = np.arange(1, len(vs)+1, dtype=float)
        return ranks
    rx, ry = rank(np.array(xs)), rank(np.array(ys))
    return pearson(rx, ry)

def safe_log(x, eps=1e-6):
    return math.log(max(x, eps))

def plot_scatter_loglog(xs, ys, out_png, title="TPM correlation", label_x="Your TPM", label_y="Kallisto TPM"):
    plt.figure(figsize=(6,6))
    plt.scatter(xs, ys, s=10, alpha=0.6)
    plt.xscale("log"); plt.yscale("log")
    lim_min = min([min([v for v in xs if v>0], default=1e-6),
                   min([v for v in ys if v>0], default=1e-6)])
    lim_max = max(max(xs), max(ys), 1.0)
    plt.plot([lim_min, lim_max],[lim_min, lim_max],'r--',lw=1,label="y=x")
    plt.xlabel(label_x); plt.ylabel(label_y)
    plt.title(title)
    plt.grid(ls="--", alpha=0.4)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_png, dpi=150)
    plt.close()

def plot_bland_altman(xs, ys, out_png):
    a = np.array(xs); b = np.array(ys)
    mean = (a+b)/2.0
    diff = a-b
    plt.figure(figsize=(6,4))
    plt.scatter(mean, diff, s=10, alpha=0.5)
    plt.axhline(0, color='red', ls='--')
    plt.xscale("log")
    plt.xlabel("Mean TPM (log scale)")
    plt.ylabel("Difference (Your - Kallisto)")
    plt.title("Bland–Altman: TPM Bias Check")
    plt.grid(ls="--", alpha=0.4)
    plt.tight_layout()
    plt.savefig(out_png, dpi=150)
    plt.close()

def plot_hist_log_tpm(df, out_png, title="Global Expression Distribution"):
    vals = df["tpm"].values
    plt.figure(figsize=(6,4))
    plt.hist(vals, bins=50, edgecolor="k")
    plt.xscale("log")
    plt.xlabel("TPM (log scale)")
    plt.ylabel("Transcript count")
    plt.title(title)
    plt.grid(ls="--", alpha=0.4)
    plt.tight_layout()
    plt.savefig(out_png, dpi=150)
    plt.close()

def plot_bootstrap_cv(df_bs, out_png):
    if "cv_tpm" not in df_bs.columns:
        raise ValueError("bootstrap summary lacks cv_tpm column.")
    vals = df_bs["cv_tpm"].values
    plt.figure(figsize=(6,4))
    plt.hist(vals, bins=30, edgecolor="k")
    plt.xlabel("CV (Coefficient of Variation)")
    plt.ylabel("Count")
    plt.title("Bootstrap Variability across Transcripts")
    plt.grid(ls="--", alpha=0.4)
    plt.tight_layout()
    plt.savefig(out_png, dpi=150)
    plt.close()

def top_spearman(df_aligned, topn=50):
    # topN بر اساس TPM تو
    d = df_aligned.sort_values("tpm_my", ascending=False).head(topn)
    return d["tpm_my"].corr(d["tpm_kal"], method="spearman")

def main():
    ap = argparse.ArgumentParser(description="Analyze quant outputs vs Kallisto")
    ap.add_argument("--mine", default="abundance.tsv", help="Your abundance.tsv")
    ap.add_argument("--kallisto", required=True, help="kallisto abundance.tsv")
    ap.add_argument("--bootstrap", default="abundance_bootstrap_summary.tsv", help="(optional) bootstrap summary path")
    ap.add_argument("--outdir", default="plots", help="output plots directory")
    ap.add_argument("--topn", type=int, default=50, help="Top-N for Spearman on top expressed transcripts")
    args = ap.parse_args()

    ensure_dir(args.outdir)

    # Load
    df_my  = load_abundance(args.mine)
    df_kal = load_abundance(args.kallisto)
    df_aln = align_by_target(df_my, df_kal, suffix_a="my", suffix_b="kal")

    # Correlations
    xs = df_aln["tpm_my"].tolist()
    ys = df_aln["tpm_kal"].tolist()
    p = pearson(xs, ys)
    s = spearman(xs, ys)
    s_top = top_spearman(df_aln, topn=args.topn)

    # Plots
    plot_scatter_loglog(xs, ys, os.path.join(args.outdir, "tpm_scatter.png"),
                        title=f"TPM correlation (Pearson={p:.3f})",
                        label_x="Your TPM", label_y="Kallisto TPM")
    plot_bland_altman(xs, ys, os.path.join(args.outdir, "tpm_bland_altman.png"))
    plot_hist_log_tpm(df_my, os.path.join(args.outdir, "tpm_hist_mine.png"),
                      title="Global Expression Distribution (Your quant)")

    # Bootstrap (optional)
    bs_done = False
    if os.path.exists(args.bootstrap):
        try:
            df_bs = pd.read_csv(args.bootstrap, sep="\t")
            plot_bootstrap_cv(df_bs, os.path.join(args.outdir, "bootstrap_cv_hist.png"))
            bs_done = True
        except Exception as e:
            print(f"[warn] bootstrap summary not parsed: {e}")

    # Save aligned table (for reproducibility)
    df_aln_out = df_aln[["target_id","tpm_my","tpm_kal"]].copy()
    df_aln_out.to_csv(os.path.join(args.outdir, "aligned_tpm.tsv"), sep="\t", index=False)

    # Summary
    with open(os.path.join(args.outdir, "analysis_summary.txt"), "w") as w:
        w.write("=== Quant Analysis Summary ===\n")
        w.write(f"Your file: {args.mine}\n")
        w.write(f"Kallisto : {args.kallisto}\n")
        w.write(f"Aligned transcripts: {len(df_aln)}\n\n")
        w.write(f"Pearson (all):  {p:.4f}\n")
        w.write(f"Spearman (all): {s:.4f}\n")
        w.write(f"Spearman (top {args.topn} by your TPM): {s_top:.4f}\n")
        w.write("\nGenerated plots:\n")
        w.write(" - plots/tpm_scatter.png\n")
        w.write(" - plots/tpm_bland_altman.png\n")
        w.write(" - plots/tpm_hist_mine.png\n")
        if bs_done:
            w.write(" - plots/bootstrap_cv_hist.png\n")
        else:
            w.write(" - (bootstrap plot skipped: file not found or parse error)\n")
        w.write("\nAligned TPM table: plots/aligned_tpm.tsv\n")
    print("✅ Done. See:", args.outdir)

if __name__ == "__main__":
    main()
