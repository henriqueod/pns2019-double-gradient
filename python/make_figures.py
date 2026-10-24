"""Figures 1, 2 and S1 of the manuscript, read from the CSVs written by the R scripts.

Usage (from the repository root, after running the R pipeline):
    python python/make_figures.py [results_dir] [figures_dir]

Defaults: results_dir = ./results, figures_dir = ./figures.

Inputs
  Figure 1: results/anymh/03b_prevalencia_por_quintil.csv, 07b_PR_quintis_blocoB.csv,
            08_forma_do_gradiente.csv                                  (R/02_main_analysis.R)
  Figure 2: results/anymh/06_blocos_coorte.csv                            (R/02_main_analysis.R)
            results/depression_cascade_domain_rank/*_PR_income.csv       (R/04_depression_cascade.R)
            results/detection_gap_domain_rank/03_lacuna_RII_renda.csv    (R/05_detection_gap.R)
            results/referral_failure_blockD/33_falha_encaminhamento_modelos.csv
                                                                         (R/08_referral_failure_blockD.R)
  Figure S1: results/forest_anymh/50c_importancia_media_5_particoes_anymh.csv (R/07_random_forest.R)
"""
import os
import sys

import logging
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D

RES = sys.argv[1] if len(sys.argv) > 1 else "results"
FIGDIR = sys.argv[2] if len(sys.argv) > 2 else "figures"
os.makedirs(FIGDIR, exist_ok=True)
# Carlito is metrically compatible with Calibri (the manuscript font); matplotlib falls back
# to the next available family.
logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
plt.rcParams.update({"font.family": ["Carlito", "Calibri", "DejaVu Sans"], "font.size": 10})

BLUE, RED = "#2f6f9f", "#b5452f"


def rd(*parts):
    return pd.read_csv(os.path.join(RES, *parts))


# ---------------------------------------------------------------------------
# Figure 1. Prevalence by quintile of the income rank, with block-B PRs
# ---------------------------------------------------------------------------
def figure1():
    q = rd("anymh", "03b_prevalencia_por_quintil.csv")
    pr = rd("anymh", "07b_PR_quintis_blocoB.csv")
    shape = rd("anymh", "08_forma_do_gradiente.csv")
    labs = {"psychotherapy": "A. Current psychotherapy", "high_need_gap": "B. High need gap"}
    cols = {"psychotherapy": BLUE, "high_need_gap": RED}
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.8), sharey=False)
    for ax, (k, t) in zip(axes, labs.items()):
        d = q[q.desfecho == k].sort_values("income5").reset_index(drop=True)
        x = range(1, 6)
        y = d.prevalencia * 100
        lo = d.ci_l * 100
        hi = d.ci_u * 100
        ax.errorbar(x, y, yerr=[y - lo, hi - y], fmt="o-", color=cols[k], capsize=4, lw=1.6, ms=6)
        p = pr[pr.desfecho == k].sort_values("termo").reset_index(drop=True)
        for i in range(5):
            if i == 0:
                s = "ref."
            else:
                r = p.iloc[i - 1]
                s = f"PR {r.est:.2f}\n({r.lower:.2f}–{r.upper:.2f})"
            ax.annotate(s, (i + 1, hi[i]), textcoords="offset points", xytext=(0, 6), ha="center",
                        fontsize=7.5, color="#333")
        ax.set_xticks(list(x))
        ax.set_xticklabels(["Q1\n(poorest)", "Q2", "Q3", "Q4", "Q5\n(richest)"])
        ax.set_title(t, loc="left", fontsize=11, fontweight="bold")
        ax.set_ylabel("Weighted prevalence (%)")
        ax.set_ylim(0, max(hi) * 1.28)
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", color="#e5e5e5", lw=0.6)
        ax.set_axisbelow(True)
        # Departure from linearity: quintile term conditional on the linear rank (income5 = quintiles).
        pv = shape[(shape.desfecho == k) & (shape.resolucao == "income5")
                   & shape.teste.str.startswith("desvio")].p.iloc[0]
        ptxt = "p < 0.001" if pv < 0.001 else (f"p = {pv:.3f}" if pv < 0.01 else f"p = {pv:.2f}")
        ax.text(0.02, 0.03, f"Departure from linearity {ptxt}", transform=ax.transAxes, fontsize=8,
                color="#555")
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "fig1_shape.png"), dpi=300)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 2. Forest plot of RIIs across care indicators (minimal vs adjusted model)
# ---------------------------------------------------------------------------
def two_models(df, min_label, adj_label, col_model="model", est="PR", lo="lower", hi="upper"):
    a = df[df[col_model] == min_label].iloc[0]
    b = df[df[col_model] == adj_label].iloc[0]
    return (a[est], a[lo], a[hi], b[est], b[lo], b[hi], int(a["n_domain"]) if "n_domain" in df else None)


def figure2():
    blk = rd("anymh", "06_blocos_coorte.csv")
    blk = blk[blk.analise == "coorte"]

    def cohort(outcome):
        b = blk[(blk.desfecho == outcome) & (blk.bloco == "B")].iloc[0]
        d = blk[(blk.desfecho == outcome) & (blk.bloco == "D")].iloc[0]
        return (b.RII, b.RII_lower, b.RII_upper, d.RII, d.RII_lower, d.RII_upper, int(b.n))

    cas = "depression_cascade_domain_rank"
    spec = two_models(rd(cas, "specialist_reached_PR_income.csv"), "minimum", "expanded")
    hid = rd(cas, "hidden_remission_gap_PR_income.csv")
    hid_sev = two_models(hid[hid.variable == "hidden_remission_severe"], "minimum", "expanded")
    hid_sh = two_models(hid[hid.variable == "hidden_remission_selfharm"], "minimum", "expanded")
    hid_any = two_models(hid[hid.variable == "hidden_remission_gap"], "minimum", "expanded")
    struct = two_models(rd(cas, "structural_barrier_PR_income.csv"), "minimum", "expanded")
    crit = two_models(rd(cas, "critical_gap_PR_income.csv"), "minimum", "expanded")
    # Referral failure comes from R/08: minimal model (M0) and block-D model (M1), the same
    # covariates as the other adjusted rows. The R/04 "expanded" model also adjusted for
    # public-sector care and is now a sensitivity model (M2 in R/08).
    r8 = rd("referral_failure_blockD", "33_falha_encaminhamento_modelos.csv")
    m0 = r8[r8.modelo == "M0_minimo"].iloc[0]
    m1 = r8[r8.modelo == "M1_blocoD"].iloc[0]
    ref = (m0.RII, m0.lower, m0.upper, m1.RII, m1.lower, m1.upper, int(m0.n))
    det = rd("detection_gap_domain_rank", "03_lacuna_RII_renda.csv")
    det = det[(det.dominio == "D1_todos_phq10") & (det.desfecho == "undiagnosed") & (det["rank"] == "dominio")]
    dmin = det[det.modelo == "minimo"].iloc[0]
    dexp = det[det.modelo == "expandido"].iloc[0]
    undiag = (dmin.RII, dmin.lower, dmin.upper, dexp.RII, dexp.lower, dexp.upper, int(dmin.n))

    def n(v):
        return f"{v:,}"

    psy, hng = cohort("psychotherapy"), cohort("high_need_gap")
    rows = [
        ("A", "Current psychotherapy", f"Any diagnosed disorder (n={n(psy[6])})", *psy[:6]),
        ("A", "Reached a mental health specialist", f"Diagnosed depression, recent care (n={n(spec[6])})", *spec[:6]),
        ("A", "Undiagnosed depression (RII < 1 = more\nfrequent among the poorest)", f"PHQ-9 ≥10 (n={n(undiag[6])})", *undiag[:6]),
        ("B", "High need gap", f"Any diagnosed disorder (n={n(hng[6])})", *hng[:6]),
        ("B", "Perceived remission, PHQ-9 ≥15, untreated", f"Perceived remission (n={n(hid_sev[6])})", *hid_sev[:6]),
        ("B", "Perceived remission, thoughts of death, untreated", f"Perceived remission (n={n(hid_sh[6])})", *hid_sh[:6]),
        ("B", "Structural barrier as main reason", f"Not in regular follow-up (n={n(struct[6])})", *struct[:6]),
        ("C", "Perceived remission, PHQ-9 ≥10, untreated", f"Perceived remission (n={n(hid_any[6])})", *hid_any[:6]),
        ("C", "Untreated", f"Diagnosed, PHQ-9 ≥10, thoughts of death (n={n(crit[6])})", *crit[:6]),
        ("C", "No specialist visit after referral", f"Referred (n={n(ref[6])})", *ref[:6]),
    ]
    fig, ax = plt.subplots(figsize=(9, 6.6))
    ypos, labels, heads = [], [], []
    y = 0
    pat = {"A": "Pattern A — access to care resources", "B": "Pattern B — unmet need", "C": "No clear gradient"}
    cur = None
    for r in rows:
        if r[0] != cur:
            cur = r[0]
            y -= 0.4
            heads.append((y, pat[cur]))
            y -= 0.9
        ypos.append(y)
        labels.append(r)
        y -= 1.0
    for (yy, t) in heads:
        ax.text(0.052, yy, t, fontsize=10, fontweight="bold", va="center", transform=ax.get_yaxis_transform())
    for yy, r in zip(ypos, labels):
        ax.errorbar(r[3], yy + 0.17, xerr=[[r[3] - r[4]], [r[5] - r[3]]], fmt="o", color=BLUE, capsize=3, ms=5, lw=1.3)
        ax.errorbar(r[6], yy - 0.17, xerr=[[r[6] - r[7]], [r[8] - r[6]]], fmt="s", color=RED, capsize=3, ms=5,
                    lw=1.3, mfc="white")
    ax.axvline(1, color="#444", lw=0.8, ls="--")
    ax.set_xscale("log")
    ax.set_xticks([0.1, 0.2, 0.5, 1, 2, 4])
    ax.set_xticklabels(["0.1", "0.2", "0.5", "1", "2", "4"])
    ax.set_xlim(0.07, 4.5)
    ax.set_yticks(ypos)
    ax.set_yticklabels([f"{r[1]}\n{r[2]}" for r in labels], fontsize=8.5)
    ax.set_ylim(min(ypos) - 0.8, 0.2)
    ax.set_xlabel("Relative index of inequality (richest vs poorest; log scale)")
    ax.spines[["top", "right"]].set_visible(False)
    h = [Line2D([], [], marker="o", color=BLUE, ls="", label="Minimal model (age, sex, region)"),
         Line2D([], [], marker="s", color=RED, mfc="white", ls="",
                label="Adjusted model (+ race, education, urban, health plan)")]
    ax.legend(handles=h, loc="upper center", bbox_to_anchor=(0.45, -0.09), ncol=2, fontsize=8.5, frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "fig2_forest.png"), dpi=300)
    plt.close(fig)
    return rows


# ---------------------------------------------------------------------------
# Figure S1. Permutation importance, random forest (mean ± SD over five holdouts)
# ---------------------------------------------------------------------------
def figure_s1():
    d = rd("forest_anymh", "50c_importancia_media_5_particoes_anymh.csv")
    lab = {"income_rank_mental": "Income rank", "education3": "Education",
           "multiple_mental_dx": "Multiple mental disorder diagnoses", "age": "Age",
           "self_report": "Self-respondent (vs. proxy)", "depression": "Depression diagnosis",
           "schizophrenia": "Schizophrenia diagnosis", "region": "Region",
           "other_mental": "Other mental illness diagnosis", "family_support0": "Lack of family support",
           "friend_support0": "Lack of friend support", "bipolar": "Bipolar disorder diagnosis",
           "race3": "Race/skin color", "sex": "Sex", "urban": "Urban residence", "ocd": "OCD diagnosis"}
    d = d.iloc[::-1]
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.barh([lab.get(v, v) for v in d.variable], d.importance * 1000, xerr=d.sd * 1000, color=RED, alpha=.85,
            error_kw=dict(ecolor="#333", capsize=3, lw=1))
    ax.set_xlabel("Permutation importance (×10⁻³), mean ± SD across five holdouts")
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", color="#e5e5e5", lw=.6)
    ax.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "figS1_importance.png"), dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    figure1()
    figure2()
    figure_s1()
    print(f"Figures written to {FIGDIR}/")
