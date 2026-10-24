"""Graphical abstract: care indicators across five household per-capita income bands (PNS 2019).
Reads the CSVs produced by R/03_income_bands.R (results/income_bands_fine).
Usage (from the repository root): python python/graphical_abstract.py [income_bands_dir] [figures_dir]"""
import sys
import logging
import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, NullLocator

import os
D = sys.argv[1] if len(sys.argv) > 1 else os.path.join("results", "income_bands_fine")
FIGDIR = sys.argv[2] if len(sys.argv) > 2 else "figures"
os.makedirs(FIGDIR, exist_ok=True)
a0 = pd.read_csv(f"{D}/A0_contagens_classes_e_faixas.csv"); a0 = a0[a0.por == "income_band"]
b1 = pd.read_csv(f"{D}/B1_depressao_prevalencia_por_faixa.csv"); b1 = b1[b1.por == "income_band"]
a1 = pd.read_csv(f"{D}/A1_distribuicao_faixas_tabela1.csv")
share = a1[(a1.grupo == "coorte") & (a1.variavel == "income_band")].pct.round(0).astype(int).tolist()
BANDS = ["ate_meio_SM", "meio_a_1_SM", "1_a_2_SM", "2_a_5_SM", "mais_5_SM"]

def series(df, name):
    d = df[df.desfecho == name].set_index("cat").loc[BANDS]
    return (d.prevalencia * 100).values, d.eventos.values

rows = [  # label, sublabel, source, name, group
    ("Diagnosed with depression", "PHQ-9 ≥10 screen-positive", b1, "diagnosticado_depressao_entre_PHQ10", "A"),
    ("Current psychotherapy", "any diagnosed disorder", a0, "psychotherapy", "A"),
    ("Reached a mental health specialist", "diagnosed depression, recent care", b1, "chegou_ao_especialista", "A"),
    ("High need gap", "any diagnosed disorder", a0, "high_need_gap", "B"),
    ("Structural barrier to follow-up", "depression, not in regular follow-up", b1, "barreira_estrutural_motivo", "B"),
    ("Feels remitted · thoughts of death · untreated", "reports remission as the reason", b1, "remissao_percebida_ideacao_sem_tratamento", "B"),
    ("Feels remitted · PHQ-9 ≥15 · untreated", "reports remission as the reason", b1, "remissao_percebida_PHQ15_sem_tratamento", "B"),
]

logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
plt.rcParams.update({"font.family": ["Carlito", "Calibri", "DejaVu Sans"], "font.size": 10})
BLUE, ORANGE = "#2a78d6", "#eb6834"
INK, INK2, MUTED, GRID, SURF = "#0b0b0b", "#52514e", "#8a8984", "#e6e5e1", "#fcfcfb"

fig = plt.figure(figsize=(10.5, 6.9), dpi=300)
ax = fig.add_axes([0.075, 0.165, 0.47, 0.64])
x = np.arange(5)
YMIN, YMAX = 11, 200
ax.set_yscale("log"); ax.set_ylim(YMIN, YMAX); ax.set_xlim(-0.12, 4.12)
ax.axhspan(100, YMAX, color=BLUE, alpha=0.05, lw=0)
ax.axhspan(YMIN, 100, color=ORANGE, alpha=0.05, lw=0)
ax.axhline(100, color=INK2, lw=1.0, zorder=1)
for y in [12.5, 25, 50, 150]:
    ax.axhline(y, color=GRID, lw=0.7, zorder=0)

ends = []
for lab, sub, df, name, g in rows:
    v, ev = series(df, name)
    idx = v / v[0] * 100
    c = BLUE if g == "A" else ORANGE
    ax.plot(x, idx, color=c, lw=2.3, solid_capstyle="round", zorder=3)
    for xi, yi, e in zip(x, idx, ev):
        small = e < 20
        ax.scatter([xi], [yi], s=44, zorder=4, color="white" if small else c,
                   edgecolor=c if small else "white", linewidth=1.6)
    ends.append([idx[-1], lab, sub, v, ev, g])

# de-collide right-hand labels in log space
MIN_GAP = {"A": 0.085, "B": 0.125}  # in log10 units
for grp in ("A", "B"):
    items = sorted([e for e in ends if e[5] == grp], key=lambda e: -e[0])
    pos = [np.log10(e[0]) for e in items]
    for i in range(1, len(pos)):
        if pos[i - 1] - pos[i] < MIN_GAP[grp]:
            pos[i] = pos[i - 1] - MIN_GAP[grp]
    if grp == "A":  # keep block below the top and above the 100 line
        shift = max(0, pos[0] - np.log10(170)); pos = [p - shift for p in pos]
        low = np.log10(112)
        if pos[-1] < low: pos = [p + (low - pos[-1]) for p in pos]
    for e, p in zip(items, pos):
        e.append(10 ** p)

tx = 4.22
for yend, lab, sub, v, ev, g, ty in ends:
    c = BLUE if g == "A" else ORANGE
    ax.annotate("", xy=(4.08, yend), xytext=(tx - 0.03, ty), annotation_clip=False,
                arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.6, shrinkA=0, shrinkB=0))
    ratio = v[-1] / v[0]
    dag = "†" if ev[-1] < 20 else ""
    ax.text(tx, ty, f"×{ratio:.2f}{dag}", color=c, fontsize=13, fontweight="bold", va="center", ha="left", clip_on=False)
    ax.text(tx + 0.66, ty * 1.02, lab, color=INK, fontsize=9.6, fontweight="bold", va="bottom", ha="left", clip_on=False)
    ax.text(tx + 0.66, ty * 0.985, f"{v[0]:.1f}% → {v[-1]:.1f}%  ·  {sub}", color=INK2, fontsize=8.4,
            va="top", ha="left", clip_on=False)

# group headers sit just above the first label of their group
for grp, title, col in (("A", "ACCESS TO CARE RESOURCES", BLUE), ("B", "UNMET NEED", ORANGE)):
    top = max(e[6] for e in ends if e[5] == grp)
    ax.text(tx, top * 10 ** 0.075, title, color=col, fontsize=9, fontweight="bold", va="center", clip_on=False)

ax.yaxis.set_major_locator(FixedLocator([12.5, 25, 50, 100, 150]))
ax.yaxis.set_minor_locator(NullLocator())
ax.set_yticklabels(["12.5", "25", "50", "100", "150"], color=INK2)
ax.set_ylabel("Prevalence relative to the poorest band (poorest = 100, log scale)", color=INK2)
ax.set_xticks(x)
ax.set_xticklabels([f"≤½\n({share[0]}%)", f">½–1\n({share[1]}%)", f">1–2\n({share[2]}%)",
                    f">2–5\n({share[3]}%)", f">5\n({share[4]}%)"], color=INK2)
ax.set_xlabel("Household per-capita income, minimum wages (share of the diagnosed cohort)", color=INK2, labelpad=8)
ax.tick_params(length=0, pad=6)
for s in ax.spines.values(): s.set_visible(False)
ax.text(0.03, 176, "more frequent among the better-off", color=BLUE, fontsize=9, style="italic", va="center")
ax.text(0.03, 11.9, "more frequent among the poorest", color=ORANGE, fontsize=9, style="italic", va="center")

fig.text(0.075, 0.955, "Up the income ladder, access to care rises — and unmet need falls by nearly two-thirds or more",
         fontsize=16, fontweight="bold", color=INK, ha="left")
fig.text(0.075, 0.905,
         "Adults with a diagnosed mental disorder or PHQ-9 ≥10, Brazil, PNS 2019. From the poorest band (≤½ minimum wage per person) to the richest (>5),\n"
         "access indicators rise 1.2–1.7-fold, while unmet-need indicators fall to between about one third and one sixth of their level among the poorest.",
         fontsize=10, color=INK2, ha="left", va="top", linespacing=1.35)
fig.text(0.075, 0.018,
         "Weighted prevalence by band; ratios are crude (richest ÷ poorest band). Adjusted relative indices of inequality: Table 6. "
         "Hollow markers and † = fewer than 20 events in the band (imprecise).\nUntreated = no current psychotherapy, medication or integrative practice for depression. "
         "Feels remitted = gives \"no longer depressed\" as the main reason for not attending regularly.",
         fontsize=7.6, color=MUTED, ha="left", va="bottom", linespacing=1.3)
fig.savefig(os.path.join(FIGDIR, "fig_graphical_abstract.png"), dpi=300, facecolor=SURF)
fig.savefig(os.path.join(FIGDIR, "fig_graphical_abstract.pdf"), facecolor=SURF)
for e in ends: print(e[1], [round(z, 1) for z in e[3]], "ratio", round(e[3][-1] / e[3][0], 2), "events", list(e[4]))
