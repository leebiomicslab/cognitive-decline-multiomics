import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'config'))
from paths import *
"""
Purpose:          Path diagram for sequential mediation structures
Manuscript:       Methods - Mediation analysis
Figure/Table:     Fig. 8a
"""


"""
Publication-quality serial mediation path diagram for 4 alternative models.
Aesthetic: white background, muted palette, thin lines, journal-ready.
"""
import os
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

matplotlib.rcParams.update({
    "pdf.fonttype": 42,
    "ps.fonttype":  42,
    "font.family":  "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "figure.facecolor": "white",
    "axes.facecolor":   "white",
})

# ?? Palette ???????????????????????????????????????????????????????????????????
C_MIC  = "#c0745a"   # warm terracotta ??Microbiome
C_MET  = "#c59a45"   # muted amber    ??帠-Butyrobetaine
C_PROT = "#4e8f85"   # slate teal     ??Protein
C_MOCA = "#5a6e82"   # steel blue     ??MoCA
C_TXT  = "#1a1a1a"
C_EDGE = "#444444"
C_DIRE = "#999999"   # direct path (c') ??visible but subordinate

# ?? Helper: draw one panel ????????????????????????????????????????????????????
def draw_panel(ax, panel_label, nodes, edges, ind, lo, hi, p):
    ax.set_xlim(0, 10)
    ax.set_ylim(-1.5, 6.5)
    ax.axis("off")

    ax.text(-0.02, 1.02, panel_label, transform=ax.transAxes,
            fontsize=13, fontweight="bold", color=C_TXT, va="top")

    box_w, box_h = 2.4, 0.85
    box_patches = {}
    for name, (cx, cy) in nodes.items():
        color = {
            "Microbiome\nCISS":   C_MIC,
            "帠-Butyrobetaine\n(GBB)": C_MET,
            "Protein\nCISS":      C_PROT,
            "MoCA\nScore":        C_MOCA,
        }.get(name, "#cccccc")

        rect = FancyBboxPatch(
            (cx - box_w/2, cy - box_h/2), box_w, box_h,
            boxstyle="round,pad=0.15",
            facecolor=color, edgecolor="white",
            linewidth=0, zorder=3, alpha=0.88
        )
        ax.add_patch(rect)
        box_patches[name] = (cx, cy)

        ax.text(cx, cy, name, ha="center", va="center",
                fontsize=11, fontweight="bold", color="white",
                zorder=4, multialignment="center",
                linespacing=1.25)

    fig_w_in, fig_h_in = ax.get_figure().get_size_inches()
    x_range = ax.get_xlim()[1] - ax.get_xlim()[0]
    y_range = ax.get_ylim()[1] - ax.get_ylim()[0]
    ax_w_in = ax.get_position().width  * fig_w_in
    ax_h_in = ax.get_position().height * fig_h_in
    pts_per_xu = ax_w_in / x_range * 72   
    pts_per_yu = ax_h_in / y_range * 72   

    for (src, dst, label, is_direct) in edges:
        x1, y1 = nodes[src]
        x2, y2 = nodes[dst]
        dx, dy = x2 - x1, y2 - y1
        L_xu = (dx**2 + dy**2) ** 0.5   

        ux, uy = dx / (L_xu + 1e-9), dy / (L_xu + 1e-9)
        shrink_x = abs(ux) * (box_w  / 2) * pts_per_xu
        shrink_y = abs(uy) * (box_h  / 2) * pts_per_yu
        shrink = (shrink_x**2 + shrink_y**2) ** 0.5

        lc  = C_DIRE if is_direct else C_EDGE
        lw  = 1.2    if is_direct else 1.5
        ls  = (0, (5, 4)) if is_direct else "solid"

        ax.annotate(
            "", xy=(x2, y2), xytext=(x1, y1),
            arrowprops=dict(
                arrowstyle="-|>",
                mutation_scale=12,
                lw=lw, color=lc, linestyle=ls,
                shrinkA=shrink + 12,
                shrinkB=shrink + 12,
                connectionstyle="arc3,rad=0.0"
            ),
            zorder=2
        )

        mx = (x1 + x2) / 2
        my = (y1 + y2) / 2
        if is_direct:
            my -= 0.55
        elif abs(dx) < 0.5:
            mx -= 0.85
        elif abs(dy) < 0.5:
            my += 0.38
        else:
            mx += 0.85

        ax.text(mx, my, label, ha="center", va="center",
                fontsize=10, color=lc if is_direct else C_TXT,
                fontweight="bold",
                bbox=dict(facecolor="white", edgecolor="none",
                          alpha=0.9, pad=1.5),
                zorder=5)

    summary = (
        f"Serial indirect effect = {ind:.4f}\n"
        f"95% CI: {lo:.4f} to {hi:.4f}\n"
        f"p = {p:.4f}" + ("  *" if p < 0.05 else " (ns)")
    )
    ax.text(5.0, 3.0, summary,
            ha="center", va="center", fontsize=10.5, color=C_TXT,
            fontweight="bold",
            linespacing=1.6,
            bbox=dict(facecolor="#fdfdfd", edgecolor="#999999",
                      boxstyle="round,pad=0.8", linewidth=1.0, alpha=0.95),
            zorder=6)

# ?? Data & Layouts ????????????????????????????????????????????????????????????
MET = "帠-Butyrobetaine\n(GBB)"
M_CISS = "Microbiome\nCISS"
P_CISS = "Protein\nCISS"
MOCA = "MoCA\nScore"

fig, axs = plt.subplots(4, 1, figsize=(8.5, 18), facecolor="white")

# Model 1
n1 = {M_CISS: (2.0, 1.2), MET: (2.0, 4.8), P_CISS: (8.0, 4.8), MOCA: (8.0, 1.2)}
e1 = [(M_CISS, MET, "帣 = -0.182", False), (MET, P_CISS, "帣 = -0.105", False), (P_CISS, MOCA, "帣 = -0.172", False), (M_CISS, MOCA, "帣 = -0.201 (direct)", True)]
draw_panel(axs[0], "A  Primary model: Microbiome ??帠-Butyrobetaine ??Protein ??MoCA", n1, e1, -0.0033, -0.0080, -0.0004, 0.0196)

# Model 2
n2 = {P_CISS: (2.0, 1.2), MET: (2.0, 4.8), M_CISS: (8.0, 4.8), MOCA: (8.0, 1.2)}
e2 = [(P_CISS, MET, "帣 = -0.113", False), (MET, M_CISS, "帣 = -0.178", False), (M_CISS, MOCA, "帣 = -0.201", False), (P_CISS, MOCA, "帣 = -0.172 (direct)", True)]
draw_panel(axs[1], "B  Alternative 1: Protein ??帠-Butyrobetaine ??Microbiome ??MoCA", n2, e2, -0.0040, -0.0091, -0.0007, 0.0152)

# Model 3
n3 = {M_CISS: (2.0, 1.2), P_CISS: (2.0, 4.8), MET: (8.0, 4.8), MOCA: (8.0, 1.2)}
e3 = [(M_CISS, P_CISS, "帣 = 0.055", False), (P_CISS, MET, "帣 = -0.103", False), (MET, MOCA, "帣 = 0.014", False), (M_CISS, MOCA, "帣 = -0.201 (direct)", True)]
draw_panel(axs[2], "C  Alternative 2: Microbiome ??Protein ??帠-Butyrobetaine ??MoCA", n3, e3, -0.0001, -0.0011, 0.0007, 0.8420)

# Model 4
n4 = {MOCA: (2.0, 1.2), P_CISS: (2.0, 4.8), MET: (8.0, 4.8), M_CISS: (8.0, 1.2)}
e4 = [(MOCA, P_CISS, "帣 = -0.219", False), (P_CISS, MET, "帣 = -0.102", False), (MET, M_CISS, "帣 = -0.166", False), (MOCA, M_CISS, "帣 = -0.243 (direct)", True)]
draw_panel(axs[3], "D  Reverse Sensitivity: MoCA ??Protein ??帠-Butyrobetaine ??Microbiome", n4, e4, -0.0037, -0.0092, -0.0003, 0.0296)

note = (
    "Standardised path coefficients (帣) shown. "
    "Serial indirect effects estimated by nonparametric bootstrapping (5 000 iterations). "
    "CI = bias-corrected 95% confidence interval. "
    "Adjusted for: age, sex, education. "
)
fig.text(0.5, 0.005, note, ha="center", va="bottom", fontsize=8,
         color="#555555", wrap=True, style="italic", transform=fig.transFigure)

plt.subplots_adjust(hspace=0.2, top=0.98, bottom=0.03, left=0.02, right=0.98)

out_dir  = str(FIGURES_DIR / "metabolic")
os.makedirs(out_dir, exist_ok=True)
out_path = os.path.join(out_dir, "v2_sequential_mediation_TMAO_4models.png")
fig.savefig(out_path, dpi=600, bbox_inches="tight", facecolor="white")
print(f"Saved: {out_path}")
plt.close()

