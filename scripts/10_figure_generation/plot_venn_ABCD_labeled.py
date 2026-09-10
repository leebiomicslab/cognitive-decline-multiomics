import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'config'))
from paths import *
"""
Purpose:          Venn diagram of DA taxa overlap across Models A-D
Manuscript:       Results - Microbiome DA robustness
Figure/Table:     Fig. 5a
"""


"""
plot_venn_ABCD_labeled.py
Draw the 4-way ABCD microbiota Venn diagram with taxa labels for every
non-empty region, laid out as an annotated table below the Venn plot.
"""

import os
import itertools
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
from venn import venn

# ?? Data ?????????????????????????????????????????????????????????????????????
def get_sig_taxa(model_path):
    df = pd.read_csv(os.path.join(model_path, 'microbiota_da_results_final.csv'))
    return set(df[df['Union_Sig'] == True]['taxon'].tolist())

base = str(MICRO_ASSOC_DIR)
sig = {
    'A': get_sig_taxa(os.path.join(base, 'Model_A_Main')),
    'B': get_sig_taxa(os.path.join(base, 'Model_B_Sens1_Systemic')),
    'C': get_sig_taxa(os.path.join(base, 'Model_C_Sens2_MicrobiomeDrugs')),
    'D': get_sig_taxa(os.path.join(base, 'Model_D_Sens3_TxSummary')),
}

# ?? Compute all 15 non-empty intersections ????????????????????????????????????
keys = ['A', 'B', 'C', 'D']
labels_map = {
    'A':    'A only',
    'B':    'B only',
    'C':    'C only',
    'D':    'D only',
    'AB':   'A ??B',
    'AC':   'A ??C',
    'AD':   'A ??D',
    'BC':   'B ??C',
    'BD':   'B ??D',
    'CD':   'C ??D',
    'ABC':  'A ??B ??C',
    'ABD':  'A ??B ??D',
    'ACD':  'A ??C ??D',
    'BCD':  'B ??C ??D',
    'ABCD': 'A ??B ??C ??D',
}

regions = {}
for r in range(1, len(keys) + 1):
    for combo in itertools.combinations(keys, r):
        # taxa IN all combo keys but NOT in any other key
        in_set  = set.intersection(*[sig[k] for k in combo])
        out_set = set.union(*[sig[k] for k in keys if k not in combo]) if len(combo) < len(keys) else set()
        exclusive = in_set - out_set
        key = ''.join(combo)
        regions[key] = sorted(exclusive)

# Keep only non-empty regions, ordered largest?mallest
non_empty = [(k, v) for k, v in regions.items() if v]
non_empty.sort(key=lambda x: -len(x[1]))

# ?? Layout ???????????????????????????????????????????????????????????????????
n_cols = 4  # columns in the annotation table
n_rows_table = -(-len(non_empty) // n_cols)  # ceiling division

fig = plt.figure(figsize=(20, 10 + n_rows_table * 1.8))
gs = GridSpec(2, 1, height_ratios=[5, n_rows_table * 1.2], hspace=0.08)

# ?? Top panel: Venn ???????????????????????????????????????????????????????????
ax_venn = fig.add_subplot(gs[0])
dataset_abcd = {
    'Model A\n(Main)':            sig['A'],
    'Model B\n(Systemic Sens.)':  sig['B'],
    'Model C\n(Microbiome Drugs)':sig['C'],
    'Model D\n(TxSummary Sens.)': sig['D'],
}
venn(dataset_abcd, ax=ax_venn)
ax_venn.set_title('Significant Taxa Overlap (Models A, B, C, D)',
                  fontsize=16, fontweight='bold', pad=14)
leg = ax_venn.get_legend()
if leg:
    leg.set_bbox_to_anchor((1.18, 1.0))

# ?? Bottom panel: annotation table ???????????????????????????????????????????
ax_tbl = fig.add_subplot(gs[1])
ax_tbl.axis('off')

# Divider line between panels
fig.add_artist(plt.Line2D([0.04, 0.96], [gs[1].get_position(fig).y1,
                                          gs[1].get_position(fig).y1],
                           transform=fig.transFigure,
                           color='#CCCCCC', linewidth=1.5))

# Colour palette matching the Venn circles
region_colors = {
    'A':    '#7B68EE',
    'B':    '#6BAED6',
    'C':    '#74C476',
    'D':    '#FDD835',
    'AB':   '#5A7EC7',
    'AC':   '#6BAE91',
    'AD':   '#B5A0DC',
    'BC':   '#68BEB8',
    'BD':   '#A8D08D',
    'CD':   '#AEDB91',
    'ABC':  '#6CADB0',
    'ABD':  '#8989CC',
    'ACD':  '#8EB899',
    'BCD':  '#82C985',
    'ABCD': '#9E9E9E',
}

# Draw each non-empty region as a box with taxa list
cell_w = 1.0 / n_cols
cell_h = 1.0 / (n_rows_table + 0.5)

for idx, (rid, taxa) in enumerate(non_empty):
    col = idx % n_cols
    row = idx // n_cols

    x0 = col * cell_w + 0.01
    y0 = 1.0 - (row + 1) * cell_h + 0.01

    color = region_colors.get(rid, '#DDDDDD')
    title = labels_map.get(rid, rid)
    n_taxa = len(taxa)

    # Header box
    ax_tbl.text(x0 + cell_w * 0.5, y0 + cell_h - 0.02,
                f'{title}  (n={n_taxa})',
                transform=ax_tbl.transAxes,
                fontsize=10, fontweight='bold', ha='center', va='top',
                color='white',
                bbox=dict(boxstyle='round,pad=0.35', facecolor=color,
                          edgecolor='none', alpha=0.9))

    # Taxa list
    taxa_str = '\n'.join(f'??{t}' for t in taxa)
    ax_tbl.text(x0 + 0.01, y0 + cell_h - 0.055,
                taxa_str,
                transform=ax_tbl.transAxes,
                fontsize=8.5, ha='left', va='top',
                color='#333333', family='Arial')

    # Cell border
    rect = mpatches.FancyBboxPatch(
        (x0, y0), cell_w - 0.02, cell_h - 0.02,
        boxstyle='round,pad=0.01',
        linewidth=1, edgecolor=color, facecolor='none',
        transform=ax_tbl.transAxes, clip_on=False
    )
    ax_tbl.add_patch(rect)

# ?? Save ?????????????????????????????????????????????????????????????????????
out_dir = str(FIGURES_DIR / "microbiota")
os.makedirs(out_dir, exist_ok=True)
out_path = os.path.join(out_dir, 'v2_microbiota_sens_venn_ABCD_labeled.png')
plt.savefig(out_path, dpi=600, bbox_inches='tight', facecolor='white')
plt.close()
print(f'Saved labeled ABCD Venn to {out_path}')

