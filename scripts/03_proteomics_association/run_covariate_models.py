"""
Purpose:          Protein-wise covariate-adjusted ANCOVA (OLS with HC3 robust SE)
                  across four sensitivity models (Models A-D). Computes omnibus
                  F-test for group effect, BH-FDR correction, and adjusted
                  marginal means for candidate proteins.
Manuscript:       Methods — Proteomic association analysis
Figure/Table:     Fig. 3a-c; Supplementary Fig. 5; Supplementary Table 2
Input:            PROTEOMICS_LOG2, VASCULAR_CSV, MEDICATION_XLSX,
                  GROUP_CSV (k-means_3group.csv), COMMON_SAMPLES_CSV
Output:           results/proteomics_covariate_analysis/<Model>/*
Main dependencies: pandas, numpy, statsmodels, matplotlib, seaborn, adjustText
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'config'))
from paths import *

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import statsmodels.api as sm
import statsmodels.formula.api as smf
from adjustText import adjust_text
from statsmodels.stats.multitest import multipletests

# High-impact journal style settings (Ultra-clear version)
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 25
plt.rcParams['axes.linewidth'] = 1
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.titleweight'] = 'bold'
plt.rcParams['axes.titlesize'] = 26
plt.rcParams['axes.labelsize'] = 22
plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18
plt.rcParams['legend.fontsize'] = 18
plt.rcParams['legend.title_fontsize'] = 20

# Paths (from config/paths.py)
prot_file = PROTEOMICS_LOG2
vasc_file = VASCULAR_CSV
med_file = MEDICATION_XLSX
group_file = GROUP_CSV
common_file = COMMON_SAMPLES_CSV
out_dir_base = PROT_ASSOC_DIR



def load_data():
    print("Loading datasets...")
    prot = pd.read_csv(prot_file, index_col=0)
    vasc = pd.read_csv(vasc_file)
    med = pd.read_excel(med_file)
    group = pd.read_csv(group_file)
    common = pd.read_csv(common_file)

    for df in [vasc, med, group, common]:
        if "ID" in df.columns:
            df.rename(columns={"ID": "id"}, inplace=True)
        if "id" in df.columns:
            df["id"] = df["id"].astype(str).str.strip()
    prot.index = prot.index.astype(str).str.strip()

    # Rename specific proteins for cleaner figures
    rename_proteins = {
        "Putative uncharacterized protein encoded by LINC00469": "LINC00469",
        "NFX1-type zinc finger-containing protein 1": "ZNFX1",
        "CFAP410&CPXM1&PHB2": "PHB2",
        "KIF2A&KIF2B": "KIF2A",
        "Mannose-binding protein C": "MBP-C",
        "Mitogen-activated protein kinase kinase kinase 21": "MAP3K21",
        "TRAF3&Testis-specific serine/threonine-protein kinase 6": "TRAF3",
        "DLG5&SLC49A3": "DLG5",
        "DSG1&DSG4": "DSG1"
    }
    prot.rename(columns=rename_proteins, inplace=True)

    # Required columns
    target_vasc = [
        "id",
        "age",
        "gender",
        "edu",
        "bmi",
        "smoking",
        "alcohol",
        "HTN",
        "DM",
        "CVD",
        "CRP",
        "egfr",
    ]
    target_med = ["id", "Statin", "dm_med", "htn_med", "chol_med"]

    # Case insensitive extraction just in case
    vasc_cols = {col.lower().strip(): col for col in vasc.columns}
    cv_vasc = [
        vasc_cols[x.lower().strip()]
        for x in target_vasc
        if x.lower().strip() in vasc_cols
    ]

    med_cols = {col.lower().strip(): col for col in med.columns}
    cv_med = [
        med_cols[x.lower().strip()] for x in target_med if x.lower().strip() in med_cols
    ]

    clin = pd.merge(group[["id", "group"]], vasc[cv_vasc], on="id", how="inner")
    clin = pd.merge(clin, med[cv_med], on="id", how="inner")

    # Rename specifically back to our formula targets
    rename_mapping = {vc: target_vasc[i] for i, vc in enumerate(cv_vasc)}
    rename_mapping.update({mc: target_med[i] for i, mc in enumerate(cv_med)})
    clin.rename(columns=rename_mapping, inplace=True)

    clin = clin[clin["id"].isin(common["id"])]

    # Ensure categorical columns are standard int
    cat_cols = [
        "gender",
        "smoking",
        "alcohol",
        "HTN",
        "DM",
        "CVD",
        "Statin",
        "dm_med",
        "htn_med",
        "chol_med",
    ]
    num_cols = ["age", "edu", "bmi", "CRP", "egfr"]

    for c in cat_cols + num_cols:
        clin[c] = pd.to_numeric(clin[c], errors="coerce")

    clin = clin.dropna().set_index("id")

    for c in cat_cols:
        clin[c] = clin[c].astype(int)

    for c in num_cols:
        clin[c] = clin[c].astype(float)

    clin.rename(columns={"group": "Group"}, inplace=True)
    clin["Group"] = (
        clin["Group"].astype(str).astype(object)
    )  # strictly base object string

    valid_ids = clin.index.intersection(prot.index)
    clin_fn = clin.loc[valid_ids]
    prot_fn = prot.loc[valid_ids]

    print(f"Final cohort size for modeling: {len(valid_ids)}")
    return clin_fn, prot_fn


def calculate_marginal_means(model, data):
    groups = ["C1", "C2", "C3"]
    pred_data = []

    for g in groups:
        row = {"Group": g}
        for col in data.columns:
            if col in ["Group", "Expression", "id"]:
                continue
            if data[col].nunique() > 5:
                row[col] = data[col].mean()
            else:
                row[col] = data[col].mode()[0]
        pred_data.append(row)

    pred_df = pd.DataFrame(pred_data)
    for col in data.columns:
        if col in pred_df.columns:
            pred_df[col] = pred_df[col].astype(data[col].dtype)

    try:
        prediction = model.get_prediction(pred_df)
        summary = prediction.summary_frame(alpha=0.05)
        summary["Group"] = groups
        return summary
    except Exception as e:
        print("Prediction error:", e)
        return None


def run_ancova(clinical_df, prot_df, formula, model_name):
    print(f"\n[{model_name}] Running ANCOVA over {prot_df.shape[1]} proteins...")
    out_dir = os.path.join(out_dir_base, model_name)
    os.makedirs(out_dir, exist_ok=True)

    results = []
    proteins = prot_df.columns.tolist()

    for i, protein in enumerate(proteins):
        if i % 100 == 0:
            print(f"Processing {i}/{len(proteins)}...", end="\r")

        y = prot_df[protein]
        data = clinical_df.copy()
        data["Expression"] = y

        # Omit NaNs in expression
        data = data.dropna(subset=["Expression"])
        data = data[~np.isinf(data["Expression"])]

        if len(data) < 10:
            continue

        try:
            model = smf.ols(formula, data=data).fit(cov_type="HC3")
            hypotheses = 'C(Group, Treatment(reference="C1"))[T.C2] = 0, C(Group, Treatment(reference="C1"))[T.C3] = 0'
            f_test = model.f_test(hypotheses)
            group_p = f_test.pvalue

            term_c2 = 'C(Group, Treatment(reference="C1"))[T.C2]'
            term_c3 = 'C(Group, Treatment(reference="C1"))[T.C3]'
            
            p_c2 = model.pvalues.get(term_c2, np.nan)
            p_c3 = model.pvalues.get(term_c3, np.nan)

            beta_c2 = model.params[term_c2]
            ci_c2 = model.conf_int().loc[term_c2]
            beta_c3 = model.params[term_c3]
            ci_c3 = model.conf_int().loc[term_c3]

            marg_means = calculate_marginal_means(model, data)
            if marg_means is not None:
                mm_c1 = marg_means.loc[marg_means["Group"] == "C1", "mean"].values[0]
                mm_c2 = marg_means.loc[marg_means["Group"] == "C2", "mean"].values[0]
                mm_c3 = marg_means.loc[marg_means["Group"] == "C3", "mean"].values[0]
                se_c1 = marg_means.loc[marg_means["Group"] == "C1", "mean_se"].values[0]
                se_c2 = marg_means.loc[marg_means["Group"] == "C2", "mean_se"].values[0]
                se_c3 = marg_means.loc[marg_means["Group"] == "C3", "mean_se"].values[0]
            else:
                mm_c1 = mm_c2 = mm_c3 = se_c1 = se_c2 = se_c3 = np.nan

            results.append(
                {
                    "Protein": protein,
                    "P_Value": float(group_p),
                    "P_C2": p_c2,
                    "P_C3": p_c3,
                    "Beta_C2": beta_c2,
                    "Beta_C3": beta_c3,
                    "CI_C2_Low": ci_c2[0],
                    "CI_C2_High": ci_c2[1],
                    "CI_C3_Low": ci_c3[0],
                    "CI_C3_High": ci_c3[1],
                    "AdjMean_C1": mm_c1,
                    "AdjMean_C2": mm_c2,
                    "AdjMean_C3": mm_c3,
                    "SE_C1": se_c1,
                    "SE_C2": se_c2,
                    "SE_C3": se_c3,
                    "N_Samples": len(data),
                }
            )
        except Exception as e:
            print(f"Error on {protein}: {e}")
            continue

    res_df = pd.DataFrame(results)
    if res_df.empty:
        print("No models converged successfully!")
        return None, None

    # FDR
    reject, q_values, _, _ = multipletests(
        res_df["P_Value"], alpha=0.05, method="fdr_bh"
    )
    res_df["Q_Value"] = q_values
    res_df["Nominal_Sig"] = res_df["P_Value"] < 0.05
    res_df["FDR_Sig"] = res_df["Q_Value"] < 0.05
    
    _, q_c2, _, _ = multipletests(res_df["P_C2"].fillna(1), alpha=0.05, method="fdr_bh")
    _, q_c3, _, _ = multipletests(res_df["P_C3"].fillna(1), alpha=0.05, method="fdr_bh")
    res_df["Q_C2"] = q_c2
    res_df["Q_C3"] = q_c3
    res_df.sort_values("P_Value", inplace=True)

    res_df.to_csv(
        os.path.join(out_dir, "proteomics_ancova_detailed_results.csv"), index=False
    )

    nom_sig = res_df[res_df["P_Value"] < 0.05]
    nom_sig.to_csv(
        os.path.join(out_dir, "nominal_significant_proteins.csv"), index=False
    )

    return res_df, out_dir


def plot_fig1_ranked_pvals(df, out_dir):
    df = df.sort_values("P_Value", ascending=True).reset_index(drop=True)
    df['logP'] = -np.log10(df['P_Value'])
    df['rank'] = df.index + 1

    plt.figure(figsize=(14, 9))
    sns.set_style("white")

    sig_mask = df['P_Value'] < 0.05
    color_sig = "#D73027"
    color_ns = "#BDBDBD"

    plt.scatter(df.loc[~sig_mask, 'rank'], df.loc[~sig_mask, 'logP'],
                alpha=0.3, s=15, color=color_ns, label='NS', rasterized=True)

    plt.scatter(df.loc[sig_mask, 'rank'], df.loc[sig_mask, 'logP'],
                alpha=0.9, s=50, color=color_sig, label='P < 0.05',
                edgecolors='white', linewidths=0.8, zorder=3)

    sig_df = df[sig_mask].copy()
    n_sig = len(sig_df)

    label_x_pos = df['rank'].max() * 1.15
    if n_sig > 0:
        y_logp_max = sig_df['logP'].max()
        y_logp_min = 0.15
        y_positions = np.linspace(y_logp_max, y_logp_min, n_sig)

        for idx, (i, row) in enumerate(sig_df.iterrows()):
            p_name = str(row['Protein']).split(';')[0]
            target_y = y_positions[idx]

            plt.annotate(p_name,
                         xy=(row['rank'], row['logP']),
                         xytext=(label_x_pos, target_y),
                         fontsize=14, fontweight='bold', family='Arial',
                         va='center', ha='left',
                         arrowprops=dict(arrowstyle='->', color='#555555', lw=1.2, alpha=0.35,
                                         connectionstyle="arc3,rad=0.05"))

    plt.xlim(-50, df['rank'].max() * 1.7)
    plt.axhline(-np.log10(0.05), color='#4575B4', linestyle='--', linewidth=1.2, label='P=0.05', alpha=0.8)

    handles, labels = plt.gca().get_legend_handles_labels()
    plt.figlegend(handles, labels, loc='lower right',
                  bbox_to_anchor=(0.99, 0.02), fontsize=16,
                  frameon=True, framealpha=0.9, edgecolor='#D0D0D0')

    plt.xlabel('Protein Rank (by P-value)', labelpad=12, fontweight='bold')
    plt.ylabel(r'-$\log_{10}$(P-value)', labelpad=12, fontweight='bold')
    plt.title('Ranked Significance of Proteomic Differences', pad=25, fontsize=22, fontweight='bold')

    sns.despine(trim=True, offset=10)
    plt.grid(axis='y', linestyle='--', alpha=0.2)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "Figure1_Ranked_Pvals.png"), dpi=600, bbox_inches='tight')
    plt.close()

def plot_fig2_forest(df, out_dir, ref_order=None):
    sig_mask = df["P_Value"] < 0.05
    df_sig = df[sig_mask].copy()
    if len(df_sig) == 0:
        return

    if ref_order is not None:
        def get_sort_key(row):
            prot = row["Protein"]
            if prot in ref_order:
                return ref_order.index(prot)
            return len(ref_order) + row["P_Value"]
        df_sig["sort_key"] = df_sig.apply(get_sort_key, axis=1)
        df_sig = df_sig.sort_values("sort_key", ascending=True).drop(columns=["sort_key"])
    else:
        df_sig = df_sig.sort_values("P_Value", ascending=True)

    df_sig = df_sig.iloc[::-1]
    proteins = df_sig["Protein"].tolist()
    y_pos = np.arange(len(proteins))
    figsize_h = max(10, len(df_sig) * 0.55)
    plt.figure(figsize=(14, figsize_h))
    
    for i, row in enumerate(df_sig.itertuples()):
        ci_lo = row.CI_C2_Low
        ci_hi = row.CI_C2_High
        beta = row.Beta_C2
        if ci_lo > 0 or ci_hi < 0:
            marker, fillstyle = "o", "full"
        else:
            marker, fillstyle = "o", "none"
        plt.errorbar(
            beta, i + 0.15, xerr=[[beta - ci_lo], [ci_hi - beta]], fmt=marker, color="#eca362",
            capsize=3, alpha=0.9, markersize=8, markerfacecolor="#eca362" if fillstyle == "full" else "white", markeredgecolor="#eca362"
        )

    for i, row in enumerate(df_sig.itertuples()):
        ci_lo = row.CI_C3_Low
        ci_hi = row.CI_C3_High
        beta = row.Beta_C3
        if ci_lo > 0 or ci_hi < 0:
            marker, fillstyle = "o", "full"
        else:
            marker, fillstyle = "o", "none"
        plt.errorbar(
            beta, i - 0.15, xerr=[[beta - ci_lo], [ci_hi - beta]], fmt=marker, color="#6fcbbc",
            capsize=3, alpha=0.9, markersize=8, markerfacecolor="#6fcbbc" if fillstyle == "full" else "white", markeredgecolor="#6fcbbc"
        )
        
    plt.axvline(0, color="grey", linestyle="--", linewidth=1)
    
    # Increase font size for feature names
    plt.yticks(y_pos, proteins, fontsize=24)
    plt.xlabel("Adjusted log2 Fold Change (Beta)", fontsize=18, fontweight='bold')
    plt.title("Effect sizes of candidate proteins (F-test P < 0.05)", fontsize=20, fontweight='bold', pad=15)

    plt.grid(axis="x", linestyle=":", alpha=0.6)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "Figure2_Forest_Plot.png"), dpi=600)
    plt.close()

    # Create standalone legend
    fig_leg = plt.figure(figsize=(6, 3))
    ax_leg = fig_leg.add_subplot(111)
    ax_leg.axis('off')
    
    handles = [
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="#eca362", markeredgecolor="#eca362", markersize=10),
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="#6fcbbc", markeredgecolor="#6fcbbc", markersize=10),
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="gray", markeredgecolor="gray", markersize=10),
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="white", markeredgecolor="gray", markersize=10)
    ]
    labels = ["C2 vs C1", "C3 vs C1", "95% CI does not cross 0", "95% CI crosses 0"]
    
    ax_leg.legend(handles, labels, loc='center', fontsize=18)
    fig_leg.savefig(os.path.join(out_dir, "Figure2_Forest_Plot_Legend.png"), dpi=600, bbox_inches='tight')
    plt.close(fig_leg)

def plot_fig3_marginal_means(df, out_dir, ref_order=None):
    sig_mask = df['P_Value'] < 0.05
    df_sig = df[sig_mask].copy()
    n_sig = len(df_sig)
    if n_sig == 0: return

    if ref_order is not None:
        def get_sort_key(row):
            prot = row["Protein"]
            if prot in ref_order:
                return ref_order.index(prot)
            return len(ref_order) + row["P_Value"]
        df_sig["sort_key"] = df_sig.apply(get_sort_key, axis=1)
        df_sig = df_sig.sort_values("sort_key", ascending=True).drop(columns=["sort_key"])
    else:
        df_sig = df_sig.sort_values("P_Value", ascending=True)
    cols = 5
    rows = int(np.ceil(n_sig / cols))
    figsize_h = max(6, rows * 3.5)
    figs, axes = plt.subplots(rows, cols, figsize=(20, figsize_h))
    if n_sig == 1: axes = [axes]
    else: axes = axes.flatten()
    groups = ['C1', 'C2', 'C3']
    colors = ['#ee7a5b', '#eca362', '#6fcbbc']
    for i in range(len(axes)):
        ax = axes[i]
        if i < n_sig:
            row = df_sig.iloc[i]
            means = [float(row['AdjMean_C1']), float(row['AdjMean_C2']), float(row['AdjMean_C3'])]
            es = [float(row['SE_C1']), float(row['SE_C2']), float(row['SE_C3'])]

            ax.plot([0, 1, 2], means, linestyle='--', color='black', alpha=0.6, zorder=0, linewidth=2)
            ax.errorbar([0, 1, 2], means, yerr=es, fmt='none', capsize=5, color='black', zorder=1)
            for g_idx, g in enumerate(groups):
                ax.plot(g_idx, means[g_idx], 'o', markersize=10, color=colors[g_idx], alpha=0.9, zorder=2)

            ax.set_xticks([0, 1, 2])
            ax.set_xticklabels(groups)

            p_name = str(row['Protein'])
            p_name = p_name[:27] + "..." if len(p_name) > 30 else p_name
            ax.set_title(f"{p_name}\nP={row['P_Value']:.4f}", fontsize=18)
            if i % cols == 0: ax.set_ylabel('Adj. Expr')
        else:
            ax.axis('off')
    plt.suptitle('Adjusted Marginal Means (candidate proteins)', fontsize=22, y=0.995, fontweight='bold')
    plt.tight_layout(rect=[0, 0.02, 1, 0.98])
    plt.savefig(os.path.join(out_dir, "Figure3_Marginal_Means.png"), dpi=600)
    plt.close()


def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)
    clin_df, prot_df = load_data()

    formulas = {
        "Main_Model": 'Expression ~ C(Group, Treatment(reference="C1")) + age + C(gender) + edu + bmi + C(smoking) + C(alcohol) + C(HTN) + C(DM) + C(CVD)',
        "Sensitivity_1_Systemic": 'Expression ~ C(Group, Treatment(reference="C1")) + age + C(gender) + edu + bmi + C(smoking) + C(alcohol) + C(HTN) + C(DM) + C(CVD) + CRP + egfr',
        "Sensitivity_2_Statin": 'Expression ~ C(Group, Treatment(reference="C1")) + age + C(gender) + edu + bmi + C(smoking) + C(alcohol) + C(HTN) + C(DM) + C(CVD) + C(Statin)',
        "Sensitivity_3_Medications": 'Expression ~ C(Group, Treatment(reference="C1")) + age + C(gender) + edu + bmi + C(smoking) + C(alcohol) + C(HTN) + C(DM) + C(CVD) + C(dm_med) + C(htn_med) + C(chol_med)',
    }

    ref_order = None
    for m_name, form in formulas.items():
        res_df, m_dir = run_ancova(clin_df, prot_df, form, m_name)
        if m_name == "Main_Model" and res_df is not None:
            ref_order = res_df[res_df["P_Value"] < 0.05].sort_values("P_Value", ascending=True)["Protein"].tolist()
            
        plot_fig1_ranked_pvals(res_df, m_dir)
        plot_fig2_forest(res_df, m_dir, ref_order)
        plot_fig3_marginal_means(res_df, m_dir, ref_order)
        print(f"[{m_name}] Completed and visually constructed!")
        
    print("Generating cross-model sensitivity summary...")
    main_df = pd.read_csv(os.path.join(out_dir_base, "Main_Model", "proteomics_ancova_detailed_results.csv"))
    main_cands = main_df[main_df["Nominal_Sig"] == True].copy()
    
    sens1_df = pd.read_csv(os.path.join(out_dir_base, "Sensitivity_1_Systemic", "proteomics_ancova_detailed_results.csv"))
    sens2_df = pd.read_csv(os.path.join(out_dir_base, "Sensitivity_2_Statin", "proteomics_ancova_detailed_results.csv"))
    sens3_df = pd.read_csv(os.path.join(out_dir_base, "Sensitivity_3_Medications", "proteomics_ancova_detailed_results.csv"))
    
    summary_data = []
    for _, row in main_cands.iterrows():
        prot = row["Protein"]
        main_p = row["P_Value"]
        main_b2 = row["Beta_C2"]
        main_b3 = row["Beta_C3"]
        
        s1 = sens1_df[sens1_df["Protein"] == prot].iloc[0] if prot in sens1_df["Protein"].values else None
        s2 = sens2_df[sens2_df["Protein"] == prot].iloc[0] if prot in sens2_df["Protein"].values else None
        s3 = sens3_df[sens3_df["Protein"] == prot].iloc[0] if prot in sens3_df["Protein"].values else None
        
        s1_p = s1["P_Value"] if s1 is not None else np.nan
        s2_p = s2["P_Value"] if s2 is not None else np.nan
        s3_p = s3["P_Value"] if s3 is not None else np.nan
        
        s1_b2 = s1["Beta_C2"] if s1 is not None else np.nan
        s1_b3 = s1["Beta_C3"] if s1 is not None else np.nan
        s2_b2 = s2["Beta_C2"] if s2 is not None else np.nan
        s2_b3 = s2["Beta_C3"] if s2 is not None else np.nan
        s3_b2 = s3["Beta_C2"] if s3 is not None else np.nan
        s3_b3 = s3["Beta_C3"] if s3 is not None else np.nan
        
        # Check direction stability
        dirs_b2 = [np.sign(b) for b in [main_b2, s1_b2, s2_b2, s3_b2] if not np.isnan(b)]
        dirs_b3 = [np.sign(b) for b in [main_b3, s1_b3, s2_b3, s3_b3] if not np.isnan(b)]
        
        dir_stable = (len(set(dirs_b2)) <= 1) and (len(set(dirs_b3)) <= 1)
        
        # Models retained
        retained = []
        if main_p < 0.05: retained.append("Main")
        if s1_p < 0.05: retained.append("Systemic")
        if s2_p < 0.05: retained.append("Statin")
        if s3_p < 0.05: retained.append("Medications")
        
        summary_data.append({
            "Protein": prot,
            "Main P": main_p,
            "Main Beta_C2": main_b2,
            "Main Beta_C3": main_b3,
            "Systemic P": s1_p,
            "Statin P": s2_p,
            "Medication P": s3_p,
            "Direction stable": dir_stable,
            "Models retained": ", ".join(retained)
        })
        
    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv(os.path.join(out_dir_base, "Cross_Model_Sensitivity_Summary.csv"), index=False)
    print("Cross-model sensitivity summary generated!")


if __name__ == "__main__":
    main()
