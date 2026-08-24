import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Set overall publication style
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['axes.titlesize'] = 14
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['xtick.labelsize'] = 10
plt.rcParams['ytick.labelsize'] = 10
plt.rcParams['legend.fontsize'] = 10

figures_dir = r"c:\Users\enosh\Desktop\BDSA PROJECT\figures"
tables_dir = r"c:\Users\enosh\Desktop\BDSA PROJECT\results\tables"
os.makedirs(figures_dir, exist_ok=True)

df_h2 = pd.read_csv(os.path.join(tables_dir, "snp_heritability_summary.csv"))
df_bs = pd.read_csv(os.path.join(tables_dir, "between_sex_genetic_correlation.csv"))
df_pw = pd.read_csv(os.path.join(tables_dir, "within_sex_pairwise_genetic_correlation.csv"))
df_cross = pd.read_csv(os.path.join(tables_dir, "cross_sex_cross_trait_correlation.csv"))

# Color palettes
female_color = '#E64B35'  # Red/Coral
male_color = '#4DBBD5'    # Cyan/Blue

# ==========================================
# FIGURE 1: HERITABILITY COMPARISON BY SEX
# ==========================================
fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)

for i, trans in enumerate(["IRNT (quantile)", "Raw (untransformed)"]):
    ax = axes[i]
    sub = df_h2[df_h2['Transformation'] == trans].copy()
    
    # Sort traits consistently
    traits_order = ["Testosterone", "HDL Cholesterol", "Triglycerides", "Glycated Haemoglobin (HbA1c)", "Diastolic Blood Pressure", "Systolic Blood Pressure", "LDL Cholesterol"]
    sub['Trait_Cat'] = pd.Categorical(sub['Trait'], categories=traits_order, ordered=True)
    sub = sub.sort_values('Trait_Cat')
    
    y_pos = np.arange(len(traits_order))
    width = 0.35
    
    sub_f = sub[sub['Sex'] == 'Female'].set_index('Trait').reindex(traits_order)
    sub_m = sub[sub['Sex'] == 'Male'].set_index('Trait').reindex(traits_order)
    
    # Plot bars with error bars
    ax.barh(y_pos - width/2, sub_f['h2'], width, xerr=sub_f['h2_SE']*1.96,
            label='Female', color=female_color, alpha=0.85, capsize=4, edgecolor='black', linewidth=0.8)
    ax.barh(y_pos + width/2, sub_m['h2'], width, xerr=sub_m['h2_SE']*1.96,
            label='Male', color=male_color, alpha=0.85, capsize=4, edgecolor='black', linewidth=0.8)
    
    ax.set_yticks(y_pos)
    if i == 0:
        ax.set_yticklabels(traits_order, fontweight='semibold')
    ax.set_xlabel(r'SNP Heritability ($h^2 \pm 95\%\text{ CI}$)')
    ax.set_title(f"SNP Heritability - {trans}", fontweight='bold', pad=12)
    ax.grid(axis='x', linestyle='--', alpha=0.5)
    ax.legend(loc='lower right', frameon=True)
    ax.set_xlim(0, 0.32)

plt.tight_layout()
fig.savefig(os.path.join(figures_dir, "fig1_heritability_comparison.png"), bbox_inches='tight')
fig.savefig(os.path.join(figures_dir, "fig1_heritability_comparison.pdf"), bbox_inches='tight')
plt.close(fig)
print("Fig 1 saved.")

# ==========================================
# FIGURE 2: BETWEEN-SEX GENETIC CORRELATION FOREST PLOT
# ==========================================
fig, ax = plt.subplots(figsize=(10, 6.5))

bs_irnt = df_bs[df_bs['Transformation'] == 'IRNT (quantile)'].copy()
bs_irnt = bs_irnt.sort_values(by='rg', ascending=True).reset_index(drop=True)

y_indices = np.arange(len(bs_irnt))

# Reference line at rg = 1.0 and rg = 0.0
ax.axvline(1.0, color='gray', linestyle='--', linewidth=1.2, label='Complete Concordance ($r_g = 1$)')
ax.axvline(0.0, color='black', linestyle=':', linewidth=1.0, alpha=0.6)

# Error bars and point estimates
for idx, row in bs_irnt.iterrows():
    color = '#D95F02' if row['Trait'] == 'Testosterone' else '#1F78B4'
    ax.errorbar(row['rg'], idx, xerr=row['SE']*1.96, fmt='o', color=color,
                ecolor=color, elinewidth=2.2, capsize=6, markersize=9, zorder=4)
    
    # Annotate values
    p_text = f"P < 1e-15" if row['P_rg_neq_0'] < 1e-15 else f"P = {row['P_rg_neq_0']:.3f}"
    diff_text = f"P_diff = {row['P_diff_from_1']:.2e}" if row['P_diff_from_1'] < 0.05 else "P_diff = NS"
    ax.text(row['rg'], idx + 0.22, f"rg = {row['rg']:.3f} ({p_text})", 
            ha='center', va='bottom', fontsize=9, fontweight='semibold', color=color)

ax.set_yticks(y_indices)
ax.set_yticklabels(bs_irnt['Trait'], fontweight='bold', fontsize=11)
ax.set_xlabel(r'Between-Sex Genetic Correlation ($r_g \pm 95\%\text{ CI}$)', fontsize=12, labelpad=8)
ax.set_title('Between-Sex Genetic Architecture Concordance (Male vs Female)', fontsize=14, fontweight='bold', pad=15)
ax.set_xlim(-0.25, 1.25)
ax.grid(axis='x', linestyle='--', alpha=0.5)

# Highlight box around testosterone
ax.axhspan(-0.5, 0.5, color='#FEE0D2', alpha=0.4, zorder=1)
ax.text(0.65, 0.0, 'Sex-Specific Genetic Regulation\n(Near-Zero Correlation: $r_g = 0.025$)', 
        fontsize=10, fontweight='bold', color='#B30000', va='center')

plt.tight_layout()
fig.savefig(os.path.join(figures_dir, "fig2_between_sex_genetic_correlation.png"), bbox_inches='tight')
fig.savefig(os.path.join(figures_dir, "fig2_between_sex_genetic_correlation.pdf"), bbox_inches='tight')
plt.close(fig)
print("Fig 2 saved.")

# ==========================================
# FIGURE 3: GENETIC CORRELATION HEATMAPS (MALES VS FEMALES)
# ==========================================
traits_list = ['testo', 'HDL', 'LDL', 'TGL', 'GHb', 'SBP', 'DBP']
traits_labels = ['Testosterone', 'HDL', 'LDL', 'Triglycerides', 'HbA1c', 'SBP', 'DBP']

def build_matrix(sex, trans='IRNT'):
    mat = pd.DataFrame(np.eye(len(traits_list)), index=traits_list, columns=traits_list)
    sub = df_pw[(df_pw['Sex'] == sex) & (df_pw['Transformation'] == trans)]
    for _, r in sub.iterrows():
        t1, t2, val = r['Trait_1_Code'], r['Trait_2_Code'], r['rg']
        if t1 in mat.index and t2 in mat.columns:
            mat.loc[t1, t2] = val
            mat.loc[t2, t1] = val
    mat.index = traits_labels
    mat.columns = traits_labels
    return mat

mat_m = build_matrix('Male', 'IRNT')
mat_f = build_matrix('Female', 'IRNT')

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))

sns.heatmap(mat_m, annot=True, fmt=".2f", cmap="coolwarm", vmin=-0.7, vmax=1.0, 
            cbar=False, ax=ax1, square=True, linewidths=1, linecolor='white')
ax1.set_title("Male Pairwise Genetic Correlations ($r_g$)", fontweight='bold', fontsize=14, pad=12)

sns.heatmap(mat_f, annot=True, fmt=".2f", cmap="coolwarm", vmin=-0.7, vmax=1.0, 
            cbar_kws={'label': 'Genetic Correlation ($r_g$)'}, ax=ax2, square=True, linewidths=1, linecolor='white')
ax2.set_title("Female Pairwise Genetic Correlations ($r_g$)", fontweight='bold', fontsize=14, pad=12)

plt.tight_layout()
fig.savefig(os.path.join(figures_dir, "fig3_genetic_correlation_matrices.png"), bbox_inches='tight')
fig.savefig(os.path.join(figures_dir, "fig3_genetic_correlation_matrices.pdf"), bbox_inches='tight')
plt.close(fig)
print("Fig 3 saved.")

# ==========================================
# FIGURE 4: TESTOSTERONE VS METABOLIC TRAITS (SEX DIVERGENCE)
# ==========================================
testo_sub = df_pw[(df_pw['Trait_1_Code'] == 'testo') & (df_pw['Transformation'] == 'IRNT')].copy()
metabolic_order = ["HDL Cholesterol", "Triglycerides", "Glycated Haemoglobin (HbA1c)", "LDL Cholesterol", "Diastolic Blood Pressure", "Systolic Blood Pressure"]

testo_m = testo_sub[testo_sub['Sex'] == 'Male'].set_index('Trait_2').reindex(metabolic_order)
testo_f = testo_sub[testo_sub['Sex'] == 'Female'].set_index('Trait_2').reindex(metabolic_order)

fig, ax = plt.subplots(figsize=(11, 6))

y_pos = np.arange(len(metabolic_order))
width = 0.35

rects_m = ax.barh(y_pos - width/2, testo_m['rg'], width, xerr=testo_m['SE']*1.96,
                  label='Male Testosterone', color=male_color, alpha=0.9, capsize=4, edgecolor='black', linewidth=0.8)
rects_f = ax.barh(y_pos + width/2, testo_f['rg'], width, xerr=testo_f['SE']*1.96,
                  label='Female Testosterone', color=female_color, alpha=0.9, capsize=4, edgecolor='black', linewidth=0.8)

ax.axvline(0, color='black', linestyle='-', linewidth=0.8)
ax.set_yticks(y_pos)
ax.set_yticklabels(metabolic_order, fontweight='bold', fontsize=11)
ax.set_xlabel(r'Genetic Correlation with Testosterone ($r_g \pm 95\%\text{ CI}$)', fontsize=12, labelpad=8)
ax.set_title('Divergent Genetic Relationships: Testosterone vs. Metabolic Traits by Sex', fontsize=14, fontweight='bold', pad=15)
ax.legend(loc='lower left', frameon=True, fontsize=11)
ax.grid(axis='x', linestyle='--', alpha=0.5)
ax.set_xlim(-0.55, 0.45)

plt.tight_layout()
fig.savefig(os.path.join(figures_dir, "fig4_testosterone_metabolic_cross_correlation.png"), bbox_inches='tight')
fig.savefig(os.path.join(figures_dir, "fig4_testosterone_metabolic_cross_correlation.pdf"), bbox_inches='tight')
plt.close(fig)
print("Fig 4 saved.")

# ==========================================
# FIGURE 5: CROSS-SEX CROSS-TRAIT HEATMAP
# ==========================================
# Male Testosterone vs Female Metabolic Traits AND Female Testosterone vs Male Metabolic Traits
cross_irnt = df_cross[df_cross['Transformation'] == 'IRNT'].copy()
matrix_data = {
    'Male Testosterone vs Female Traits': cross_irnt[cross_irnt['Sex_1'] == 'Male'].set_index('Metabolic_Trait')['rg'],
    'Female Testosterone vs Male Traits': cross_irnt[cross_irnt['Sex_1'] == 'Female'].set_index('Metabolic_Trait')['rg']
}
df_cross_mat = pd.DataFrame(matrix_data).reindex(metabolic_order)

fig, ax = plt.subplots(figsize=(8, 6))
sns.heatmap(df_cross_mat, annot=True, fmt=".3f", cmap="coolwarm", vmin=-0.4, vmax=0.4,
            linewidths=1.2, linecolor='white', cbar_kws={'label': 'Cross-Sex Genetic Correlation ($r_g$)'}, ax=ax)
ax.set_title("Cross-Sex Cross-Trait Genetic Correlations\n(Testosterone vs. Opposite-Sex Metabolic Phenotypes)", fontweight='bold', fontsize=13, pad=15)
ax.set_ylabel("Metabolic Trait", fontweight='semibold')

plt.tight_layout()
fig.savefig(os.path.join(figures_dir, "fig5_cross_sex_cross_trait_heatmap.png"), bbox_inches='tight')
fig.savefig(os.path.join(figures_dir, "fig5_cross_sex_cross_trait_heatmap.pdf"), bbox_inches='tight')
plt.close(fig)
print("Fig 5 saved.")

print("All publication figures successfully created!")
