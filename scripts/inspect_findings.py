import pandas as pd
import os

tables_dir = r"c:\Users\enosh\Desktop\BDSA PROJECT\results\tables"

df_h2 = pd.read_csv(os.path.join(tables_dir, "snp_heritability_summary.csv"))
df_bs = pd.read_csv(os.path.join(tables_dir, "between_sex_genetic_correlation.csv"))
df_pw = pd.read_csv(os.path.join(tables_dir, "within_sex_pairwise_genetic_correlation.csv"))
df_cross = pd.read_csv(os.path.join(tables_dir, "cross_sex_cross_trait_correlation.csv"))

print("=== 1. BETWEEN-SEX GENETIC CORRELATIONS (MALE vs FEMALE) ===")
print(df_bs[['Trait', 'Transformation', 'rg', 'SE', 'Z', 'P_rg_neq_0', 'P_diff_from_1']].to_string(index=False))

print("\n=== 2. SNP HERITABILITY SUMMARY ===")
print(df_h2[['Trait', 'Sex', 'Transformation', 'h2', 'h2_SE', 'Lambda_GC', 'Intercept', 'Ratio']].to_string(index=False))

print("\n=== 3. TESTOSTERONE VS METABOLIC TRAITS WITHIN-SEX ===")
testo_pw = df_pw[df_pw['Trait_1_Code'] == 'testo'][['Trait_2', 'Sex', 'Transformation', 'rg', 'SE', 'P']]
print(testo_pw.sort_values(by=['Transformation', 'Trait_2', 'Sex']).to_string(index=False))

