import os
import glob
import re
import pandas as pd
import numpy as np

base_dir = r"c:\Users\enosh\Desktop\BDSA PROJECT"
tables_dir = os.path.join(base_dir, "results", "tables")
os.makedirs(tables_dir, exist_ok=True)

# Trait display names mapping
TRAIT_MAP = {
    'testo': 'Testosterone',
    'HDL': 'HDL Cholesterol',
    'LDL': 'LDL Cholesterol',
    'TGL': 'Triglycerides',
    'GHb': 'Glycated Haemoglobin (HbA1c)',
    'SBP': 'Systolic Blood Pressure',
    'DBP': 'Diastolic Blood Pressure'
}

# 1. PARSE HERITABILITY LOGS
h2_dir = os.path.join(base_dir, "results", "heritability")
h2_rows = []

for log_file in glob.glob(os.path.join(h2_dir, "*.log")):
    fname = os.path.basename(log_file)
    with open(log_file, 'r', encoding='utf-8', errors='ignore') as f:
        text = f.read()
    
    h2_match = re.search(r'Total Observed scale h2:\s*([\d.]+)\s*\(([\d.]+)\)', text)
    lambda_match = re.search(r'Lambda GC:\s*([\d.]+)', text)
    mean_chi2_match = re.search(r'Mean Chi\^2:\s*([\d.]+)', text)
    intercept_match = re.search(r'Intercept:\s*([\d.]+)\s*\(([\d.]+)\)', text)
    ratio_match = re.search(r'Ratio:\s*([-\d.]+)\s*\(([\d.]+)\)', text)
    snps_match = re.search(r'After merging with regression SNP LD,\s*(\d+)\s*SNPs remain', text)
    
    if h2_match:
        # Parse trait, sex, transform from filename (e.g. testo_m_irnt_v3.log)
        parts = fname.replace(".log", "").replace("_v3", "").split("_")
        trait_code = parts[0]
        sex_code = parts[1] if len(parts) > 1 else 'unknown'
        trans_code = parts[2] if len(parts) > 2 else 'unknown'
        
        sex = "Female" if sex_code == 'f' else ("Male" if sex_code == 'm' else "Both")
        trans = "IRNT (quantile)" if trans_code == 'irnt' else "Raw (untransformed)"
        trait_name = TRAIT_MAP.get(trait_code, trait_code)
        
        h2 = float(h2_match.group(1))
        h2_se = float(h2_match.group(2))
        z_score = h2 / h2_se if h2_se > 0 else np.nan
        # 2-sided p-value from normal approximation
        from scipy.stats import norm
        pval = 2 * (1 - norm.cdf(abs(z_score))) if not np.isnan(z_score) else np.nan
        
        ci_lower = max(0.0, h2 - 1.96 * h2_se)
        ci_upper = h2 + 1.96 * h2_se
        
        h2_rows.append({
            'Trait_Code': trait_code,
            'Trait': trait_name,
            'Sex': sex,
            'Transformation': trans,
            'SNPs': int(snps_match.group(1)) if snps_match else np.nan,
            'h2': h2,
            'h2_SE': h2_se,
            'h2_95CI_Lower': ci_lower,
            'h2_95CI_Upper': ci_upper,
            'Z_Score': z_score,
            'P_Value': pval,
            'Lambda_GC': float(lambda_match.group(1)) if lambda_match else np.nan,
            'Mean_Chi2': float(mean_chi2_match.group(1)) if mean_chi2_match else np.nan,
            'Intercept': float(intercept_match.group(1)) if intercept_match else np.nan,
            'Intercept_SE': float(intercept_match.group(2)) if intercept_match else np.nan,
            'Ratio': float(ratio_match.group(1)) if ratio_match and float(ratio_match.group(1)) != -1 else np.nan,
            'Ratio_SE': float(ratio_match.group(2)) if ratio_match else np.nan,
            'Filename': fname
        })

df_h2 = pd.DataFrame(h2_rows).sort_values(by=['Trait', 'Transformation', 'Sex'])
h2_csv_path = os.path.join(tables_dir, "snp_heritability_summary.csv")
df_h2.to_csv(h2_csv_path, index=False)
print(f"Saved {len(df_h2)} heritability records to {h2_csv_path}")

# 2. PARSE BETWEEN-SEX GENETIC CORRELATION LOGS
between_sex_dir = os.path.join(base_dir, "results", "genetic_correlations", "between_sex")
bs_rows = []

for log_file in glob.glob(os.path.join(between_sex_dir, "*.log")):
    fname = os.path.basename(log_file)
    with open(log_file, 'r', encoding='utf-8', errors='ignore') as f:
        text = f.read()
    
    if "Summary of Genetic Correlation Results" in text:
        summary_part = text.split("Summary of Genetic Correlation Results")[1].strip().splitlines()
        if len(summary_part) > 1:
            data_line = summary_part[1].split()
            if len(data_line) >= 6:
                rg = float(data_line[2])
                se = float(data_line[3])
                z = float(data_line[4])
                p = float(data_line[5])
                h2_m = float(data_line[6]) if len(data_line) > 6 else np.nan
                h2_m_se = float(data_line[7]) if len(data_line) > 7 else np.nan
                
                parts = fname.replace(".log", "").split("_")
                trait_code = parts[0]
                trans_code = parts[1]
                
                trans = "IRNT (quantile)" if trans_code == 'irnt' else "Raw (untransformed)"
                trait_name = TRAIT_MAP.get(trait_code, trait_code)
                
                # Test for sex difference: H0: rg = 1. Z_diff = (1 - rg) / se
                z_diff = (1.0 - rg) / se if se > 0 else np.nan
                from scipy.stats import norm
                p_diff = 2 * (1 - norm.cdf(abs(z_diff))) if not np.isnan(z_diff) else np.nan
                
                bs_rows.append({
                    'Trait_Code': trait_code,
                    'Trait': trait_name,
                    'Transformation': trans,
                    'rg': rg,
                    'SE': se,
                    'rg_95CI_Lower': rg - 1.96 * se,
                    'rg_95CI_Upper': min(1.0, rg + 1.96 * se),
                    'Z': z,
                    'P_rg_neq_0': p,
                    'Z_diff_from_1': z_diff,
                    'P_diff_from_1': p_diff,
                    'h2_Female_Obs': h2_m,
                    'h2_Female_SE': h2_m_se,
                    'Filename': fname
                })

df_bs = pd.DataFrame(bs_rows).sort_values(by=['Transformation', 'rg'], ascending=[True, False])
bs_csv_path = os.path.join(tables_dir, "between_sex_genetic_correlation.csv")
df_bs.to_csv(bs_csv_path, index=False)
print(f"Saved {len(df_bs)} between-sex genetic correlations to {bs_csv_path}")

# 3. PARSE CROSS-SEX CROSS-TRAIT LOGS
cross_dir = os.path.join(base_dir, "results", "genetic_correlations", "cross_sex_cross_trait")
cross_rows = []

for log_file in glob.glob(os.path.join(cross_dir, "*.log")):
    fname = os.path.basename(log_file)
    with open(log_file, 'r', encoding='utf-8', errors='ignore') as f:
        text = f.read()
    
    if "Summary of Genetic Correlation Results" in text:
        summary_part = text.split("Summary of Genetic Correlation Results")[1].strip().splitlines()
        if len(summary_part) > 1:
            data_line = summary_part[1].split()
            if len(data_line) >= 6:
                rg = float(data_line[2])
                se = float(data_line[3])
                z = float(data_line[4])
                p = float(data_line[5])
                
                # Filename: testo_M_vs_HDL_F_irnt.log or testo_F_vs_HDL_M_irnt.log
                fname_clean = fname.replace(".log", "")
                parts = fname_clean.split("_")
                t1 = parts[0]
                s1 = parts[1]
                t2 = parts[3]
                s2 = parts[4]
                trans = parts[5]
                
                cross_rows.append({
                    'Trait_1': f"{TRAIT_MAP.get(t1, t1)} ({s1})",
                    'Trait_2': f"{TRAIT_MAP.get(t2, t2)} ({s2})",
                    'Sex_1': 'Male' if s1 == 'M' else 'Female',
                    'Sex_2': 'Female' if s2 == 'F' else 'Male',
                    'Metabolic_Trait': TRAIT_MAP.get(t2, t2),
                    'Transformation': "IRNT" if trans == 'irnt' else "Raw",
                    'rg': rg,
                    'SE': se,
                    'rg_95CI_Lower': rg - 1.96 * se,
                    'rg_95CI_Upper': rg + 1.96 * se,
                    'Z': z,
                    'P': p,
                    'Filename': fname
                })

df_cross = pd.DataFrame(cross_rows).sort_values(by=['Transformation', 'Sex_1', 'Metabolic_Trait'])
cross_csv_path = os.path.join(tables_dir, "cross_sex_cross_trait_correlation.csv")
df_cross.to_csv(cross_csv_path, index=False)
print(f"Saved {len(df_cross)} cross-sex cross-trait genetic correlations to {cross_csv_path}")

# 4. PARSE WITHIN-SEX PAIRWISE LOGS
pairwise_rows = []
for p_dir in [os.path.join(base_dir, "results", "genetic_correlations", "pairwise"),
              os.path.join(base_dir, "results", "genetic_correlations")]:
    for log_file in glob.glob(os.path.join(p_dir, "*.log")):
        fname = os.path.basename(log_file)
        # Skip between-sex or cross-sex filenames
        if "_M_vs_F" in fname or "_vs_" in fname or fname.startswith("testom"):
            continue
        with open(log_file, 'r', encoding='utf-8', errors='ignore') as f:
            text = f.read()
        
        if "Summary of Genetic Correlation Results" in text:
            summary_part = text.split("Summary of Genetic Correlation Results")[1].strip().splitlines()
            for line in summary_part[1:]:
                parts = line.split()
                if len(parts) >= 6 and not line.startswith("Analysis") and not line.startswith("Total"):
                    p1 = os.path.basename(parts[0]).replace(".sumstats.gz", "").replace("_v3", "")
                    p2 = os.path.basename(parts[1]).replace(".sumstats.gz", "").replace("_v3", "")
                    try:
                        rg = float(parts[2])
                        se = float(parts[3])
                        z = float(parts[4])
                        p = float(parts[5])
                        
                        p1_parts = p1.split("_")
                        p2_parts = p2.split("_")
                        
                        t1_code = p1_parts[0]
                        s1_code = p1_parts[1] if len(p1_parts) > 1 else 'unknown'
                        tr1_code = p1_parts[2] if len(p1_parts) > 2 else 'unknown'
                        
                        t2_code = p2_parts[0]
                        s2_code = p2_parts[1] if len(p2_parts) > 1 else 'unknown'
                        tr2_code = p2_parts[2] if len(p2_parts) > 2 else 'unknown'
                        
                        if s1_code == s2_code:
                            sex_label = "Male" if s1_code == 'm' else ("Female" if s1_code == 'f' else s1_code)
                            trans_label = "IRNT" if tr1_code == 'irnt' else "Raw"
                            
                            pairwise_rows.append({
                                'Trait_1_Code': t1_code,
                                'Trait_1': TRAIT_MAP.get(t1_code, t1_code),
                                'Trait_2_Code': t2_code,
                                'Trait_2': TRAIT_MAP.get(t2_code, t2_code),
                                'Sex': sex_label,
                                'Transformation': trans_label,
                                'rg': rg,
                                'SE': se,
                                'rg_95CI_Lower': rg - 1.96 * se,
                                'rg_95CI_Upper': rg + 1.96 * se,
                                'Z': z,
                                'P': p,
                                'Source_Log': fname
                            })
                    except Exception:
                        continue

df_pw = pd.DataFrame(pairwise_rows).drop_duplicates(subset=['Trait_1_Code', 'Trait_2_Code', 'Sex', 'Transformation'])
pw_csv_path = os.path.join(tables_dir, "within_sex_pairwise_genetic_correlation.csv")
df_pw.to_csv(pw_csv_path, index=False)
print(f"Saved {len(df_pw)} within-sex pairwise genetic correlations to {pw_csv_path}")

print("All summary tables successfully compiled!")
