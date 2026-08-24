import os
import glob
import re
import pandas as pd
import numpy as np

# Base directories
base_dir = r"c:\Users\enosh\Desktop\BDSA PROJECT"
h2_dir = os.path.join(base_dir, "results", "heritability")
gencor_between_sex_dir = os.path.join(base_dir, "results", "genetic_correlations", "between_sex")
gencor_cross_dir = os.path.join(base_dir, "results", "genetic_correlations", "cross_sex_cross_trait")
gencor_pairwise_dir = os.path.join(base_dir, "results", "genetic_correlations", "pairwise")
gencor_base_dir = os.path.join(base_dir, "results", "genetic_correlations")
tables_dir = os.path.join(base_dir, "results", "tables")
os.makedirs(tables_dir, exist_ok=True)

# Helper function to parse h2 log
def parse_h2_log(filepath):
    with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
        text = f.read()
    
    # Extract Total Observed scale h2: 0.1877 (0.0242)
    h2_match = re.search(r'Total Observed scale h2:\s*([\d.]+)\s*\(([\d.]+)\)', text)
    lambda_match = re.search(r'Lambda GC:\s*([\d.]+)', text)
    mean_chi2_match = re.search(r'Mean Chi\^2:\s*([\d.]+)', text)
    intercept_match = re.search(r'Intercept:\s*([\d.]+)\s*\(([\d.]+)\)', text)
    ratio_match = re.search(r'Ratio:\s*([-\d.]+)\s*\(([\d.]+)\)', text)
    
    if h2_match:
        h2 = float(h2_match.group(1))
        h2_se = float(h2_match.group(2))
        lam = float(lambda_match.group(1)) if lambda_match else np.nan
        mean_chi2 = float(mean_chi2_match.group(1)) if mean_chi2_match else np.nan
        intercept = float(intercept_match.group(1)) if intercept_match else np.nan
        intercept_se = float(intercept_match.group(2)) if intercept_match else np.nan
        ratio = float(ratio_match.group(1)) if ratio_match else np.nan
        ratio_se = float(ratio_match.group(2)) if ratio_match else np.nan
        
        return {
            'h2': h2,
            'h2_se': h2_se,
            'lambda_gc': lam,
            'mean_chi2': mean_chi2,
            'intercept': intercept,
            'intercept_se': intercept_se,
            'ratio': ratio,
            'ratio_se': ratio_se
        }
    return None

# Helper function to parse summary table from rg log
def parse_rg_log(filepath):
    with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
        text = f.read()
    
    if "Summary of Genetic Correlation Results" not in text:
        return []
    
    summary_part = text.split("Summary of Genetic Correlation Results")[1]
    lines = summary_part.strip().splitlines()
    if not lines:
        return []
    
    header_line = lines[0]
    data_lines = [l for l in lines[1:] if l.strip() and not l.startswith("Analysis finished") and not l.startswith("Total time")]
    
    results = []
    for dl in data_lines:
        parts = dl.split()
        if len(parts) >= 6:
            p1 = os.path.basename(parts[0]).replace(".sumstats.gz", "")
            p2 = os.path.basename(parts[1]).replace(".sumstats.gz", "")
            try:
                rg = float(parts[2])
                se = float(parts[3])
                z = float(parts[4])
                p = float(parts[5])
                h2_obs = float(parts[6]) if len(parts) > 6 else np.nan
                h2_obs_se = float(parts[7]) if len(parts) > 7 else np.nan
                gcov_int = float(parts[10]) if len(parts) > 10 else np.nan
                
                results.append({
                    'p1': p1,
                    'p2': p2,
                    'rg': rg,
                    'se': se,
                    'z': z,
                    'p': p,
                    'h2_obs': h2_obs,
                    'h2_obs_se': h2_obs_se,
                    'gcov_int': gcov_int,
                    'source_log': os.path.basename(filepath)
                })
            except Exception as e:
                continue
    return results

print("Parser functions compiled.")
