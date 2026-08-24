import subprocess
import os
import sys

ref_ld = os.path.normpath("reference/eur_w_ld_chr/")
w_ld = os.path.normpath("reference/eur_w_ld_chr/")

traits = ["testo", "HDL", "LDL", "TGL", "GHb", "SBP", "DBP"]
metabolic_traits = ["HDL", "LDL", "TGL", "GHb", "SBP", "DBP"]
transforms = ["irnt", "raw"]

out_dir_cross = os.path.normpath("results/genetic_correlations/cross_sex_cross_trait")
out_dir_pairwise = os.path.normpath("results/genetic_correlations/pairwise")
os.makedirs(out_dir_cross, exist_ok=True)
os.makedirs(out_dir_pairwise, exist_ok=True)

# 1. Male Testosterone vs Female Metabolic Traits
for trans in transforms:
    testo_m = os.path.normpath(f"data/processed/testo_m_{trans}_v3.sumstats.gz")
    for meta in metabolic_traits:
        meta_f = os.path.normpath(f"data/processed/{meta}_f_{trans}_v3.sumstats.gz")
        out_prefix = os.path.join(out_dir_cross, f"testo_M_vs_{meta}_F_{trans}")
        
        cmd = [
            sys.executable,
            os.path.normpath("reference/ldsc/ldsc.py"),
            "--rg", f"{testo_m},{meta_f}",
            "--ref-ld-chr", ref_ld + os.sep,
            "--w-ld-chr", w_ld + os.sep,
            "--out", out_prefix
        ]
        print(f"Running Male Testo vs Female {meta} ({trans})...")
        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode == 0:
            print(f"  -> SUCCESS: {out_prefix}.log")
        else:
            print(f"  -> ERROR: {res.stderr}")

# 2. Female Testosterone vs Male Metabolic Traits
for trans in transforms:
    testo_f = os.path.normpath(f"data/processed/testo_f_{trans}_v3.sumstats.gz")
    for meta in metabolic_traits:
        meta_m = os.path.normpath(f"data/processed/{meta}_m_{trans}_v3.sumstats.gz")
        out_prefix = os.path.join(out_dir_cross, f"testo_F_vs_{meta}_M_{trans}")
        
        cmd = [
            sys.executable,
            os.path.normpath("reference/ldsc/ldsc.py"),
            "--rg", f"{testo_f},{meta_m}",
            "--ref-ld-chr", ref_ld + os.sep,
            "--w-ld-chr", w_ld + os.sep,
            "--out", out_prefix
        ]
        print(f"Running Female Testo vs Male {meta} ({trans})...")
        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode == 0:
            print(f"  -> SUCCESS: {out_prefix}.log")
        else:
            print(f"  -> ERROR: {res.stderr}")

# 3. Full 7x7 matrix for Males (IRNT & RAW) and Females (IRNT & RAW)
for sex in ["m", "f"]:
    for trans in transforms:
        file_list = [os.path.normpath(f"data/processed/{t}_{sex}_{trans}_v3.sumstats.gz") for t in traits]
        rg_arg = ",".join(file_list)
        out_prefix = os.path.join(out_dir_pairwise, f"all_traits_{sex}_{trans}")
        
        cmd = [
            sys.executable,
            os.path.normpath("reference/ldsc/ldsc.py"),
            "--rg", rg_arg,
            "--ref-ld-chr", ref_ld + os.sep,
            "--w-ld-chr", w_ld + os.sep,
            "--out", out_prefix
        ]
        print(f"Running Full 7-trait Matrix for Sex={sex}, Trans={trans}...")
        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode == 0:
            print(f"  -> SUCCESS: {out_prefix}.log")
        else:
            print(f"  -> ERROR: {res.stderr}")

print("All cross-sex-cross-trait and pairwise genetic correlations completed!")
