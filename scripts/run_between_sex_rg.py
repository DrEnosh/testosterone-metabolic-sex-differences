import subprocess
import os
import sys

ref_ld = os.path.normpath("reference/eur_w_ld_chr/")
w_ld = os.path.normpath("reference/eur_w_ld_chr/")

traits = [
    ("testo", "Testosterone"),
    ("HDL", "HDL Cholesterol"),
    ("LDL", "LDL Cholesterol"),
    ("TGL", "Triglycerides"),
    ("GHb", "Glycated Haemoglobin (HbA1c)"),
    ("SBP", "Systolic Blood Pressure"),
    ("DBP", "Diastolic Blood Pressure")
]

transforms = ["irnt", "raw"]

out_dir = os.path.normpath("results/genetic_correlations/between_sex")
os.makedirs(out_dir, exist_ok=True)

for trait_code, trait_name in traits:
    for trans in transforms:
        p1 = os.path.normpath(f"data/processed/{trait_code}_m_{trans}_v3.sumstats.gz")
        p2 = os.path.normpath(f"data/processed/{trait_code}_f_{trans}_v3.sumstats.gz")
        out_prefix = os.path.join(out_dir, f"{trait_code}_{trans}_M_vs_F")
        
        if not os.path.exists(p1) or not os.path.exists(p2):
            print(f"Missing file for {trait_code} {trans}: {p1} or {p2}")
            continue
            
        cmd = [
            sys.executable,
            os.path.normpath("reference/ldsc/ldsc.py"),
            "--rg", f"{p1},{p2}",
            "--ref-ld-chr", ref_ld + os.sep,
            "--w-ld-chr", w_ld + os.sep,
            "--out", out_prefix
        ]
        
        print(f"Running Between-Sex rg for {trait_name} ({trans})...")
        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode == 0:
            print(f"  -> SUCCESS: {out_prefix}.log generated.")
        else:
            print(f"  -> ERROR: {res.stderr}")

print("All between-sex genetic correlations completed!")
