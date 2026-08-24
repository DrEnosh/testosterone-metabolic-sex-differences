# BDSA Project: Sex Differences in the Genetic Architecture of Testosterone and Metabolic Traits in Humans

An autonomous bioinformatics and statistical genetics pipeline utilizing **Linkage Disequilibrium Score Regression (LDSC)** to evaluate sex-stratified genetic architecture, SNP heritability ($h^2$), and cross-trait / cross-sex genetic correlations ($r_g$).

---

## 📌 Objectives
1. **SNP Heritability ($h^2$):** Estimate sex-stratified SNP heritability for 7 key phenotypes (Testosterone, HDL-C, LDL-C, Triglycerides, HbA1c, SBP, DBP) across transformations (`IRNT` and `Raw`).
2. **Between-Sex Genetic Correlation:** Quantify polygenic concordance ($r_g$) between males and females for each phenotype.
3. **Between-Trait Genetic Correlation:** Evaluate pairwise genetic correlations among all cardiometabolic traits separately in males and females.
4. **Between-Sex-Between-Trait Genetic Correlation:** Investigate how testosterone in one sex correlates with cardiometabolic phenotypes in the opposite sex.

---

## 🔬 Key Scientific Findings
- **Near-Zero Between-Sex Genetic Correlation for Testosterone:** $r_g = 0.0249 \pm 0.0510$ ($P = 0.626; P_{\text{diff from 1.0}} < 10^{-50}$), demonstrating profound sexual dimorphism and distinct polygenic control in men versus women.
- **High Genetic Concordance for Cardiometabolic Traits:** All metabolic phenotypes show high genetic stability across sexes (HbA1c $r_g = 0.969$; Triglycerides $r_g = 0.910$; HDL-C $r_g = 0.903$; LDL-C $r_g = 0.819$; SBP $r_g = 0.814$; DBP $r_g = 0.673$).
- **Sex-Divergent Metabolic Impact of Testosterone:**
  - **In Males:** Genetically higher testosterone is strongly metabolically protective (HDL: $r_g = +0.294, P = 1.38 \times 10^{-23}$; Triglycerides: $r_g = -0.369, P = 1.71 \times 10^{-19}$; HbA1c: $r_g = -0.238, P = 2.05 \times 10^{-15}$).
  - **In Females:** These protective relationships are absent ($r_g \approx 0$).

---

## 📂 Repository Structure
```
BDSA PROJECT/
├── data/
│   ├── raw/                               # Raw GWAS files
│   └── processed/                         # Munged .sumstats.gz GWAS files
│
├── reference/
│   ├── eur_w_ld_chr/                      # 1000G European LD scores (chrs 1-22)
│   └── ldsc/                              # Python 3 LDSC codebase
│
├── scripts/
│   ├── run_between_sex_rg.py              # Between-sex genetic correlation pipeline
│   ├── run_cross_sex_cross_trait_rg.py    # Cross-sex & pairwise 7x7 correlation pipeline
│   ├── compile_all_tables.py              # Log parser and summary CSV generator
│   ├── generate_all_figures.py            # Publication figure generator (PNG + PDF)
│   ├── Plots_GWAS_Publication.R           # R ggplot2 script for figure reproducibility
│   ├── run_pipeline.ps1                   # Master PowerShell execution script
│   └── run_pipeline.sh                    # Master Bash execution script
│
├── results/
│   ├── heritability/                      # 31 h2 logs & summaries
│   ├── genetic_correlations/              # LDSC correlation logs (between-sex & pairwise)
│   ├── tables/
│   │   ├── snp_heritability_summary.csv
│   │   ├── between_sex_genetic_correlation.csv
│   │   ├── within_sex_pairwise_genetic_correlation.csv
│   │   └── cross_sex_cross_trait_correlation.csv
│   └── reports/
│       └── BDSA_PROJECT_FINAL_REPORT.md   # Complete scientific final report
│
└── figures/
    ├── fig1_heritability_comparison.png (.pdf)
    ├── fig2_between_sex_genetic_correlation.png (.pdf)
    ├── fig3_genetic_correlation_matrices.png (.pdf)
    ├── fig4_testosterone_metabolic_cross_correlation.png (.pdf)
    └── fig5_cross_sex_cross_trait_heatmap.png (.pdf)
```

---

## 🚀 Execution & Reproducibility

### PowerShell (Windows):
```powershell
powershell -ExecutionPolicy Bypass -File scripts/run_pipeline.ps1
```

### Bash (Linux / macOS / WSL):
```bash
bash scripts/run_pipeline.sh
```

---

## 📊 Summary Tables
All output tables are available in `results/tables/`:
- `snp_heritability_summary.csv`
- `between_sex_genetic_correlation.csv`
- `within_sex_pairwise_genetic_correlation.csv`
- `cross_sex_cross_trait_correlation.csv`

---

## 📜 Full Report
For the complete methodology, biological discussion, and detailed statistics, refer to [results/reports/BDSA_PROJECT_FINAL_REPORT.md](results/reports/BDSA_PROJECT_FINAL_REPORT.md).
