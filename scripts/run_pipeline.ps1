# ==============================================================================
# Autonomous BDSA Project - Master Pipeline Execution Script (PowerShell)
# Project: Sex differences in genetic architecture of testosterone & metabolic traits
# ==============================================================================

Write-Host "==========================================================" -ForegroundColor Cyan
Write-Host "Starting Autonomous BDSA LDSC Analysis Pipeline" -ForegroundColor Cyan
Write-Host "==========================================================" -ForegroundColor Cyan

# 1. Environment and dependencies
Write-Host "[Step 1/5] Checking Python Environment & Dependencies..." -ForegroundColor Green
python -m pip install -q pandas numpy scipy matplotlib seaborn statsmodels requests urllib3

# 2. Run Between-Sex Genetic Correlations
Write-Host "[Step 2/5] Running Between-Sex LDSC Genetic Correlations..." -ForegroundColor Green
python scripts/run_between_sex_rg.py

# 3. Run Cross-Sex Cross-Trait & Pairwise Correlations
Write-Host "[Step 3/5] Running Cross-Sex & Pairwise Genetic Correlations..." -ForegroundColor Green
python scripts/run_cross_sex_cross_trait_rg.py

# 4. Compile Tables
Write-Host "[Step 4/5] Compiling and Structuring Result Tables..." -ForegroundColor Green
python scripts/compile_all_tables.py

# 5. Generate Publication Figures
Write-Host "[Step 5/5] Generating Publication-Quality Figures..." -ForegroundColor Green
python scripts/generate_all_figures.py

Write-Host "==========================================================" -ForegroundColor Cyan
Write-Host "Pipeline Completed Successfully! Results in results/ and figures/" -ForegroundColor Cyan
Write-Host "==========================================================" -ForegroundColor Cyan
