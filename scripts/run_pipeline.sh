#!/usr/bin/env bash
# ==============================================================================
# Autonomous BDSA Project - Master Pipeline Execution Script (Bash)
# Project: Sex differences in genetic architecture of testosterone & metabolic traits
# ==============================================================================

set -e

echo "=========================================================="
echo "Starting Autonomous BDSA LDSC Analysis Pipeline (Bash)"
echo "=========================================================="

echo "[Step 1/5] Checking Python Dependencies..."
python3 -m pip install -q pandas numpy scipy matplotlib seaborn statsmodels requests urllib3

echo "[Step 2/5] Running Between-Sex LDSC Genetic Correlations..."
python3 scripts/run_between_sex_rg.py

echo "[Step 3/5] Running Cross-Sex & Pairwise Genetic Correlations..."
python3 scripts/run_cross_sex_cross_trait_rg.py

echo "[Step 4/5] Compiling and Structuring Result Tables..."
python3 scripts/compile_all_tables.py

echo "[Step 5/5] Generating Publication-Quality Figures..."
python3 scripts/generate_all_figures.py

echo "=========================================================="
echo "Pipeline Completed Successfully! All results in results/ and figures/"
echo "=========================================================="
