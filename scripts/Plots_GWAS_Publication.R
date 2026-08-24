# ==============================================================================
# Autonomous BDSA Project: Sex differences in Testosterone & Metabolic Traits
# Publication Figures Script in R (ggplot2, patchwork, readr, dplyr, tidyr)
# ==============================================================================

suppressPackageStartupMessages({
  if (!requireNamespace("readr", quietly = TRUE)) install.packages("readr", repos="https://cloud.r-project.org")
  if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr", repos="https://cloud.r-project.org")
  if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr", repos="https://cloud.r-project.org")
  if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2", repos="https://cloud.r-project.org")
  if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork", repos="https://cloud.r-project.org")
  
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

base_dir <- "c:/Users/enosh/Desktop/BDSA PROJECT"
tables_dir <- file.path(base_dir, "results/tables")
figures_dir <- file.path(base_dir, "figures")

# Load Summary Tables
h2_df <- read_csv(file.path(tables_dir, "snp_heritability_summary.csv"), show_col_types = FALSE)
bs_df <- read_csv(file.path(tables_dir, "between_sex_genetic_correlation.csv"), show_col_types = FALSE)
pw_df <- read_csv(file.path(tables_dir, "within_sex_pairwise_genetic_correlation.csv"), show_col_types = FALSE)

# ------------------------------------------------------------------------------
# 1. Heritability Forest Plot (R Version)
# ------------------------------------------------------------------------------
h2_irnt <- h2_df %>% filter(grepl("IRNT", Transformation))

plot_h2_r <- ggplot(h2_irnt, aes(x = h2, y = reorder(Trait, h2), color = Sex, shape = Sex)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbarh(aes(xmin = h2 - 1.96 * h2_SE, xmax = h2 + 1.96 * h2_SE),
                 position = position_dodge(width = 0.5), height = 0.25, linewidth = 0.8) +
  scale_color_manual(values = c("Female" = "#E64B35", "Male" = "#4DBBD5")) +
  theme_classic(base_size = 12) +
  labs(
    title = "SNP Heritability (h2) across Sexes (IRNT Phenotypes)",
    subtitle = "LD Score Regression (1000G EUR Baseline)",
    x = expression(paste("Observed SNP Heritability (", h^2, " ± 95% CI)")),
    y = "Phenotype"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "gray30"),
    legend.position = "bottom"
  )

ggsave(file.path(figures_dir, "fig1_heritability_R.png"), plot_h2_r, width = 8, height = 6, dpi = 300)

# ------------------------------------------------------------------------------
# 2. Between-Sex Genetic Correlation (R Version)
# ------------------------------------------------------------------------------
bs_irnt <- bs_df %>% filter(grepl("IRNT", Transformation))

plot_bs_r <- ggplot(bs_irnt, aes(x = rg, y = reorder(Trait, rg))) +
  geom_vline(xintercept = 1.0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0.0, linetype = "dotted", color = "gray70") +
  geom_errorbarh(aes(xmin = rg - 1.96 * SE, xmax = rg + 1.96 * SE),
                 height = 0.2, linewidth = 0.9, color = "#1F78B4") +
  geom_point(size = 3.5, color = "#1F78B4") +
  theme_classic(base_size = 12) +
  labs(
    title = "Between-Sex Genetic Correlation (Male vs. Female)",
    x = expression(paste("Genetic Correlation (", r[g], " ± 95% CI)")),
    y = "Phenotype"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

ggsave(file.path(figures_dir, "fig2_between_sex_rg_R.png"), plot_bs_r, width = 8, height = 5.5, dpi = 300)

message("R publication figures generated successfully!")
