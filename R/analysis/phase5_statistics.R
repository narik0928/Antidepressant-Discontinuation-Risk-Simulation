# =============================================================================
# phase5_statistics.R
# Phase 5 — Statistical Analysis (Option B: PK metrics internal analysis)
#
# Analysis 1: Spearman correlation — t50% vs AUC deficit (per drug + overall)
# Analysis 2: CYP2D6 effect size — Kruskal-Wallis + eta-squared
# Analysis 3: Drug comparison — Dunn test + Bonferroni correction
#
# Output:
#   outputs/tables/table5_correlations.csv
#   outputs/tables/table6_kruskal_cyp2d6.csv
#   outputs/tables/table7_dunn_drug_comparison.csv
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(rstatix)
})

BASE      <- "~/nari-research/pkpd-antidepressant-sim"
OUT_DIR   <- file.path(BASE, "outputs", "mc")
TABLE_DIR <- file.path(BASE, "outputs", "tables")
dir.create(TABLE_DIR, showWarnings = FALSE, recursive = TRUE)

N <- 10000

# ── 1. Load & filter ──────────────────────────────────────────────────────────
results <- read_csv(
  file.path(OUT_DIR, sprintf("mc_results_all_n%d.csv", N)),
  show_col_types = FALSE
) %>%
  filter(
    !is.na(t_50pct),        t_50pct >= 1,
    !is.na(AUC_deficit_7d), AUC_deficit_7d > 0,
    !is.na(dCdt_rel_24h),   dCdt_rel_24h < 0
  ) %>%
  mutate(
    abs_dCdt = abs(dCdt_rel_24h),
    drug     = factor(drug),
    cyp2d6   = factor(cyp2d6, levels = c("poor", "normal", "rapid"))
  )

cat(sprintf("\n총 분석 대상: %d명\n", nrow(results)))

# ── 2. Analysis 1: Spearman Correlation ───────────────────────────────────────
cat("\n══════════════════════════════════════════════════════════════\n")
cat(" ANALYSIS 1: Spearman Correlation — t50% vs AUC deficit\n")
cat("══════════════════════════════════════════════════════════════\n\n")

# 전체
overall_cor <- suppressWarnings(
  cor.test(results$t_50pct, results$AUC_deficit_7d, method = "spearman")
)
cat(sprintf("Overall: rho = %.3f, p %s\n",
            overall_cor$estimate,
            ifelse(overall_cor$p.value < 0.001, "< 0.001",
                   sprintf("= %.3f", overall_cor$p.value))))

# 약물별
drug_cor <- results %>%
  group_by(drug) %>%
  summarise(
    n     = n(),
    rho_t50_AUC = suppressWarnings(
      cor.test(t_50pct, AUC_deficit_7d, method = "spearman")$estimate
    ),
    p_t50_AUC = suppressWarnings(
      cor.test(t_50pct, AUC_deficit_7d, method = "spearman")$p.value
    ),
    rho_dCdt_AUC = suppressWarnings(
      cor.test(abs_dCdt, AUC_deficit_7d, method = "spearman")$estimate
    ),
    p_dCdt_AUC = suppressWarnings(
      cor.test(abs_dCdt, AUC_deficit_7d, method = "spearman")$p.value
    ),
    .groups = "drop"
  ) %>%
  mutate(
    across(starts_with("rho"), ~ round(.x, 3)),
    across(starts_with("p"),   ~ ifelse(.x < 0.001, "< 0.001",
                                        sprintf("%.3f", .x)))
  )

print(drug_cor, n = Inf)

t5_path <- file.path(TABLE_DIR, "table5_correlations.csv")
write_csv(drug_cor, t5_path)
cat(sprintf("\n✓ Table 5 saved: %s\n", t5_path))

# ── 3. Analysis 2: CYP2D6 Kruskal-Wallis + eta-squared ───────────────────────
cat("\n══════════════════════════════════════════════════════════════\n")
cat(" ANALYSIS 2: CYP2D6 Effect — Kruskal-Wallis + eta-squared\n")
cat("══════════════════════════════════════════════════════════════\n\n")

df_cyp <- results %>% filter(drug != "Sertraline")

# t50% 기준 CYP2D6 효과
kw_cyp <- df_cyp %>%
  group_by(drug) %>%
  do(
    kw    = kruskal_test(., t_50pct ~ cyp2d6),
    effsize = kruskal_effsize(., t_50pct ~ cyp2d6)
  ) %>%
  summarise(
    drug       = first(drug),
    H          = round(kw$statistic, 2),
    df         = kw$df,
    p_value    = ifelse(kw$p < 0.001, "< 0.001", sprintf("%.4f", kw$p)),
    eta_sq     = round(effsize$effsize, 3),
    magnitude  = effsize$magnitude,
    .groups    = "drop"
  )

print(kw_cyp, n = Inf)

t6_path <- file.path(TABLE_DIR, "table6_kruskal_cyp2d6.csv")
write_csv(kw_cyp, t6_path)
cat(sprintf("\n✓ Table 6 saved: %s\n", t6_path))

# ── 4. Analysis 3: Drug Comparison — Dunn test ───────────────────────────────
cat("\n══════════════════════════════════════════════════════════════\n")
cat(" ANALYSIS 3: Drug Comparison — Dunn test (Bonferroni)\n")
cat("══════════════════════════════════════════════════════════════\n\n")

# t50% 약물 간 비교
dunn_t50 <- results %>%
  dunn_test(t_50pct ~ drug, p.adjust.method = "bonferroni") %>%
  select(group1, group2, statistic, p, p.adj, p.adj.signif) %>%
  mutate(
    metric    = "t50% C0 (h)",
    statistic = round(statistic, 3),
    p         = ifelse(p < 0.001, "< 0.001", sprintf("%.4f", p)),
    p.adj     = ifelse(p.adj < 0.001, "< 0.001", sprintf("%.4f", p.adj))
  )

# AUC deficit 약물 간 비교
dunn_AUC <- results %>%
  dunn_test(AUC_deficit_7d ~ drug, p.adjust.method = "bonferroni") %>%
  select(group1, group2, statistic, p, p.adj, p.adj.signif) %>%
  mutate(
    metric    = "AUC deficit 7d",
    statistic = round(statistic, 3),
    p         = ifelse(p < 0.001, "< 0.001", sprintf("%.4f", p)),
    p.adj     = ifelse(p.adj < 0.001, "< 0.001", sprintf("%.4f", p.adj))
  )

dunn_all <- bind_rows(dunn_t50, dunn_AUC) %>%
  arrange(metric, group1, group2)

cat("── t50% pairwise comparisons ──\n")
print(dunn_t50 %>% select(group1, group2, statistic, p.adj, p.adj.signif),
      n = Inf)

cat("\n── AUC deficit pairwise comparisons ──\n")
print(dunn_AUC %>% select(group1, group2, statistic, p.adj, p.adj.signif),
      n = Inf)

t7_path <- file.path(TABLE_DIR, "table7_dunn_drug_comparison.csv")
write_csv(dunn_all, t7_path)
cat(sprintf("\n✓ Table 7 saved: %s\n", t7_path))

# ── 5. Summary ────────────────────────────────────────────────────────────────
cat("\n══════════════════════════════════════════════════════════════\n")
cat(" PHASE 5 COMPLETE\n")
cat("══════════════════════════════════════════════════════════════\n")
cat(sprintf("  Table 5: Spearman correlations (%d drug-level pairs)\n",
            nrow(drug_cor)))
cat(sprintf("  Table 6: CYP2D6 Kruskal-Wallis (%d drugs)\n",
            nrow(kw_cyp)))
cat(sprintf("  Table 7: Dunn test (%d pairwise comparisons × 2 metrics)\n",
            nrow(dunn_all)))
cat("\n  Saved to: outputs/tables/\n\n")