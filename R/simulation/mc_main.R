# =============================================================================
# mc_main.R
# Phase 4 Final — Table 3 & Table 4 생성
#
# Table 3: Drug Ranking Summary
#   약물별 PK discontinuation stress 지표 요약 (median [IQR])
#   컬럼: Drug | Css (ng/mL) | t50% (h) | t25% (h) |
#          Rel.dCdt @24h (/h) | AUC deficit 7d (mg/L·h) | Stress Rank
#
# Table 4: CYP2D6 Subgroup Analysis
#   CYP2D6 적용 약물(Fluoxetine, Paroxetine, Venlafaxine IR/XR)에서
#   phenotype별 t50% 및 AUC deficit 비교
#
# Output:
#   outputs/tables/table3_drug_ranking.csv
#   outputs/tables/table4_cyp2d6_subgroup.csv
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

BASE      <- "~/nari-research/pkpd-antidepressant-sim"
OUT_DIR   <- file.path(BASE, "outputs", "mc")
TABLE_DIR <- file.path(BASE, "outputs", "tables")
dir.create(TABLE_DIR, showWarnings = FALSE, recursive = TRUE)

N <- 10000

# ── 1. Load results ───────────────────────────────────────────────────────────
results <- read_csv(
  file.path(OUT_DIR, sprintf("mc_results_all_n%d.csv", N)),
  show_col_types = FALSE
)

# ── 2. Helper: median [Q1–Q3] 형식 문자열 ────────────────────────────────────
fmt_iqr <- function(x, digits = 2) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return("—")
  sprintf("%s [%s–%s]",
          round(median(x), digits),
          round(quantile(x, 0.25), digits),
          round(quantile(x, 0.75), digits))
}

# ── 3. Table 3 — Drug Ranking Summary ────────────────────────────────────────

# Drug 표시 이름
drug_labels <- c(
  "Sertraline"     = "Sertraline (SSRI, IR)",
  "Fluoxetine"     = "Fluoxetine +NFX (SSRI, IR)",
  "Paroxetine"     = "Paroxetine (SSRI, IR)",
  "Venlafaxine_IR" = "Venlafaxine IR (SNRI)",
  "Venlafaxine_XR" = "Venlafaxine XR (SNRI)"
)

# Stress Rank 기준: t50% ascending (빠른 소실 = 높은 stress = 낮은 rank 번호)
# t50% median으로 순위 결정
rank_order <- results %>%
  filter(!is.na(t_50pct)) %>%
  group_by(drug) %>%
  summarise(med_t50 = median(t_50pct), .groups = "drop") %>%
  arrange(med_t50) %>%           # 짧은 t50% = 빠른 소실 = 높은 stress
  mutate(stress_rank = row_number())

table3 <- results %>%
  filter(!is.na(t_50pct), t_50pct >= 1, !is.na(AUC_deficit_7d),
         AUC_deficit_7d > 0) %>%
  group_by(drug) %>%
  summarise(
    n_valid            = n(),
    Css_ngmL           = fmt_iqr(Css * 1000, digits = 1),
    t50_h              = fmt_iqr(t_50pct, digits = 1),
    t25_h              = fmt_iqr(t_25pct, digits = 1),
    rel_dCdt_24h       = fmt_iqr(abs(dCdt_rel_24h), digits = 4),
    AUC_deficit_7d     = fmt_iqr(AUC_deficit_7d, digits = 2),
    .groups            = "drop"
  ) %>%
  left_join(rank_order %>% select(drug, stress_rank), by = "drug") %>%
  arrange(stress_rank) %>%
  mutate(
    Drug = drug_labels[drug],
    `Stress rank` = case_when(
      stress_rank == 1 ~ "1 (highest)",
      stress_rank == 5 ~ "5 (lowest)",
      TRUE             ~ as.character(stress_rank)
    )
  ) %>%
  select(
    Drug,
    `n (valid)`          = n_valid,
    `Css, ng/mL`         = Css_ngmL,
    `t50% C0, h`         = t50_h,
    `t25% C0, h`         = t25_h,
    `|Rel. dC/dt| @24h, h⁻¹` = rel_dCdt_24h,
    `AUC deficit 7d, mg/L·h`  = AUC_deficit_7d,
    `Stress rank`
  )

cat("\n══════════════════════════════════════════════════════════════\n")
cat(" TABLE 3 — Drug Ranking Summary  (median [Q1–Q3])\n")
cat("══════════════════════════════════════════════════════════════\n\n")
print(table3, n = Inf, width = Inf)

# 저장
t3_path <- file.path(TABLE_DIR, "table3_drug_ranking.csv")
write_csv(table3, t3_path)
cat(sprintf("\n✓ Table 3 saved: %s\n", t3_path))
cat(paste0(
  "\nNote: n (valid) < 10000 reflects patients excluded by t50% < 1h filter\n",
  "(extreme rapid metabolizers with near-instantaneous post-discontinuation decline).\n"
))

# ── 4. Table 4 — CYP2D6 Subgroup ─────────────────────────────────────────────
# Sertraline 제외 (CYP2D6 N/A)
# phenotype 순서: poor → normal → rapid

cyp2d6_order <- c("poor", "normal", "rapid")
cyp2d6_labels <- c(
  "poor"   = "Poor (CL × 0.30)",
  "normal" = "Normal (CL × 1.00)",
  "rapid"  = "Rapid (CL × 1.60)"
)

table4 <- results %>%
  filter(
    drug != "Sertraline",
    !is.na(t_50pct), t_50pct >= 1,
    !is.na(AUC_deficit_7d), AUC_deficit_7d > 0,
    cyp2d6 %in% cyp2d6_order
  ) %>%
  mutate(
    cyp2d6 = factor(cyp2d6, levels = cyp2d6_order),
    Drug   = drug_labels[drug]
  ) %>%
  group_by(Drug, cyp2d6) %>%
  summarise(
    n                = n(),
    `t50% C0, h`     = fmt_iqr(t_50pct, digits = 1),
    `AUC deficit 7d` = fmt_iqr(AUC_deficit_7d, digits = 2),
    `|Rel. dC/dt| @24h` = fmt_iqr(abs(dCdt_rel_24h), digits = 4),
    .groups          = "drop"
  ) %>%
  mutate(`CYP2D6 phenotype` = cyp2d6_labels[as.character(cyp2d6)]) %>%
  select(Drug, `CYP2D6 phenotype`, n,
         `t50% C0, h`, `AUC deficit 7d`, `|Rel. dC/dt| @24h`) %>%
  arrange(Drug, match(str_extract(`CYP2D6 phenotype`, "^\\w+"),
                      c("Poor", "Normal", "Rapid")))

cat("\n══════════════════════════════════════════════════════════════\n")
cat(" TABLE 4 — CYP2D6 Subgroup Analysis  (median [Q1–Q3])\n")
cat("══════════════════════════════════════════════════════════════\n\n")
print(table4, n = Inf, width = Inf)

t4_path <- file.path(TABLE_DIR, "table4_cyp2d6_subgroup.csv")
write_csv(table4, t4_path)
cat(sprintf("\n✓ Table 4 saved: %s\n", t4_path))

# ── 5. 완료 요약 ──────────────────────────────────────────────────────────────
cat("\n══════════════════════════════════════════════════════════════\n")
cat(" PHASE 4 COMPLETE\n")
cat("══════════════════════════════════════════════════════════════\n")
cat(sprintf("  Table 3: %d drugs × 7 metrics\n", nrow(table3)))
cat(sprintf("  Table 4: %d drug × phenotype combinations\n", nrow(table4)))
cat("\n  Output files:\n")
cat(sprintf("  %s\n", t3_path))
cat(sprintf("  %s\n", t4_path))
cat("\n▶ 다음 단계: Phase 5 통계 분석\n\n")