# =============================================================================
# stability_check.R
# Monte Carlo Stability Check — n=1,000 / 5,000 / 10,000 수렴 확인
#
# 방법: 기존 n=10,000 결과에서 서브샘플링 (재시뮬 불필요)
#       각 n에서 10회 반복 → 핵심 지표 median 분포 확인
#       CV < 5% → 수렴 판정 (논문 Methods 정당화)
#
# Output: outputs/tables/stability_check.csv
#         outputs/figures/stability_check_plot.png
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

BASE      <- "~/nari-research/pkpd-antidepressant-sim"
OUT_DIR   <- file.path(BASE, "outputs", "mc")
TABLE_DIR <- file.path(BASE, "outputs", "tables")
FIG_DIR   <- file.path(BASE, "outputs", "figures")
dir.create(TABLE_DIR, showWarnings = FALSE, recursive = TRUE)

set.seed(2025)

# ── 1. Load full results ──────────────────────────────────────────────────────
results <- read_csv(
  file.path(OUT_DIR, "mc_results_all_n10000.csv"),
  show_col_types = FALSE
) %>%
  filter(!is.na(t_50pct), t_50pct >= 1,
         !is.na(AUC_deficit_7d), AUC_deficit_7d > 0)

# ── 2. Stability check parameters ────────────────────────────────────────────
n_sizes  <- c(1000, 5000, 10000)
n_reps   <- 10          # 반복 횟수 per n
drugs    <- unique(results$drug)

# ── 3. Subsampling loop ───────────────────────────────────────────────────────
stability_rows <- list()

for (n_target in n_sizes) {
  for (rep in seq_len(n_reps)) {
    
    # 약물별 n_target명 샘플링
    sample_df <- results %>%
      group_by(drug) %>%
      group_modify(~ slice_sample(.x, n = min(n_target, nrow(.x)),
                                  replace = FALSE)) %>%
      ungroup()
    
    # 핵심 지표 계산
    metrics <- sample_df %>%
      group_by(drug) %>%
      summarise(
        med_t50      = median(t_50pct),
        med_AUC      = median(AUC_deficit_7d),
        med_dCdt     = median(abs(dCdt_rel_24h)),
        .groups      = "drop"
      ) %>%
      mutate(n_size = n_target, rep = rep)
    
    stability_rows[[length(stability_rows) + 1]] <- metrics
  }
}

stab_df <- bind_rows(stability_rows)

# ── 4. CV 계산 (핵심 수렴 지표) ───────────────────────────────────────────────
# CV = SD / mean × 100% — 약물별, n별
cv_table <- stab_df %>%
  group_by(drug, n_size) %>%
  summarise(
    CV_t50_pct  = round(sd(med_t50)  / mean(med_t50)  * 100, 2),
    CV_AUC_pct  = round(sd(med_AUC)  / mean(med_AUC)  * 100, 2),
    CV_dCdt_pct = round(sd(med_dCdt) / mean(med_dCdt) * 100, 2),
    .groups     = "drop"
  ) %>%
  arrange(drug, n_size)

cat("\n══════════════════════════════════════════════════════════════\n")
cat(" STABILITY CHECK — CV of median estimates across 10 replicates\n")
cat("══════════════════════════════════════════════════════════════\n\n")
print(cv_table, n = Inf)

# ── 5. Pass/Fail 판정 (CV < 5% = PASS) ───────────────────────────────────────
max_cv <- cv_table %>%
  filter(n_size == 10000) %>%
  summarise(max_cv = max(CV_t50_pct, CV_AUC_pct, CV_dCdt_pct)) %>%
  pull(max_cv)

cat(sprintf("\n── n=10,000 최대 CV: %.2f%%", max_cv))
if (max_cv < 5) {
  cat(" → ✅ PASS (< 5%): n=10,000 충분히 안정적\n")
} else {
  cat(" → ⚠ CONCERN (≥ 5%): n 증가 검토 필요\n")
}

# ── 6. Save CSV ───────────────────────────────────────────────────────────────
csv_path <- file.path(TABLE_DIR, "stability_check.csv")
write_csv(cv_table, csv_path)
cat(sprintf("\n✓ Stability table saved: %s\n", csv_path))

# ── 7. Plot: median t50% convergence by drug ──────────────────────────────────
drug_colors <- c(
  "Fluoxetine"     = "#2D6A4F",
  "Sertraline"     = "#2166AC",
  "Paroxetine"     = "#D95F02",
  "Venlafaxine_IR" = "#7B2D8B",
  "Venlafaxine_XR" = "#264653"
)

# n=10,000 전체 median (참조값)
ref_df <- results %>%
  group_by(drug) %>%
  summarise(ref_t50 = median(t_50pct), .groups = "drop")

p_stab <- stab_df %>%
  left_join(ref_df, by = "drug") %>%
  mutate(
    n_label = factor(paste0("n = ", scales::comma(n_size)),
                     levels = paste0("n = ", scales::comma(n_sizes))),
    drug    = factor(drug, levels = names(drug_colors))
  ) %>%
  ggplot(aes(x = n_label, y = med_t50, colour = drug, group = drug)) +
  # 반복 샘플 점
  geom_jitter(width = 0.08, alpha = 0.5, size = 1.5) +
  # 평균선
  stat_summary(fun = mean, geom = "line",
               linewidth = 0.8, linetype = "solid") +
  stat_summary(fun = mean, geom = "point",
               size = 3, shape = 21, fill = "white", stroke = 1.2) +
  # n=10,000 참조 수평선
  geom_hline(data = ref_df,
             aes(yintercept = ref_t50, colour = drug),
             linetype = "dashed", linewidth = 0.4, alpha = 0.6) +
  scale_colour_manual(values = drug_colors, name = "Drug") +
  scale_y_log10(
    breaks = c(5, 10, 24, 72, 168, 360, 720),
    labels = c("5h","10h","1d","3d","7d","15d","30d")
  ) +
  facet_wrap(~ drug, scales = "free_y", ncol = 5) +
  labs(
    x       = "Sample size",
    y       = expression(italic(t)[50]*"% C"[0]*"  median  (log scale)"),
    caption = paste0(
      "Points = median from each of 10 random subsamples. ",
      "Line = mean across replicates. ",
      "Dashed = full n=10,000 reference.\n",
      "CV < 5% at n=10,000 indicates sufficient convergence."
    )
  ) +
  theme_bw(base_size = 10) +
  theme(
    legend.position  = "none",
    strip.text       = element_text(face = "bold", size = 9),
    panel.grid.minor = element_blank(),
    plot.caption     = element_text(size = 7.5, colour = "grey40",
                                    hjust = 0, lineheight = 1.3)
  )

plot_path <- file.path(FIG_DIR, "stability_check_plot.png")
ggsave(plot_path, p_stab, width = 13, height = 4, dpi = 300)
cat(sprintf("✓ Stability plot saved: %s\n", plot_path))

cat("\n▶ Stability check 완료 — Phase 5 진행 가능\n\n")