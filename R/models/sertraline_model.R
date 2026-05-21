# ============================================================
# Phase 3: Sertraline 1-Compartment PK Model
# Author : Nari Kim
# v5     : 시나리오별 Css 50% 기준선 색깔별 표시
# ============================================================

library(rxode2)
library(ggplot2)
library(dplyr)

# ── 1. 모델 정의 ───────────────────────────────────────────
sertraline_model <- rxode2({
  C           = A_c / Vd
  d/dt(A_gut) = -Ka * A_gut
  d/dt(A_c)   =  Ka * A_gut - (CL / Vd) * A_c
})

# ── 2. 파라미터 ────────────────────────────────────────────
params <- c(Ka = 0.50, Vd = 1400, CL = 37.3)

dose      <- 50
t_ss_days <- 14
t_ss      <- t_ss_days * 24
t_post    <- 30 * 24
t_total   <- t_ss + t_post
obs_seq   <- seq(0, t_total, by = 0.25)

t_half <- 0.693 * 1400 / 37.3
cat(sprintf("t½ = %.1f h (%.1f days)\n", t_half, t_half / 24))

# ── 3. 폴더 생성 ───────────────────────────────────────────
dir.create(
  "~/nari-research/pkpd-antidepressant-sim/figures",
  showWarnings = FALSE, recursive = TRUE
)

# ── 4. 이벤트 테이블 ───────────────────────────────────────
disc_abrupt  <- t_ss
disc_taper2w <- t_ss + 168
disc_taper4w <- t_ss + 504

ev_abrupt <- et(amt = dose, cmt = "A_gut",
                ii = 24, addl = t_ss_days - 1, time = 0) %>%
  et(obs_seq)

ev_taper2w <- et(amt = dose, cmt = "A_gut",
                 ii = 24, addl = t_ss_days - 1, time = 0) %>%
  et(amt = 25, cmt = "A_gut", ii = 24, addl = 6, time = t_ss) %>%
  et(obs_seq)

ev_taper4w <- et(amt = dose, cmt = "A_gut",
                 ii = 24, addl = t_ss_days - 1, time = 0) %>%
  et(amt = 37.5, cmt = "A_gut", ii = 24, addl = 6, time = t_ss) %>%
  et(amt = 25,   cmt = "A_gut", ii = 24, addl = 6, time = t_ss + 168) %>%
  et(amt = 12.5, cmt = "A_gut", ii = 24, addl = 6, time = t_ss + 336) %>%
  et(obs_seq)

# ── 5. 시뮬레이션 ──────────────────────────────────────────
sim_abrupt  <- rxSolve(sertraline_model, params, ev_abrupt,  c(A_gut=0, A_c=0))
sim_taper2w <- rxSolve(sertraline_model, params, ev_taper2w, c(A_gut=0, A_c=0))
sim_taper4w <- rxSolve(sertraline_model, params, ev_taper4w, c(A_gut=0, A_c=0))

# ── 6. 중단 후 구간만 추출 ─────────────────────────────────
sim_all <- bind_rows(
  as.data.frame(sim_abrupt) %>%
    filter(time >= disc_abrupt) %>%
    mutate(scenario = "Abrupt",
           time_rel = (time - disc_abrupt) / 24),
  as.data.frame(sim_taper2w) %>%
    filter(time >= disc_taper2w) %>%
    mutate(scenario = "Tapering 2w",
           time_rel = (time - disc_taper2w) / 24),
  as.data.frame(sim_taper4w) %>%
    filter(time >= disc_taper4w) %>%
    mutate(scenario = "Tapering 4w",
           time_rel = (time - disc_taper4w) / 24)
) %>%
  mutate(scenario = factor(scenario,
                           levels = c("Abrupt", "Tapering 2w", "Tapering 4w")))

# ── 7. 시나리오별 Css 및 50% threshold 계산 ───────────────
get_css <- function(sim_df, disc_time) {
  as.data.frame(sim_df) %>%
    filter(time >= disc_time - 1, time <= disc_time) %>%
    arrange(time) %>%
    slice(1) %>%
    pull(C)
}

Css_abrupt  <- get_css(sim_abrupt,  disc_abrupt)
Css_taper2w <- get_css(sim_taper2w, disc_taper2w)
Css_taper4w <- get_css(sim_taper4w, disc_taper4w)

thr_abrupt  <- Css_abrupt  * 0.50
thr_taper2w <- Css_taper2w * 0.50
thr_taper4w <- Css_taper4w * 0.50

cat(sprintf("Css — Abrupt: %.5f | 2w: %.5f | 4w: %.5f mg/L\n",
            Css_abrupt, Css_taper2w, Css_taper4w))
cat(sprintf("50%% threshold — Abrupt: %.5f | 2w: %.5f | 4w: %.5f mg/L\n",
            thr_abrupt, thr_taper2w, thr_taper4w))

# ── 8. dC/dt 계산 ──────────────────────────────────────────
calc_dcdt <- function(df, scen) {
  df %>%
    filter(scenario == scen, time_rel >= 0, time_rel <= 1) %>%
    arrange(time_rel) %>%
    mutate(dCdt = c(NA, diff(C) / diff(time_rel * 24))) %>%
    summarise(mean_dCdt = mean(dCdt, na.rm = TRUE)) %>%
    pull(mean_dCdt)
}

dCdt_abrupt  <- calc_dcdt(sim_all, "Abrupt")
dCdt_taper2w <- calc_dcdt(sim_all, "Tapering 2w")
dCdt_taper4w <- calc_dcdt(sim_all, "Tapering 4w")

cat(sprintf("dC/dt (첫 24h) — Abrupt: %.5f | 2w: %.5f | 4w: %.5f mg/L/h\n",
            dCdt_abrupt, dCdt_taper2w, dCdt_taper4w))

# ── 9. 시각화 ──────────────────────────────────────────────
colors <- c(
  "Abrupt"      = "#E24B4A",
  "Tapering 2w" = "#378ADD",
  "Tapering 4w" = "#1D9E75"
)

p <- ggplot(sim_all,
            aes(x = time_rel, y = C,
                color = scenario, linetype = scenario)) +
  geom_line(linewidth = 1.0) +
  
  # 시나리오별 Css 50% 기준선 (각자 색깔로 표시)
  geom_hline(yintercept = thr_abrupt,
             color = "#E24B4A", linetype = "dotted", linewidth = 0.6) +
  geom_hline(yintercept = thr_taper2w,
             color = "#378ADD", linetype = "dotted", linewidth = 0.6) +
  geom_hline(yintercept = thr_taper4w,
             color = "#1D9E75", linetype = "dotted", linewidth = 0.6) +
  
  # 기준선 라벨
  annotate("text", x = 29.5, y = thr_abrupt  * 1.20,
           label = "50% Css (A)",   size = 2.8,
           color = "#E24B4A", hjust = 1) +
  annotate("text", x = 29.5, y = thr_taper2w * 1.20,
           label = "50% Css (2w)", size = 2.8,
           color = "#378ADD", hjust = 1) +
  annotate("text", x = 29.5, y = thr_taper4w * 1.20,
           label = "50% Css (4w)", size = 2.8,
           color = "#1D9E75", hjust = 1) +
  
  scale_color_manual(values = colors) +
  scale_linetype_manual(values = c(
    "Abrupt"      = "solid",
    "Tapering 2w" = "dashed",
    "Tapering 4w" = "dotdash"
  )) +
  labs(
    title    = "Sertraline — Plasma concentration after discontinuation",
    subtitle = sprintf(
      "Ka=0.50/h, Vd=1400L, CL=37.3 L/h | t\u00bd=%.1fh (%.1f days) | FDA label",
      t_half, t_half / 24),
    x        = "Days after discontinuation",
    y        = "Plasma concentration (mg/L)",
    color    = "Scenario",
    linetype = "Scenario",
    caption  = paste0(
      "Shown from time of discontinuation after 14-day dosing period.\n",
      "Dotted lines = 50% of each scenario's Css at time of discontinuation.")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position  = "bottom",
    plot.title       = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

print(p)

# ── 10. 저장 ───────────────────────────────────────────────
ggsave(
  filename = "~/nari-research/pkpd-antidepressant-sim/figures/sertraline_discontinuation.png",
  plot = p, width = 10, height = 6, dpi = 300
)
cat("✅ Sertraline 모델 v5 완료.\n")