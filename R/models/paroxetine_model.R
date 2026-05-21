# ============================================================
# Phase 3: Paroxetine 1-Compartment PK Model
# Author : Nari Kim
# Date   : 2025-05
# Model  : 1-compartment, first-order absorption
# Params : Ka=0.580/h, Vd=609L, CL=20.1 L/h, t½=21h
# Source : FDA label (Paxil); DrugBank Vd; Li 2022 Tmax
# Note   : IR only (SR 제거 — 중단 후 소실 패턴 동일)
# ============================================================

library(rxode2)
library(ggplot2)
library(dplyr)

# ── 1. 모델 정의 ───────────────────────────────────────────
paroxetine_model <- rxode2({
  C           = A_c / Vd
  d/dt(A_gut) = -Ka * A_gut
  d/dt(A_c)   =  Ka * A_gut - (CL / Vd) * A_c
})

# ── 2. 파라미터 (v4 검증값) ────────────────────────────────
# CL = 0.693 × 609 / 21 = 20.1 L/h (FDA label t½ 역산)
params <- c(
  Ka = 0.580,    # 1/h  | Tmax 5.2h 역산
  Vd = 609,      # L    | 8.7 L/kg × 70kg
  CL = 20.1      # L/h  | FDA label t½=21h 역산
)

# ── 3. 시간 설정 ───────────────────────────────────────────
dose      <- 30       # mg (표준 치료 용량)
t_ss_days <- 14       # steady-state 복용 기간 (t½=21h, SS~4일이나 14일로 안전하게)
t_ss      <- t_ss_days * 24    # 336h

disc_abrupt  <- t_ss
disc_taper2w <- t_ss + 168     # +7일
disc_taper4w <- t_ss + 504     # +21일

t_post  <- 14 * 24             # 중단 후 14일 (t½=21h → 5×21=105h면 소실, 14일로 여유)
t_total <- disc_taper4w + t_post
obs_seq <- seq(0, t_total, by = 0.25)

# ── 4. 정보 출력 ───────────────────────────────────────────
t_half <- 0.693 * 609 / 20.1
cat(sprintf("t½ = %.1f h (%.1f days)\n", t_half, t_half / 24))
cat(sprintf("Tmax = %.1f h (목표 5.2h)\n",
            log(0.580 / (20.1/609)) / (0.580 - 20.1/609)))

# ── 5. 폴더 생성 ───────────────────────────────────────────
dir.create(
  "~/nari-research/pkpd-antidepressant-sim/figures",
  showWarnings = FALSE, recursive = TRUE
)

# ── 6. 이벤트 테이블 ───────────────────────────────────────
ev_abrupt <- et(amt = dose, cmt = "A_gut",
                ii = 24, addl = t_ss_days - 1, time = 0) %>%
  et(obs_seq)

ev_taper2w <- et(amt = dose, cmt = "A_gut",
                 ii = 24, addl = t_ss_days - 1, time = 0) %>%
  et(amt = dose * 0.5, cmt = "A_gut", ii = 24, addl = 6, time = t_ss) %>%
  et(obs_seq)

ev_taper4w <- et(amt = dose, cmt = "A_gut",
                 ii = 24, addl = t_ss_days - 1, time = 0) %>%
  et(amt = dose * 0.75, cmt = "A_gut", ii = 24, addl = 6, time = t_ss) %>%
  et(amt = dose * 0.50, cmt = "A_gut", ii = 24, addl = 6, time = t_ss + 168) %>%
  et(amt = dose * 0.25, cmt = "A_gut", ii = 24, addl = 6, time = t_ss + 336) %>%
  et(obs_seq)

# ── 7. 시뮬레이션 ──────────────────────────────────────────
inits <- c(A_gut = 0, A_c = 0)

sim_abrupt  <- rxSolve(paroxetine_model, params, ev_abrupt,  inits)
sim_taper2w <- rxSolve(paroxetine_model, params, ev_taper2w, inits)
sim_taper4w <- rxSolve(paroxetine_model, params, ev_taper4w, inits)

# ── 8. 중단 후 구간만 추출 ─────────────────────────────────
disc_times <- c(
  "Abrupt"      = disc_abrupt,
  "Tapering 2w" = disc_taper2w,
  "Tapering 4w" = disc_taper4w
)

sim_all <- bind_rows(
  as.data.frame(sim_abrupt) %>%
    filter(time >= disc_abrupt) %>%
    mutate(scenario  = "Abrupt",
           time_rel  = (time - disc_abrupt) / 24),
  as.data.frame(sim_taper2w) %>%
    filter(time >= disc_taper2w) %>%
    mutate(scenario  = "Tapering 2w",
           time_rel  = (time - disc_taper2w) / 24),
  as.data.frame(sim_taper4w) %>%
    filter(time >= disc_taper4w) %>%
    mutate(scenario  = "Tapering 4w",
           time_rel  = (time - disc_taper4w) / 24)
) %>%
  mutate(scenario = factor(scenario,
                           levels = c("Abrupt", "Tapering 2w", "Tapering 4w")))

# ── 9. Css 및 50% threshold 계산 ──────────────────────────
get_css <- function(df, scen) {
  df %>%
    filter(scenario == scen, abs(time_rel) < 0.2) %>%
    arrange(abs(time_rel)) %>%
    slice(1) %>%
    pull(C)
}

Css_abrupt  <- get_css(sim_all, "Abrupt")
Css_taper2w <- get_css(sim_all, "Tapering 2w")
Css_taper4w <- get_css(sim_all, "Tapering 4w")

cat(sprintf("Css — Abrupt: %.5f | 2w: %.5f | 4w: %.5f mg/L\n",
            Css_abrupt, Css_taper2w, Css_taper4w))

# ── 10. dC/dt 계산 ─────────────────────────────────────────
calc_dcdt <- function(df, scen) {
  df %>%
    filter(scenario == scen, time_rel >= 0, time_rel <= 1) %>%
    arrange(time_rel) %>%
    mutate(dCdt = c(NA, diff(C) / diff(time_rel * 24))) %>%
    summarise(mean_dCdt = mean(dCdt, na.rm = TRUE)) %>%
    pull(mean_dCdt)
}

dCdt_a  <- calc_dcdt(sim_all, "Abrupt")
dCdt_2w <- calc_dcdt(sim_all, "Tapering 2w")
dCdt_4w <- calc_dcdt(sim_all, "Tapering 4w")

cat(sprintf("dC/dt (첫 24h) — Abrupt:%.5f | 2w:%.5f | 4w:%.5f mg/L/h\n",
            dCdt_a, dCdt_2w, dCdt_4w))

# ── 11. 시각화 ─────────────────────────────────────────────
col_scen <- c("Abrupt"="#E24B4A","Tapering 2w"="#378ADD","Tapering 4w"="#1D9E75")
lty_scen <- c("Abrupt"="solid","Tapering 2w"="dashed","Tapering 4w"="dotdash")

p <- ggplot(sim_all,
            aes(x = time_rel, y = C,
                color = scenario, linetype = scenario)) +
  geom_line(linewidth = 1.0) +
  
  # 시나리오별 50% 기준선
  geom_hline(yintercept = Css_abrupt  * 0.5,
             color="#E24B4A", linetype="dotted", linewidth=0.5) +
  geom_hline(yintercept = Css_taper2w * 0.5,
             color="#378ADD", linetype="dotted", linewidth=0.5) +
  geom_hline(yintercept = Css_taper4w * 0.5,
             color="#1D9E75", linetype="dotted", linewidth=0.5) +
  
  annotate("text", x = 13.8, y = Css_abrupt  * 0.5 * 1.25,
           label = "50% Css (A)",   size = 2.8, color="#E24B4A", hjust=1) +
  annotate("text", x = 13.8, y = Css_taper2w * 0.5 * 1.25,
           label = "50% Css (2w)", size = 2.8, color="#378ADD", hjust=1) +
  annotate("text", x = 13.8, y = Css_taper4w * 0.5 * 1.25,
           label = "50% Css (4w)", size = 2.8, color="#1D9E75", hjust=1) +
  
  scale_color_manual(values = col_scen) +
  scale_linetype_manual(values = lty_scen) +
  coord_cartesian(xlim = c(0, 14)) +
  labs(
    title    = "Paroxetine — Plasma concentration after discontinuation",
    subtitle = sprintf(
      "Ka=0.58/h, Vd=609L, CL=20.1 L/h | t\u00bd=%.1fh (%.1f days) | FDA label",
      t_half, t_half / 24),
    x        = "Days after discontinuation",
    y        = "Plasma concentration (mg/L)",
    color    = "Scenario", linetype = "Scenario",
    caption  = paste0(
      "Shown from time of discontinuation after 14-day dosing period.\n",
      "Dotted lines = 50% of each scenario's Css at time of discontinuation.")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position  = "bottom",
    plot.title       = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.key.width = unit(1.8, "cm")
  )

print(p)

ggsave(
  "~/nari-research/pkpd-antidepressant-sim/figures/paroxetine_discontinuation.png",
  plot = p, width = 10, height = 6, dpi = 300
)
cat("✅ Paroxetine 모델 완료.\n")