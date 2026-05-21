# ============================================================
# Phase 3: Fluoxetine 1-Compartment + Norfluoxetine Model
# Author : Nari Kim
# v2     : 패널 제목 단축, 범례 분리 (왼쪽만)
# ============================================================

library(rxode2)
library(ggplot2)
library(dplyr)

if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
library(patchwork)

# ── 1. 모델 정의 ───────────────────────────────────────────
fluoxetine_model <- rxode2({
  C_flx   = A_flx / Vd
  C_nfx   = A_nfx / Vd_nfx
  C_total = C_flx + C_nfx
  d/dt(A_gut) = -Ka * A_gut
  d/dt(A_flx) =  Ka * A_gut - (CL / Vd) * A_flx
  d/dt(A_nfx) =  Fm * CL * C_flx - (CL_nfx / Vd_nfx) * A_nfx
})

# ── 2. 파라미터 ────────────────────────────────────────────
params <- c(
  Ka      = 0.72,
  Vd      = 2310,
  CL      = 10.5,
  Vd_nfx  = 2310,
  CL_nfx  = 7.18,
  Fm      = 0.72
)

# ── 3. 시간 설정 ───────────────────────────────────────────
dose      <- 20
t_ss_days <- 60
t_ss      <- t_ss_days * 24
disc_abrupt  <- t_ss
disc_taper2w <- t_ss + 168
disc_taper4w <- t_ss + 504
t_post    <- 60 * 24
t_total   <- disc_taper4w + t_post
obs_seq   <- seq(0, t_total, by = 1)

t_half_flx <- 0.693 * 2310 / 10.5
t_half_nfx <- 0.693 * 2310 / 7.18
cat(sprintf("FLX t½=%.1fd | NFX t½=%.1fd\n",
            t_half_flx/24, t_half_nfx/24))

dir.create(
  "~/nari-research/pkpd-antidepressant-sim/figures",
  showWarnings = FALSE, recursive = TRUE
)

# ── 4. 이벤트 테이블 ───────────────────────────────────────
inits <- c(A_gut = 0, A_flx = 0, A_nfx = 0)

ev_abrupt <- et(amt = dose, cmt = "A_gut",
                ii = 24, addl = t_ss_days - 1, time = 0) %>%
  et(obs_seq)

ev_taper2w <- et(amt = dose, cmt = "A_gut",
                 ii = 24, addl = t_ss_days - 1, time = 0) %>%
  et(amt = dose * 0.5,  cmt = "A_gut", ii = 24, addl = 6, time = t_ss) %>%
  et(obs_seq)

ev_taper4w <- et(amt = dose, cmt = "A_gut",
                 ii = 24, addl = t_ss_days - 1, time = 0) %>%
  et(amt = dose * 0.75, cmt = "A_gut", ii = 24, addl = 6, time = t_ss) %>%
  et(amt = dose * 0.50, cmt = "A_gut", ii = 24, addl = 6, time = t_ss + 168) %>%
  et(amt = dose * 0.25, cmt = "A_gut", ii = 24, addl = 6, time = t_ss + 336) %>%
  et(obs_seq)

# ── 5. 시뮬레이션 ──────────────────────────────────────────
cat("시뮬레이션 실행 중...\n")
sim_abrupt  <- rxSolve(fluoxetine_model, params, ev_abrupt,  inits)
sim_taper2w <- rxSolve(fluoxetine_model, params, ev_taper2w, inits)
sim_taper4w <- rxSolve(fluoxetine_model, params, ev_taper4w, inits)
cat("완료!\n")

# ── 6. 중단 후 구간 추출 ───────────────────────────────────
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

# ── 7. Css 계산 ────────────────────────────────────────────
get_css <- function(df, scen, col) {
  df %>%
    filter(scenario == scen, abs(time_rel) < 0.5) %>%
    arrange(abs(time_rel)) %>%
    slice(1) %>%
    pull({{ col }})
}

Css_flx <- get_css(sim_all, "Abrupt", C_flx)
Css_nfx <- get_css(sim_all, "Abrupt", C_nfx)
Css_a   <- get_css(sim_all, "Abrupt",      C_total)
Css_2w  <- get_css(sim_all, "Tapering 2w", C_total)
Css_4w  <- get_css(sim_all, "Tapering 4w", C_total)

cat(sprintf("Css FLX:%.5f NFX:%.5f Total:%.5f mg/L\n",
            Css_flx, Css_nfx, Css_a))

# ── 8. 테마 설정 ───────────────────────────────────────────
col_scen <- c("Abrupt"="#E24B4A","Tapering 2w"="#378ADD","Tapering 4w"="#1D9E75")
lty_scen <- c("Abrupt"="solid","Tapering 2w"="dashed","Tapering 4w"="dotdash")

theme_left <- theme_minimal(base_size = 11) +
  theme(legend.position  = "bottom",
        legend.key.width = unit(1.6, "cm"),
        legend.text      = element_text(size = 10),
        plot.title       = element_text(face = "bold", size = 11),
        panel.grid.minor = element_blank(),
        plot.margin      = margin(8, 8, 4, 80))

theme_right <- theme_minimal(base_size = 11) +
  theme(legend.position  = "none",
        plot.title       = element_text(face = "bold", size = 11),
        panel.grid.minor = element_blank(),
        plot.margin      = margin(8, 8, 4, 80))

# ── 9. Panel A: FLX vs NFX (Abrupt, 왼쪽 — 범례 있음) ──────
df_abrupt <- sim_all %>% filter(scenario == "Abrupt")

pA <- ggplot(df_abrupt, aes(x = time_rel)) +
  geom_line(aes(y = C_flx,   color = "Fluoxetine"),
            linewidth = 1.0, linetype = "solid") +
  geom_line(aes(y = C_nfx,   color = "Norfluoxetine"),
            linewidth = 1.0, linetype = "dashed") +
  geom_line(aes(y = C_total, color = "FLX+NFX total"),
            linewidth = 0.7, linetype = "dotdash") +
  scale_color_manual(
    values = c("Fluoxetine"    = "#E24B4A",
               "Norfluoxetine" = "#378ADD",
               "FLX+NFX total" = "#888780"),
    name = NULL
  ) +
  coord_cartesian(xlim = c(0, 60)) +
  labs(title = "A. FLX vs NFX — Abrupt",
       x     = "Days after discontinuation",
       y     = "Plasma concentration (mg/L)") +
  theme_left

# ── 10. Panel B: 시나리오 비교 (오른쪽 — 범례 없음) ────────
pB <- ggplot(sim_all, aes(x = time_rel, y = C_total,
                          color = scenario, linetype = scenario)) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = Css_a  * 0.5,
             color="#E24B4A", linetype="dotted", linewidth=0.5) +
  geom_hline(yintercept = Css_2w * 0.5,
             color="#378ADD", linetype="dotted", linewidth=0.5) +
  geom_hline(yintercept = Css_4w * 0.5,
             color="#1D9E75", linetype="dotted", linewidth=0.5) +
  annotate("text", x=15, y=Css_a *0.5*1.15,
           label="50% Css (A)",  size=2.8, color="#E24B4A", hjust=0) +
  annotate("text", x=15, y=Css_2w*0.5*1.15,
           label="50% Css (2w)", size=2.8, color="#378ADD", hjust=0) +
  annotate("text", x=15, y=Css_4w*0.5*1.15,
           label="50% Css (4w)", size=2.8, color="#1D9E75", hjust=0) +
  scale_color_manual(values = col_scen) +
  scale_linetype_manual(values = lty_scen) +
  coord_cartesian(xlim = c(0, 60)) +
  labs(title = "B. FLX+NFX — Scenarios",
       x     = "Days after discontinuation",
       y     = "Total active concentration (mg/L)") +
  theme_right

# ── 11. 조합 및 저장 ───────────────────────────────────────
fig <- (pA | pB) +
  plot_annotation(
    title    = "Fluoxetine — Plasma concentration after discontinuation",
    subtitle = sprintf(
      "Ka=0.72/h, Vd=2310L, CL=10.5 L/h | FLX t\u00bd=%.1fd | NFX t\u00bd=%.1fd | FDA label",
      t_half_flx/24, t_half_nfx/24),
    caption  = "Shown from time of discontinuation after 60-day dosing period.\nDotted lines = 50% of each scenario's Css (total FLX+NFX).",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9),
      plot.caption  = element_text(size = 9)
    )
  )

print(fig)

ggsave(
  "~/nari-research/pkpd-antidepressant-sim/figures/fluoxetine_discontinuation.png",
  fig, width = 13, height = 6.5, dpi = 300
)
cat("✅ Fluoxetine 모델 v2 완료.\n")