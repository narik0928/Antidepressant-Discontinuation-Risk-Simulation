# ============================================================
# Phase 3: Venlafaxine IR vs XR PK Model
# Author : Nari Kim
# v5     : 범례 왼쪽 패널에만 표시, 오른쪽 없앰
#          guide_area() 제거 (버전 호환성 문제)
# ============================================================

library(rxode2)
library(ggplot2)
library(dplyr)

if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
library(patchwork)

# ── 1. 모델 ────────────────────────────────────────────────
venlafaxine_model <- rxode2({
  C_ven   = A_ven / Vd_ven
  C_odv   = A_odv / Vd_odv
  C_total = C_ven + C_odv
  d/dt(A_gut) = -Ka * A_gut
  d/dt(A_ven) =  Ka * A_gut - (CL_ven / Vd_ven) * A_ven
  d/dt(A_odv) =  Fm * CL_ven * C_ven - (CL_odv / Vd_odv) * A_odv
})

# ── 2. 파라미터 ────────────────────────────────────────────
params_IR <- c(Ka=1.50, Vd_ven=525, CL_ven=91, Vd_odv=399, CL_odv=28, Fm=0.80)
params_XR <- c(Ka=0.25, Vd_ven=525, CL_ven=91, Vd_odv=399, CL_odv=28, Fm=0.80)

# ── 3. 시간 설정 ───────────────────────────────────────────
dose      <- 75
t_ss_days <- 7
t_ss      <- t_ss_days * 24
disc_abrupt  <- t_ss
disc_taper2w <- t_ss + 168
disc_taper4w <- t_ss + 504
t_total   <- disc_taper4w + 7 * 24
obs_seq   <- seq(0, t_total, by = 0.25)

dir.create(
  "~/nari-research/pkpd-antidepressant-sim/figures",
  showWarnings = FALSE, recursive = TRUE
)

# ── 4. 이벤트 테이블 ───────────────────────────────────────
make_events <- function(scenario) {
  base <- et(amt = dose, cmt = "A_gut", ii = 24, addl = t_ss_days - 1, time = 0)
  if (scenario == "Abrupt") return(base %>% et(obs_seq))
  if (scenario == "Tapering 2w")
    return(base %>%
             et(amt = dose * 0.5,  cmt = "A_gut", ii = 24, addl = 6, time = t_ss) %>%
             et(obs_seq))
  if (scenario == "Tapering 4w")
    return(base %>%
             et(amt = dose * 0.75, cmt = "A_gut", ii = 24, addl = 6, time = t_ss) %>%
             et(amt = dose * 0.50, cmt = "A_gut", ii = 24, addl = 6, time = t_ss + 168) %>%
             et(amt = dose * 0.25, cmt = "A_gut", ii = 24, addl = 6, time = t_ss + 336) %>%
             et(obs_seq))
}

# ── 5. 시뮬레이션 ──────────────────────────────────────────
inits      <- c(A_gut = 0, A_ven = 0, A_odv = 0)
scenarios  <- c("Abrupt", "Tapering 2w", "Tapering 4w")
disc_times <- c(Abrupt=disc_abrupt, "Tapering 2w"=disc_taper2w,
                "Tapering 4w"=disc_taper4w)

run_sim <- function(params, label) {
  bind_rows(lapply(scenarios, function(scen) {
    sim <- rxSolve(venlafaxine_model, params, make_events(scen), inits)
    as.data.frame(sim) %>%
      mutate(scenario    = scen,
             formulation = label,
             time_rel    = (time - disc_times[scen]) / 24)
  }))
}

cat("IR 시뮬레이션...\n"); sim_IR <- run_sim(params_IR, "IR")
cat("XR 시뮬레이션...\n"); sim_XR <- run_sim(params_XR, "XR")

sim_all <- bind_rows(sim_IR, sim_XR) %>%
  mutate(
    scenario    = factor(scenario,
                         levels = c("Abrupt", "Tapering 2w", "Tapering 4w")),
    formulation = factor(formulation, levels = c("IR", "XR"))
  )

# ── 6. Css / 50% 기준선 ────────────────────────────────────
get_css <- function(df, scen, form, col) {
  df %>%
    filter(scenario == scen, formulation == form, abs(time_rel) < 0.2) %>%
    arrange(abs(time_rel)) %>% slice(1) %>% pull({{ col }})
}

get_thresholds <- function(form_label) {
  data.frame(
    yint  = c(
      get_css(sim_all, "Abrupt",      form_label, C_total) * 0.5,
      get_css(sim_all, "Tapering 2w", form_label, C_total) * 0.5,
      get_css(sim_all, "Tapering 4w", form_label, C_total) * 0.5
    ),
    color = c("#E24B4A", "#378ADD", "#1D9E75")
  )
}

thr_IR <- get_thresholds("IR")
thr_XR <- get_thresholds("XR")

cat(sprintf("Css total — IR Abrupt: %.4f | XR Abrupt: %.4f mg/L\n",
            thr_IR$yint[1] * 2, thr_XR$yint[1] * 2))

# ── 7. 색상/선 설정 ────────────────────────────────────────
col_scen <- c("Abrupt"="#E24B4A","Tapering 2w"="#378ADD","Tapering 4w"="#1D9E75")
lty_scen <- c("Abrupt"="solid","Tapering 2w"="dashed","Tapering 4w"="dotdash")

# ── 8. 패널 생성 함수 ──────────────────────────────────────
make_panel <- function(form_label, panel_title, thr_df, show_legend) {
  df <- sim_all %>%
    filter(formulation == form_label, time_rel >= 0, time_rel <= 5)
  
  p <- ggplot(df, aes(x = time_rel, y = C_total,
                      color = scenario, linetype = scenario)) +
    geom_line(linewidth = 0.9) +
    geom_hline(yintercept = thr_df$yint[1], color = thr_df$color[1],
               linetype = "dotted", linewidth = 0.5) +
    geom_hline(yintercept = thr_df$yint[2], color = thr_df$color[2],
               linetype = "dotted", linewidth = 0.5) +
    geom_hline(yintercept = thr_df$yint[3], color = thr_df$color[3],
               linetype = "dotted", linewidth = 0.5) +
    scale_color_manual(values = col_scen, name = NULL) +
    scale_linetype_manual(values = lty_scen, name = NULL) +
    labs(title = panel_title,
         x = "Days after discontinuation",
         y = "Total active concentration (mg/L)") +
    theme_minimal(base_size = 11) +
    theme(
      plot.title       = element_text(face = "bold", size = 11),
      panel.grid.minor = element_blank(),
      plot.margin      = margin(8, 12, 4, 20)
    )
  
  if (show_legend) {
    # 왼쪽 패널: 범례 하단 표시
    p <- p + theme(
      legend.position  = "bottom",
      legend.key.width = unit(1.8, "cm"),
      legend.text      = element_text(size = 10),
      legend.title     = element_text(size = 10)
    )
  } else {
    # 오른쪽 패널: 범례 없음
    p <- p + theme(legend.position = "none")
  }
  p
}

pA <- make_panel("IR", "A. IR  (Ka = 1.50/h, Tmax ~ 2h)",  thr_IR, show_legend = TRUE)
pB <- make_panel("XR", "B. XR  (Ka = 0.25/h, Tmax ~ 5.5h)", thr_XR, show_legend = FALSE)

# ── 9. 조합 및 저장 ────────────────────────────────────────
fig <- (pA | pB) +
  plot_annotation(
    title    = "Venlafaxine IR vs XR — VEN+ODV total active moiety after discontinuation",
    subtitle = "Vd_ven=525L, CL_ven=91 L/h | Vd_odv=399L, CL_odv=28 L/h | Fm=0.80 | FDA label (Effexor XR 2012)",
    caption  = "Shown from time of discontinuation after 7-day dosing period.\nDotted lines = 50% of each scenario's Css at time of discontinuation.",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9),
      plot.caption  = element_text(size = 9)
    )
  )

print(fig)

ggsave(
  "~/nari-research/pkpd-antidepressant-sim/figures/venlafaxine_post_discontinuation.png",
  fig, width = 12, height = 6.5, dpi = 300
)
cat("✅ Venlafaxine 모델 v5 완료.\n")