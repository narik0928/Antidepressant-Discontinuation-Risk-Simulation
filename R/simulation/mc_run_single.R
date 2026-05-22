# =============================================================================
# mc_run_single.R
# Phase 4 — Monte Carlo Simulation Runner
#
# Purpose : mc_params.R 출력 파라미터로 환자별 PK 시뮬레이션 실행 →
#           PK discontinuation metrics 추출
#
# Inputs  : outputs/mc/mc_params_all_n<N>.csv
# Outputs : outputs/mc/mc_results_<drug>_n<N>.csv  (약물별)
#           outputs/mc/mc_results_all_n<N>.csv      (전체 합본)
#
# Metrics (per patient):
#   Css            : 중단 시점 steady-state 농도 (mg/L)
#   dCdt_rel_24h   : relative dC/dt 0–24h = (C24h - Css) / (24 × Css)  [/h]
#   dCdt_rel_72h   : relative dC/dt 0–72h = (C72h - Css) / (72 × Css)  [/h]
#   pct_drop_24h   : 24h 후 Css 대비 감소율 (%)
#   pct_drop_72h   : 72h 후 Css 대비 감소율 (%)
#   t_75pct        : Css 75% 도달 시간 (h)
#   t_50pct        : Css 50% 도달 시간 (h)
#   t_25pct        : Css 25% 도달 시간 (h)
#   t_10pct        : Css 10% 도달 시간 (h)
#   AUC_deficit_7d : 중단 후 7일간 누적 노출 손실 (mg/L×h)
#
# NOTE: Fluoxetine, Venlafaxine → Ctotal = Cp + Cm (unweighted sum, base-case)
#       Limitation: equal PD potency assumed for parent and active metabolite.
#
# 실행: source("R/simulation/mc_run_single.R")
# =============================================================================

suppressPackageStartupMessages({
  library(rxode2)
  library(tidyverse)
})

# ── 0. Setup ──────────────────────────────────────────────────────────────────
BASE       <- "~/nari-research/pkpd-antidepressant-sim"
OUT_DIR    <- file.path(BASE, "outputs", "mc")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

set.seed(2025)   # 고정 seed (mc_params.R과 동일)

N          <- 10000   # 총 환자 수 (테스트 시 100으로 변경)
BATCH_SIZE <- 1000    # 배치 크기 (메모리 관리)

source(file.path(BASE, "data", "pk_parameters", "pk_parameters.R"))

# ── 1. ODE Models ─────────────────────────────────────────────────────────────

ode_1cmt <- rxode2({
  d/dt(depot) = -Ka * depot
  d/dt(Cp)    =  Ka * depot / Vd - (CL / Vd) * Cp
})

ode_parent_met <- rxode2({
  d/dt(depot) = -Ka  * depot
  d/dt(Cp)    =  Ka  * depot / Vd_p - (CL_p / Vd_p) * Cp
  d/dt(Cm)    =  Fm  * CL_p  * Cp   / Vd_m - (CL_m  / Vd_m) * Cm
  Ctotal      =  Cp + Cm
})

# ── 2. Drug Configuration ─────────────────────────────────────────────────────

drug_config <- list(
  
  Sertraline = list(
    model    = "1cmt",
    Ka       = pk_params$sertraline$Ka,        # 0.50 /h
    F        = pk_params$sertraline$F,          # 1.0
    Fm       = NA,
    dose_mg  = 50,
    tau_h    = 24,
    n_doses  = 14,
    obs_h    = 30 * 24,    # 중단 후 관찰 (h)
    dt       = 1.0,        # 시간 해상도 (h)
    conc_col = "Cp"
  ),
  
  Fluoxetine = list(
    model    = "parent_met",
    Ka       = pk_params$fluoxetine$Ka,                   # 0.72 /h
    F        = pk_params$fluoxetine$F,                     # 0.70
    Fm       = pk_params$fluoxetine$norfluoxetine$Fm,      # 0.72
    dose_mg  = 20,
    tau_h    = 24,
    n_doses  = 60,
    obs_h    = 60 * 24,
    dt       = 1.0,
    conc_col = "Ctotal"
  ),
  
  Paroxetine = list(
    model    = "1cmt",
    Ka       = pk_params$paroxetine$Ka,        # 0.580 /h
    F        = pk_params$paroxetine$F,          # 0.50
    Fm       = NA,
    dose_mg  = 20,
    tau_h    = 24,
    n_doses  = 14,
    obs_h    = 14 * 24,
    dt       = 0.5,
    conc_col = "Cp"
  ),
  
  Venlafaxine_IR = list(
    model    = "parent_met",
    Ka       = pk_params$venlafaxine_IR$Ka,              # 1.50 /h
    F        = pk_params$venlafaxine_IR$F,               # 0.45
    Fm       = pk_params$venlafaxine_IR$ODV$Fm,          # 0.80
    dose_mg  = 25,     # 25mg TID = 75mg/day (XR 75mg QD와 total daily dose 매칭)
    tau_h    = 8,      # IR: TID
    n_doses  = 21,
    obs_h    = 7 * 24,
    dt       = 0.25,
    conc_col = "Ctotal"
  ),
  
  Venlafaxine_XR = list(
    model    = "parent_met",
    Ka       = pk_params$venlafaxine_XR$Ka,              # 0.25 /h
    F        = pk_params$venlafaxine_XR$F,               # 0.45
    Fm       = pk_params$venlafaxine_XR$ODV$Fm,          # 0.80
    dose_mg  = 75,
    tau_h    = 24,     # XR: QD
    n_doses  = 7,
    obs_h    = 7 * 24,
    dt       = 0.25,
    conc_col = "Ctotal"
  )
)

# ── 3. Helper Functions ───────────────────────────────────────────────────────

# 선형 보간으로 threshold 도달 시간 계산
# concs가 threshold 아래로 처음 내려가는 시점을 보간
find_threshold_time <- function(times, concs, threshold) {
  idx <- which(concs <= threshold)[1]
  if (is.na(idx) || idx == 1) return(NA_real_)   # 관찰 기간 내 미도달
  t1 <- times[idx - 1]; t2 <- times[idx]
  c1 <- concs[idx - 1]; c2 <- concs[idx]
  t1 + (threshold - c1) * (t2 - t1) / (c2 - c1)
}

# 사다리꼴 적분 (AUC)
trapz <- function(x, y) {
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}

# 환자 1명 PK metrics 추출
# .x : group_modify 내부 data (post_time, conc, Css 컬럼 포함)
extract_patient_metrics <- function(.x, Css, obs_h) {
  
  times <- .x$post_time
  concs <- .x$conc
  
  # 유효성 검사
  if (is.na(Css) || Css <= 0 || length(times) < 3) {
    return(tibble(
      Css = NA_real_, dCdt_rel_24h = NA_real_, dCdt_rel_72h = NA_real_,
      pct_drop_24h = NA_real_, pct_drop_72h = NA_real_,
      t_75pct = NA_real_, t_50pct = NA_real_,
      t_25pct = NA_real_, t_10pct = NA_real_,
      AUC_deficit_7d = NA_real_
    ))
  }
  
  # 특정 시점 농도 (선형 보간)
  C_24h <- approx(times, concs, xout = 24)$y
  C_72h <- approx(times, concs, xout = 72)$y
  
  # Relative dC/dt (Css 정규화) — 음수: 농도 감소
  dCdt_rel_24h <- if (!is.na(C_24h)) (C_24h - Css) / (24 * Css) else NA_real_
  dCdt_rel_72h <- if (!is.na(C_72h)) (C_72h - Css) / (72 * Css) else NA_real_
  
  # % Css drop — 양수: 감소
  pct_drop_24h <- if (!is.na(C_24h)) (Css - C_24h) / Css * 100 else NA_real_
  pct_drop_72h <- if (!is.na(C_72h)) (Css - C_72h) / Css * 100 else NA_real_
  
  # Threshold times (linear interpolation)
  # "operational PK thresholds" — NOT clinical withdrawal thresholds
  t_75pct <- find_threshold_time(times, concs, 0.75 * Css)
  t_50pct <- find_threshold_time(times, concs, 0.50 * Css)
  t_25pct <- find_threshold_time(times, concs, 0.25 * Css)
  t_10pct <- find_threshold_time(times, concs, 0.10 * Css)
  
  # AUC deficit (최대 7일 = 168h, 단 관찰 기간이 짧은 약물은 obs_h로 제한)
  # Venlafaxine obs_h=120h → min(168, 120) = 120h 기준으로 계산
  # AUC_expected: 중단 없이 Css 유지됐을 경우 기대 노출
  # AUC_actual:   실제 시뮬레이션 노출
  # Deficit > 0 → 노출 손실
  obs_max       <- min(168, obs_h)
  obs_7d        <- min(obs_max, max(times, na.rm = TRUE))
  post_7d    <- .x %>% filter(post_time <= obs_7d)
  AUC_actual    <- trapz(post_7d$post_time, post_7d$conc)
  AUC_expected  <- Css * obs_7d
  AUC_deficit_7d <- AUC_expected - AUC_actual
  
  tibble(
    Css            = round(Css, 6),
    dCdt_rel_24h   = round(dCdt_rel_24h, 6),
    dCdt_rel_72h   = round(dCdt_rel_72h, 6),
    pct_drop_24h   = round(pct_drop_24h, 3),
    pct_drop_72h   = round(pct_drop_72h, 3),
    t_75pct        = round(t_75pct, 3),
    t_50pct        = round(t_50pct, 3),
    t_25pct        = round(t_25pct, 3),
    t_10pct        = round(t_10pct, 3),
    AUC_deficit_7d = round(AUC_deficit_7d, 4)
  )
}

# ── 4. 배치 시뮬레이션 함수 ───────────────────────────────────────────────────

run_batch <- function(drug_name, batch_df, cfg) {
  
  t_disc <- cfg$n_doses * cfg$tau_h   # 중단 시점 (h)
  
  # 이벤트 테이블 (모든 환자 동일)
  ev <- et(amt  = cfg$F * cfg$dose_mg,
           ii   = cfg$tau_h,
           addl = cfg$n_doses - 1,
           time = 0) %>%
    et(seq(0, t_disc + cfg$obs_h, by = cfg$dt))
  
  # 환자별 파라미터 함수 정의
  if (cfg$model == "1cmt") {
    model <- ode_1cmt
    inits <- c(depot = 0, Cp = 0)
    get_params <- function(row) {
      c(Ka = cfg$Ka,
        Vd = row$Vd_sampled,
        CL = row$CL_sampled)
    }
  } else {
    model <- ode_parent_met
    inits <- c(depot = 0, Cp = 0, Cm = 0)
    get_params <- function(row) {
      c(Ka   = cfg$Ka,
        Vd_p  = row$Vd_sampled,
        CL_p  = row$CL_sampled,
        Fm    = cfg$Fm,
        Vd_m  = row$Vd_met_sampled,
        CL_m  = row$CL_met_sampled)
    }
  }
  
  # 환자별 단일 시뮬레이션 (lapply)
  # Cloud 실행 시 lapply → mclapply 교체로 병렬화 가능
  sim_list <- lapply(seq_len(nrow(batch_df)), function(i) {
    row   <- batch_df[i, , drop = FALSE]
    sim_i <- rxSolve(model, get_params(row), ev, inits = inits) %>%
      as.data.frame()
    sim_i$id <- row$id
    sim_i
  })
  
  sim <- bind_rows(sim_list) %>%
    rename_with(~ "conc", .cols = all_of(cfg$conc_col))
  
  # Css: 마지막 dosing interval 평균
  css_df <- sim %>%
    filter(time >= (cfg$n_doses - 1) * cfg$tau_h,
           time <= t_disc) %>%
    group_by(id) %>%
    summarise(Css = mean(conc, na.rm = TRUE), .groups = "drop")
  
  # 중단 후 구간 추출 + Css 합치기
  post_df <- sim %>%
    filter(time >= t_disc) %>%
    mutate(post_time = time - t_disc) %>%
    select(id, post_time, conc) %>%
    left_join(css_df, by = "id")
  
  # 환자별 metrics 추출
  metrics <- post_df %>%
    group_by(id) %>%
    group_modify(~ extract_patient_metrics(.x,
                                           unique(.x$Css),
                                           cfg$obs_h)) %>%
    ungroup()
  
  # patient info + metrics 합치기
  batch_df %>%
    select(id, drug, cyp2d6, CL_sampled, Vd_sampled) %>%
    left_join(metrics, by = "id")
}

# ── 5. 전체 실행 ──────────────────────────────────────────────────────────────

params_all <- read_csv(
  file.path(OUT_DIR, sprintf("mc_params_all_n%d.csv", N)),
  show_col_types = FALSE
)

results_list <- list()

for (drug_name in names(drug_config)) {
  
  cfg         <- drug_config[[drug_name]]
  drug_params <- params_all %>% filter(drug == drug_name)
  n_batches   <- ceiling(nrow(drug_params) / BATCH_SIZE)
  
  cat(sprintf("\n▶ %s: %d patients / %d batches\n",
              drug_name, nrow(drug_params), n_batches))
  
  t_start      <- proc.time()["elapsed"]
  drug_results <- vector("list", n_batches)
  
  for (b in seq_len(n_batches)) {
    idx_start       <- (b - 1) * BATCH_SIZE + 1
    idx_end         <- min(b * BATCH_SIZE, nrow(drug_params))
    batch_df        <- drug_params[idx_start:idx_end, ]
    drug_results[[b]] <- run_batch(drug_name, batch_df, cfg)
    cat(sprintf("  batch %02d/%02d done (n=%d)\n",
                b, n_batches, nrow(batch_df)))
  }
  
  drug_df  <- bind_rows(drug_results)
  elapsed  <- round(proc.time()["elapsed"] - t_start, 1)
  
  # 약물별 저장
  out_path <- file.path(
    OUT_DIR,
    sprintf("mc_results_%s_n%d.csv", tolower(gsub("_", "", drug_name)), N)
  )
  write_csv(drug_df, out_path)
  cat(sprintf("✓ Saved: %s  [%.1f sec]\n", out_path, elapsed))
  
  results_list[[drug_name]] <- drug_df
}

# ── 6. 전체 합본 저장 ─────────────────────────────────────────────────────────

results_all <- bind_rows(results_list)
all_path    <- file.path(OUT_DIR, sprintf("mc_results_all_n%d.csv", N))
write_csv(results_all, all_path)

cat(sprintf("\n✅ mc_run_single.R 완료\n"))
cat(sprintf("   Total rows : %d\n", nrow(results_all)))
cat(sprintf("   Saved      : %s\n", all_path))
cat("▶ 다음 단계: mc_main.R\n")

# ── 7. 빠른 검증 ──────────────────────────────────────────────────────────────
cat("\n── Quick validation (median metrics by drug) ──\n")
results_all %>%
  group_by(drug) %>%
  summarise(
    Css_ngmL   = round(median(Css,            na.rm = TRUE) * 1000, 1),
    t50_h      = round(median(t_50pct,        na.rm = TRUE), 1),
    pct24      = round(median(pct_drop_24h,   na.rm = TRUE), 1),
    AUCdef_7d  = round(median(AUC_deficit_7d, na.rm = TRUE), 3),
    .groups    = "drop"
  ) %>%
  rename(
    `Css median (ng/mL)` = Css_ngmL,
    `t50% Css (h)`       = t50_h,
    `%drop @24h`         = pct24,
    `AUC deficit 7d`     = AUCdef_7d
  ) %>%
  print()