# =============================================================================
# mc_params.R
# Phase 4 — Monte Carlo Parameter Sampling
#
# Purpose : 약물별 가상 환자 파라미터 샘플링 테이블 생성
#           CL, Vd → LogNormal IIV 적용
#           CYP2D6 phenotype → CL 보정 계수 추가 적용
#
# Sampling model:
#   CL_i  = CL_pop × CYP2D6_factor × exp(rnorm(0, CV_CL))
#   Vd_i  = Vd_pop ×                  exp(rnorm(0, CV_Vd))
#   (LogNormal: mean = pop value, CV = IIV)
#
# Output:
#   outputs/mc/mc_params_<drug>_n<N>.csv  (약물별)
#   outputs/mc/mc_params_all_n<N>.csv     (전체 합본)
#
# 실행: source("R/simulation/mc_params.R")
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

# ── 0. Setup ──────────────────────────────────────────────────────────────────
BASE    <- "~/nari-research/pkpd-antidepressant-sim"
OUT_DIR <- file.path(BASE, "outputs", "mc")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

set.seed(2025)   # 고정 seed — 재현성 필수 (논문 Methods에 명시)

N <- 100      # 전체 환자 수 (로컬 테스트 시 100으로 변경)

source(file.path(BASE, "data", "pk_parameters", "pk_parameters.R"))

# ── 1. CYP2D6 Phenotype 배정 함수 ─────────────────────────────────────────────
# 적용 약물: Fluoxetine, Paroxetine, Venlafaxine (Sertraline은 영향 제한적)
# 참고: cyp2d6 리스트는 pk_parameters.R에서 로드됨

assign_cyp2d6 <- function(n) {
  sample(
    x       = c("poor", "normal", "rapid"),
    size    = n,
    replace = TRUE,
    prob    = c(
      cyp2d6$poor_metabolizer$freq,    # 0.08
      cyp2d6$normal_metabolizer$freq,  # 0.70
      cyp2d6$rapid_metabolizer$freq    # 0.22
    )
  )
}

get_cyp2d6_factor <- function(phenotype) {
  case_when(
    phenotype == "poor"   ~ cyp2d6$poor_metabolizer$CL_factor,    # 0.30
    phenotype == "normal" ~ cyp2d6$normal_metabolizer$CL_factor,  # 1.00
    phenotype == "rapid"  ~ cyp2d6$rapid_metabolizer$CL_factor    # 1.60
  )
}

# ── 2. LogNormal 샘플링 함수 ──────────────────────────────────────────────────
# LogNormal 파라미터화: X ~ LN(μ, σ²)
# mean(X) = pop_value 유지하면서 CV 반영
# σ = sqrt(log(CV² + 1)) ≈ CV (CV 작을 때)
# μ = log(pop_value) - σ²/2

sample_lognormal <- function(n, pop_value, cv) {
  sigma <- sqrt(log(cv^2 + 1))
  mu    <- log(pop_value) - sigma^2 / 2
  rlnorm(n, meanlog = mu, sdlog = sigma)
}

# ── 3. 약물별 파라미터 샘플링 ─────────────────────────────────────────────────

## ── 3.1 Sertraline ────────────────────────────────────────────────────────────
# CYP2D6 영향 제한적 → CYP2D6 factor 미적용 (cyp2d6 = "normal" 고정)
sample_sertraline <- function(n) {
  p <- pk_params$sertraline
  tibble(
    id            = seq_len(n),
    drug          = "Sertraline",
    cyp2d6        = "normal",          # CYP2D6 영향 제한적
    CL_pop        = p$CL,
    Vd_pop        = p$Vd,
    CL_sampled    = sample_lognormal(n, p$CL, p$CV_CL),
    Vd_sampled    = sample_lognormal(n, p$Vd, p$CV_Vd),
    CL_met_sampled = NA_real_,
    Vd_met_sampled = NA_real_
  )
}

## ── 3.2 Fluoxetine ────────────────────────────────────────────────────────────
# CYP2D6 → parent CL에 적용
# Norfluoxetine CL_met도 별도 샘플링 (CV_CL_met 적용)
sample_fluoxetine <- function(n) {
  p    <- pk_params$fluoxetine
  nfx  <- pk_params$fluoxetine$norfluoxetine
  cyp  <- assign_cyp2d6(n)
  fac  <- get_cyp2d6_factor(cyp)
  
  tibble(
    id             = seq_len(n),
    drug           = "Fluoxetine",
    cyp2d6         = cyp,
    CL_pop         = p$CL,
    Vd_pop         = p$Vd,
    CL_sampled     = sample_lognormal(n, p$CL, p$CV_CL) * fac,
    Vd_sampled     = sample_lognormal(n, p$Vd, p$CV_Vd),
    CL_met_sampled = sample_lognormal(n, nfx$CL_met, nfx$CV_CL_met),
    Vd_met_sampled = rep(nfx$Vd_met, n)   # Vd_met IIV 데이터 없음 → pop값 고정
  )
}

## ── 3.3 Paroxetine ────────────────────────────────────────────────────────────
# CYP2D6 → CL에 적용 (CV_CL 높음 = 0.60, 비선형 CYP2D6 반영)
sample_paroxetine <- function(n) {
  p   <- pk_params$paroxetine
  cyp <- assign_cyp2d6(n)
  fac <- get_cyp2d6_factor(cyp)
  
  tibble(
    id             = seq_len(n),
    drug           = "Paroxetine",
    cyp2d6         = cyp,
    CL_pop         = p$CL,
    Vd_pop         = p$Vd,
    CL_sampled     = sample_lognormal(n, p$CL, p$CV_CL) * fac,
    Vd_sampled     = sample_lognormal(n, p$Vd, p$CV_Vd),
    CL_met_sampled = NA_real_,
    Vd_met_sampled = NA_real_
  )
}

## ── 3.4 Venlafaxine IR ────────────────────────────────────────────────────────
# CYP2D6 → parent CL에 적용 (ODV 형성 비율에 영향)
# ODV CL_met 별도 샘플링
sample_venlafaxine_IR <- function(n) {
  p   <- pk_params$venlafaxine_IR
  odv <- pk_params$venlafaxine_IR$ODV
  cyp <- assign_cyp2d6(n)
  fac <- get_cyp2d6_factor(cyp)
  
  tibble(
    id             = seq_len(n),
    drug           = "Venlafaxine_IR",
    cyp2d6         = cyp,
    CL_pop         = p$CL,
    Vd_pop         = p$Vd,
    CL_sampled     = sample_lognormal(n, p$CL, p$CV_CL) * fac,
    Vd_sampled     = sample_lognormal(n, p$Vd, p$CV_Vd),
    CL_met_sampled = sample_lognormal(n, odv$CL_odv, odv$CV_CL_odv),
    Vd_met_sampled = rep(odv$Vd_odv, n)
  )
}

## ── 3.5 Venlafaxine XR ────────────────────────────────────────────────────────
# PK 파라미터 IR과 동일 (Ka만 다름 — 시뮬레이션 단계에서 처리)
sample_venlafaxine_XR <- function(n) {
  p   <- pk_params$venlafaxine_XR
  odv <- pk_params$venlafaxine_XR$ODV
  cyp <- assign_cyp2d6(n)
  fac <- get_cyp2d6_factor(cyp)
  
  tibble(
    id             = seq_len(n),
    drug           = "Venlafaxine_XR",
    cyp2d6         = cyp,
    CL_pop         = p$CL,
    Vd_pop         = p$Vd,
    CL_sampled     = sample_lognormal(n, p$CL, p$CV_CL) * fac,
    Vd_sampled     = sample_lognormal(n, p$Vd, p$CV_Vd),
    CL_met_sampled = sample_lognormal(n, odv$CL_odv, odv$CV_CL_odv),
    Vd_met_sampled = rep(odv$Vd_odv, n)
  )
}

# ── 4. 전체 샘플링 실행 ───────────────────────────────────────────────────────
cat(sprintf("\n▶ Sampling %d virtual patients per drug...\n", N))

params_all <- bind_rows(
  sample_sertraline(N),
  sample_fluoxetine(N),
  sample_paroxetine(N),
  sample_venlafaxine_IR(N),
  sample_venlafaxine_XR(N)
)

cat(sprintf("✓ Total rows: %d (5 drugs × %d patients)\n", nrow(params_all), N))

# ── 5. 샘플링 검증 ────────────────────────────────────────────────────────────
# 샘플링된 CL 분포가 pop값 근처에 있는지 확인
# E[CL_sampled / CL_pop] ≈ 1.0 이어야 함

cat("\n── Sampling verification (CL_sampled / CL_pop mean) ──\n")
params_all %>%
  filter(!is.na(CL_sampled)) %>%
  group_by(drug) %>%
  summarise(
    CL_pop         = mean(CL_pop),
    CL_mean_sim    = round(mean(CL_sampled), 2),
    CL_ratio       = round(mean(CL_sampled) / mean(CL_pop), 3),
    CL_cv_obs      = round(sd(CL_sampled) / mean(CL_sampled), 3),
    .groups = "drop"
  ) %>%
  print()

cat("\n── CYP2D6 phenotype distribution ──\n")
params_all %>%
  filter(drug != "Sertraline") %>%
  count(drug, cyp2d6) %>%
  mutate(pct = round(n / N * 100, 1)) %>%
  print(n = Inf)

# ── 6. 저장 ───────────────────────────────────────────────────────────────────

# 약물별 개별 저장
walk(unique(params_all$drug), function(d) {
  df   <- params_all %>% filter(drug == d)
  path <- file.path(OUT_DIR,
                    sprintf("mc_params_%s_n%d.csv",
                            tolower(gsub("_", "", d)), N))
  write_csv(df, path)
  cat(sprintf("✓ Saved: %s\n", path))
})

# 전체 합본 저장
all_path <- file.path(OUT_DIR, sprintf("mc_params_all_n%d.csv", N))
write_csv(params_all, all_path)
cat(sprintf("✓ Saved (all): %s\n", all_path))

cat(sprintf("\n✅ mc_params.R 완료 — %d patients × 5 drugs = %d rows\n",
            N, nrow(params_all)))
cat("▶ 다음 단계: mc_run_single.R\n\n")