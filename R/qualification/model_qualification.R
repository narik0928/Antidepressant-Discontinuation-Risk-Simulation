# =============================================================================
# model_qualification.R
# Phase 3.5 — Model Qualification & Label Verification
#
# Purpose : Verify that rxode2 simulations reproduce FDA label / literature
#           reference values for Tmax, terminal t½, and steady-state Css.
#           This table is required before Monte Carlo (Phase 4).
#
# Approach: Simulate each drug at standard dose →
#             (a) Single-dose run  → Tmax, Cmax
#             (b) Multi-dose run   → Css (last dosing interval average)
#             (c) Washout run      → terminal t½ (log-linear fit)
#           Compare to label ranges → Pass / Concern flag
#
# Output  : outputs/qualification/model_qualification_table.csv
#           outputs/qualification/model_qualification_plot.png
#
# NOTE on pk_parameters.R structure (v4):
#   Single list: pk_params
#   pk_params$sertraline$Ka / $Vd / $CL / $F
#   pk_params$fluoxetine$Ka / $Vd / $CL / $F
#   pk_params$fluoxetine$norfluoxetine$Fm / $Vd_met / $CL_met
#   pk_params$paroxetine$Ka / $Vd / $CL / $F
#   pk_params$venlafaxine_IR$Ka / $Vd / $CL / $F / $ODV$Fm / $ODV$Vd_odv / $ODV$CL_odv
#   pk_params$venlafaxine_XR$Ka  (Ka만 다름; 나머지 venlafaxine_IR과 동일)
# =============================================================================

suppressPackageStartupMessages({
  library(rxode2)
  library(tidyverse)
})

# ── 0. Paths & seed ───────────────────────────────────────────────────────────
BASE    <- "~/nari-research/pkpd-antidepressant-sim"
OUT_DIR <- file.path(BASE, "outputs", "qualification")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

set.seed(2025)   # fixed seed — required for reproducibility in all phases

source(file.path(BASE, "data", "pk_parameters", "pk_parameters.R"))

# ── 1. Parameters (pulled from pk_parameters.R) ───────────────────────────────
# pk_parameters.R exposes a single list: pk_params
# Structure: pk_params$<drug>$<param>
#            pk_params$fluoxetine$norfluoxetine$<param>
#            pk_params$venlafaxine_IR$ODV$<param>
# Simulation settings (dose_mg, tau_h, etc.) are added here, NOT in pk_params.

pk <- list(
  
  sertraline = list(
    Ka      = pk_params$sertraline$Ka,    # 0.50 /h
    Vd      = pk_params$sertraline$Vd,    # 1400 L
    CL      = pk_params$sertraline$CL,    # 37.3 L/h
    F       = pk_params$sertraline$F,     # 1.0 (apparent, F 미확립)
    dose_mg = 50,     # standard dose (mg)
    tau_h   = 24,     # dosing interval (h)
    n_doses = 14,     # doses to reach SS  (≥ 5 × t½ / tau)
    wash_h  = 200     # washout obs after last dose; ≥ 5 × 26h
  ),
  
  fluoxetine = list(
    Ka      = pk_params$fluoxetine$Ka,                       # 0.72 /h
    Vd_p    = pk_params$fluoxetine$Vd,                       # 2310 L
    CL_p    = pk_params$fluoxetine$CL,                       # 10.5 L/h
    F       = pk_params$fluoxetine$F,                        # 0.70
    Fm      = pk_params$fluoxetine$norfluoxetine$Fm,          # 0.72
    Vd_m    = pk_params$fluoxetine$norfluoxetine$Vd_met,      # 2310 L
    CL_m    = pk_params$fluoxetine$norfluoxetine$CL_met,      # 7.18 L/h
    dose_mg = 20,
    tau_h   = 24,
    n_doses = 60,     # 60일: NFX 누적 필요 (t½_nfx = 223h)
    wash_h  = 1500    # ≥ 5 × 223h
  ),
  
  paroxetine = list(
    Ka      = pk_params$paroxetine$Ka,    # 0.580 /h
    Vd      = pk_params$paroxetine$Vd,    # 609 L
    CL      = pk_params$paroxetine$CL,    # 20.1 L/h
    F       = pk_params$paroxetine$F,     # 0.50
    dose_mg = 20,
    tau_h   = 24,
    n_doses = 14,
    wash_h  = 150     # ≥ 5 × 21h
  ),
  
  venlafaxine_IR = list(
    Ka      = pk_params$venlafaxine_IR$Ka,          # 1.50 /h
    Vd_p    = pk_params$venlafaxine_IR$Vd,           # 525 L
    CL_p    = pk_params$venlafaxine_IR$CL,           # 91 L/h
    F       = pk_params$venlafaxine_IR$F,            # 0.45
    Fm      = pk_params$venlafaxine_IR$ODV$Fm,       # 0.80
    Vd_m    = pk_params$venlafaxine_IR$ODV$Vd_odv,   # 399 L
    CL_m    = pk_params$venlafaxine_IR$ODV$CL_odv,   # 28 L/h
    dose_mg = 25,    # 25mg TID = 75mg/day (MC, Figure 2와 동일)
    tau_h   = 8,     # IR: TID dosing
    n_doses = 21,    # 7일 × 3회/일
    wash_h  = 80     # ≥ 5 × 11h (ODV t½ 기준)
  ),
  
  venlafaxine_XR = list(
    Ka      = pk_params$venlafaxine_XR$Ka,           # 0.25 /h
    Vd_p    = pk_params$venlafaxine_XR$Vd,           # 525 L
    CL_p    = pk_params$venlafaxine_XR$CL,           # 91 L/h
    F       = pk_params$venlafaxine_XR$F,            # 0.45
    Fm      = pk_params$venlafaxine_XR$ODV$Fm,       # 0.80
    Vd_m    = pk_params$venlafaxine_XR$ODV$Vd_odv,   # 399 L
    CL_m    = pk_params$venlafaxine_XR$ODV$CL_odv,   # 28 L/h
    dose_mg = 75,
    tau_h   = 24,     # XR: QD dosing
    n_doses = 7,
    wash_h  = 80
  )
)

# ── 2. ODE Model Definitions ──────────────────────────────────────────────────

## 1-compartment oral (Sertraline, Paroxetine)
ode_1cmt <- rxode2({
  d/dt(depot) = -Ka * depot
  d/dt(Cp)    =  Ka * depot / Vd - (CL / Vd) * Cp
})

## 1-compartment oral + active metabolite (Fluoxetine/NFX, Venlafaxine/ODV)
## Metabolite formation: Fm fraction of parent elimination
## Both Cp and Cm in mg/L; Ctotal = unweighted sum (base-case assumption)
ode_parent_met <- rxode2({
  d/dt(depot) = -Ka  * depot
  d/dt(Cp)    =  Ka  * depot / Vd_p - (CL_p / Vd_p) * Cp
  d/dt(Cm)    =  Fm  * CL_p  * Cp   / Vd_m - (CL_m  / Vd_m) * Cm
  Ctotal      =  Cp + Cm
})

# ── 3. Helper Functions ───────────────────────────────────────────────────────

# Extract Tmax (h) and Cmax (mg/L) from a simulation data frame
get_tmax_cmax <- function(df, conc_col = "Cp") {
  idx  <- which.max(df[[conc_col]])
  list(Tmax_h = round(df$time[idx], 2),
       Cmax   = round(df[[conc_col]][idx], 5))
}

# Estimate terminal t½ via log-linear regression on washout curve.
# Uses the terminal `frac` fraction of time points above `floor` × Cmax.
get_thalf_sim <- function(df, conc_col = "Cp",
                          frac = 0.35, floor = 0.005) {
  cmax <- max(df[[conc_col]], na.rm = TRUE)
  df_t <- df %>% filter(.data[[conc_col]] > cmax * floor,
                        .data[[conc_col]] > 0)
  n <- nrow(df_t)
  if (n < 10) return(NA_real_)
  terminal  <- df_t[ceiling(n * (1 - frac)) : n, ]
  fit       <- lm(log(terminal[[conc_col]]) ~ terminal$time)
  slope     <- coef(fit)[2]
  if (is.na(slope) || slope >= 0) return(NA_real_)
  round(-log(2) / slope, 1)
}

# Mean Css (mg/L) from the last dosing interval at steady state
# Uses exact time interval filter — more reliable than tail() for variable dt
get_css <- function(df, conc_col = "Cp") {
  mean(df[[conc_col]], na.rm = TRUE)
}

# Build event table: multi-dose to SS, then washout
make_events <- function(dose_F, tau_h, n_doses, wash_h, dt = 0.25) {
  # dose_F = bioavailability-adjusted dose (F × nominal dose)
  ss_end <- n_doses * tau_h
  et(amt  = dose_F,
     ii   = tau_h,
     addl = n_doses - 1,
     time = 0) %>%
    et(seq(0, ss_end + wash_h, by = dt))
}

# ── 4. Run Simulations & Extract Metrics ─────────────────────────────────────

results <- list()

# ── 4.1 Sertraline ────────────────────────────────────────────────────────────
p  <- pk$sertraline
ev <- make_events(p$F * p$dose_mg, p$tau_h, p$n_doses, p$wash_h)

sim_sert <- rxSolve(
  ode_1cmt,
  params = c(Ka = p$Ka, Vd = p$Vd, CL = p$CL),
  events = ev,
  inits  = c(depot = 0, Cp = 0)
) %>% as.data.frame()

# Single-dose subset for Tmax
sim_sert_sd <- sim_sert %>% filter(time < p$tau_h)

# Washout subset for t½ (after last dose)
t_last_dose <- (p$n_doses - 1) * p$tau_h
sim_sert_wash <- sim_sert %>% filter(time >= t_last_dose)
sim_sert_wash <- sim_sert_wash %>%
  mutate(time = time - t_last_dose)  # re-zero to post-last-dose

results[["Sertraline"]] <- list(
  Tmax_sim   = get_tmax_cmax(sim_sert_sd)$Tmax_h,
  thalf_sim  = get_thalf_sim(sim_sert_wash),
  Css_sim_ng = round(get_css(
    sim_sert %>% filter(time >= (p$n_doses - 1) * p$tau_h,
                        time <= p$n_doses * p$tau_h)) * 1000, 1),
  thalf_analytical = round(log(2) * p$Vd / p$CL, 1)
)

# ── 4.2 Fluoxetine ────────────────────────────────────────────────────────────
p  <- pk$fluoxetine
ev <- make_events(p$F * p$dose_mg, p$tau_h, p$n_doses, p$wash_h,
                  dt = 1)  # coarser dt for long simulation

sim_flx <- rxSolve(
  ode_parent_met,
  params = c(Ka = p$Ka, Vd_p = p$Vd_p, CL_p = p$CL_p,
             Fm  = p$Fm, Vd_m = p$Vd_m, CL_m = p$CL_m),
  events = ev,
  inits  = c(depot = 0, Cp = 0, Cm = 0)
) %>% as.data.frame()

sim_flx_sd    <- sim_flx %>% filter(time < p$tau_h)
t_last_dose   <- (p$n_doses - 1) * p$tau_h
sim_flx_wash  <- sim_flx %>%
  filter(time >= t_last_dose) %>%
  mutate(time = time - t_last_dose)

results[["Fluoxetine"]] <- list(
  Tmax_sim   = get_tmax_cmax(sim_flx_sd, "Cp")$Tmax_h,
  thalf_sim  = get_thalf_sim(sim_flx_wash, "Cp"),
  Css_sim_ng = round(get_css(
    sim_flx %>% filter(time >= (p$n_doses - 1) * p$tau_h,
                       time <= p$n_doses * p$tau_h), "Ctotal") * 1000, 1),
  # NOTE: Css based on Ctotal (Cp + NFX) to match Monte Carlo convention
  thalf_analytical = round(log(2) * p$Vd_p / p$CL_p, 1)
)

results[["Norfluoxetine"]] <- list(
  Tmax_sim   = NA,   # metabolite — Tmax not applicable
  thalf_sim  = get_thalf_sim(sim_flx_wash, "Cm"),
  Css_sim_ng = round(get_css(
    sim_flx %>% filter(time >= (p$n_doses - 1) * p$tau_h,
                       time <= p$n_doses * p$tau_h), "Cm") * 1000, 1),
  thalf_analytical = round(log(2) * p$Vd_m / p$CL_m, 1)
)

# ── 4.3 Paroxetine ────────────────────────────────────────────────────────────
p  <- pk$paroxetine
ev <- make_events(p$F * p$dose_mg, p$tau_h, p$n_doses, p$wash_h)

sim_par <- rxSolve(
  ode_1cmt,
  params = c(Ka = p$Ka, Vd = p$Vd, CL = p$CL),
  events = ev,
  inits  = c(depot = 0, Cp = 0)
) %>% as.data.frame()

sim_par_sd   <- sim_par %>% filter(time < p$tau_h)
t_last_dose  <- (p$n_doses - 1) * p$tau_h
sim_par_wash <- sim_par %>%
  filter(time >= t_last_dose) %>%
  mutate(time = time - t_last_dose)

results[["Paroxetine"]] <- list(
  Tmax_sim   = get_tmax_cmax(sim_par_sd)$Tmax_h,
  thalf_sim  = get_thalf_sim(sim_par_wash),
  Css_sim_ng = round(get_css(
    sim_par %>% filter(time >= (p$n_doses - 1) * p$tau_h,
                       time <= p$n_doses * p$tau_h)) * 1000, 1),
  thalf_analytical = round(log(2) * p$Vd / p$CL, 1)
)

# ── 4.4 Venlafaxine IR ────────────────────────────────────────────────────────
p  <- pk$venlafaxine_IR
ev <- make_events(p$F * p$dose_mg, p$tau_h, p$n_doses, p$wash_h)

sim_ven_IR <- rxSolve(
  ode_parent_met,
  params = c(Ka = p$Ka, Vd_p = p$Vd_p, CL_p = p$CL_p,
             Fm  = p$Fm, Vd_m = p$Vd_m, CL_m = p$CL_m),
  events = ev,
  inits  = c(depot = 0, Cp = 0, Cm = 0)
) %>% as.data.frame()

sim_ven_IR_sd   <- sim_ven_IR %>% filter(time < p$tau_h)
t_last_dose     <- (p$n_doses - 1) * p$tau_h
sim_ven_IR_wash <- sim_ven_IR %>%
  filter(time >= t_last_dose) %>%
  mutate(time = time - t_last_dose)

results[["Venlafaxine IR"]] <- list(
  Tmax_sim   = get_tmax_cmax(sim_ven_IR_sd, "Cp")$Tmax_h,
  thalf_sim  = get_thalf_sim(sim_ven_IR_wash, "Cp"),
  Css_sim_ng = round(get_css(
    sim_ven_IR %>% filter(time >= (p$n_doses - 1) * p$tau_h,
                          time <= p$n_doses * p$tau_h), "Ctotal") * 1000, 1),
  # NOTE: Css based on Ctotal (VEN + ODV) to match Monte Carlo convention
  thalf_analytical = round(log(2) * p$Vd_p / p$CL_p, 1)
)

# ── 4.5 Venlafaxine XR ───────────────────────────────────────────────────────
p  <- pk$venlafaxine_XR
ev <- make_events(p$F * p$dose_mg, p$tau_h, p$n_doses, p$wash_h)

sim_ven_XR <- rxSolve(
  ode_parent_met,
  params = c(Ka = p$Ka, Vd_p = p$Vd_p, CL_p = p$CL_p,
             Fm  = p$Fm, Vd_m = p$Vd_m, CL_m = p$CL_m),
  events = ev,
  inits  = c(depot = 0, Cp = 0, Cm = 0)
) %>% as.data.frame()

sim_ven_XR_sd   <- sim_ven_XR %>% filter(time < p$tau_h)
t_last_dose     <- (p$n_doses - 1) * p$tau_h
sim_ven_XR_wash <- sim_ven_XR %>%
  filter(time >= t_last_dose) %>%
  mutate(time = time - t_last_dose)

results[["Venlafaxine XR"]] <- list(
  Tmax_sim   = get_tmax_cmax(sim_ven_XR_sd, "Cp")$Tmax_h,
  thalf_sim  = get_thalf_sim(sim_ven_XR_wash, "Cp"),
  Css_sim_ng = round(get_css(
    sim_ven_XR %>% filter(time >= (p$n_doses - 1) * p$tau_h,
                          time <= p$n_doses * p$tau_h), "Ctotal") * 1000, 1),
  # NOTE: Css based on Ctotal (VEN + ODV) to match Monte Carlo convention
  thalf_analytical = round(log(2) * p$Vd_p / p$CL_p, 1)
)

results[["ODV"]] <- list(
  Tmax_sim   = NA,
  thalf_sim  = get_thalf_sim(sim_ven_XR_wash, "Cm"),   # from XR washout
  Css_sim_ng = NA,
  thalf_analytical = round(log(2) * pk$venlafaxine_XR$Vd_m /
                             pk$venlafaxine_XR$CL_m, 1)
)

# ── 5. Label Reference Table ──────────────────────────────────────────────────
# Sources: FDA labels (Zoloft, Prozac, Paxil, Effexor XR 2012), NZ Medicines
# Css ranges: steady-state at standard doses; highly variable across subjects.
# Venlafaxine Css: parent only (combined moiety would be higher via ODV).

label <- tribble(
  ~drug,             ~std_dose_str, ~Tmax_lo, ~Tmax_hi,
  ~thalf_lo, ~thalf_hi,
  ~Css_lo,   ~Css_hi,    ~label_source,
  
  "Sertraline",      "50 mg QD",    4.5,  8.4,   22,   36,    20,   60,
  "FDA label (Zoloft); Monfort 2024",
  
  "Fluoxetine",      "20 mg QD",    6.0,  8.0,  120,  192,    91,  604,
  "FDA label (Prozac); NZ Medicines; Ctotal lower bound = FLX parent lower bound",
  
  "Norfluoxetine",   "—",           NA,   NA,   192,  384,    NA,   NA,
  "FDA label (Prozac)",
  
  "Paroxetine",      "20 mg QD",    3.0,  8.0,    7,   65,    30,  100,
  "FDA label (Paxil); Li 2022 (Tmax only)",
  
  "Venlafaxine IR",  "25 mg TID",   1.0,  2.5,    3,    7,    40,  200,
  "FDA label (Effexor XR 2012); Ctotal (VEN+ODV); range estimated from parent+metabolite",
  
  "Venlafaxine XR",  "75 mg QD",    4.0,  8.0,    3,    7,    40,  200,
  "FDA label (Effexor XR 2012); Ctotal (VEN+ODV); range estimated from parent+metabolite",
  
  "ODV",             "—",           NA,   NA,     9,   13,    NA,   NA,
  "FDA label (Effexor XR 2012)"
)
# Note: Css ranges for venlafaxine reflect parent compound only at 75mg.
# The combined venlafaxine+ODV Css would be substantially higher.
# CL values from literature may reflect apparent CL/F — see Methods note.

# ── 6. Assemble Qualification Table ──────────────────────────────────────────

drug_order <- c("Sertraline", "Fluoxetine", "Norfluoxetine",
                "Paroxetine", "Venlafaxine IR", "Venlafaxine XR", "ODV")

qual_rows <- map_dfr(drug_order, function(d) {
  r   <- results[[d]]
  lab <- label %>% filter(drug == d)
  
  # Pass/Concern logic ─────────────────────────────────────────────────────
  # t½ is the primary check; Tmax and Css are secondary.
  
  flag_thalf <- case_when(
    is.na(r$thalf_sim)                              ~ "N/A",
    r$thalf_sim >= lab$thalf_lo &
      r$thalf_sim <= lab$thalf_hi                   ~ "PASS",
    r$thalf_sim  < lab$thalf_lo * 0.75 |
      r$thalf_sim  > lab$thalf_hi * 1.25            ~ "CONCERN",
    TRUE                                            ~ "BORDERLINE"
  )
  
  flag_tmax <- case_when(
    is.na(r$Tmax_sim) | is.na(lab$Tmax_lo)         ~ "N/A",
    r$Tmax_sim >= lab$Tmax_lo &
      r$Tmax_sim <= lab$Tmax_hi                     ~ "PASS",
    TRUE                                            ~ "BORDERLINE"
  )
  
  flag_css <- case_when(
    is.na(r$Css_sim_ng) | is.na(lab$Css_lo)        ~ "N/A",
    r$Css_sim_ng >= lab$Css_lo &
      r$Css_sim_ng <= lab$Css_hi                    ~ "PASS",
    r$Css_sim_ng < lab$Css_lo * 0.5 |
      r$Css_sim_ng > lab$Css_hi * 2                 ~ "CONCERN",
    TRUE                                            ~ "BORDERLINE"
  )
  
  overall <- case_when(
    flag_thalf == "CONCERN" | flag_css == "CONCERN" ~ "CONCERN — review",
    flag_thalf == "PASS"                            ~ "PASS",
    flag_thalf == "BORDERLINE"                      ~ "BORDERLINE",
    TRUE                                            ~ "N/A"
  )
  
  tibble(
    Drug                  = d,
    Std_dose              = lab$std_dose_str,
    Tmax_sim_h            = r$Tmax_sim,
    Tmax_label_range_h    = ifelse(is.na(lab$Tmax_lo), "—",
                                   paste0(lab$Tmax_lo, "–", lab$Tmax_hi)),
    Tmax_flag             = flag_tmax,
    thalf_sim_h           = r$thalf_sim,
    thalf_analytical_h    = r$thalf_analytical,
    thalf_label_range_h   = paste0(lab$thalf_lo, "–", lab$thalf_hi),
    thalf_flag            = flag_thalf,
    Css_sim_ngmL          = r$Css_sim_ng,
    Css_label_range_ngmL  = ifelse(is.na(lab$Css_lo), "—",
                                   paste0(lab$Css_lo, "–", lab$Css_hi)),
    Css_flag              = flag_css,
    Overall               = overall,
    Source                = lab$label_source
  )
})

# ── 7. Save CSV ───────────────────────────────────────────────────────────────
csv_path <- file.path(OUT_DIR, "model_qualification_table.csv")
write_csv(qual_rows, csv_path)
message("✓ Qualification table saved: ", csv_path)

# ── 8. Print to console ───────────────────────────────────────────────────────
cat("\n═══════════════════════════════════════════════════════\n")
cat(" MODEL QUALIFICATION TABLE  (Phase 3.5)\n")
cat("═══════════════════════════════════════════════════════\n\n")

print_cols <- c("Drug", "Tmax_sim_h", "Tmax_label_range_h", "Tmax_flag",
                "thalf_sim_h", "thalf_analytical_h", "thalf_label_range_h",
                "thalf_flag", "Css_sim_ngmL", "Css_label_range_ngmL",
                "Css_flag", "Overall")
print(qual_rows[, print_cols], n = Inf, width = Inf)

# ── 9. Visualisation: t½ comparison ──────────────────────────────────────────
# Bar chart: simulated t½ vs label range (ribbon)
# Primary metric: t½ is most critical for discontinuation analysis

plot_data <- qual_rows %>%
  filter(!is.na(thalf_sim_h)) %>%
  left_join(label %>%
              select(drug, thalf_lo, thalf_hi),
            by = c("Drug" = "drug")) %>%
  mutate(Drug = factor(Drug, levels = rev(drug_order)))

p_thalf <- ggplot(plot_data, aes(y = Drug)) +
  # Label range shading
  geom_rect(aes(xmin = thalf_lo, xmax = thalf_hi,
                ymin = as.numeric(Drug) - 0.4,
                ymax = as.numeric(Drug) + 0.4),
            fill = "#B8D4E8", alpha = 0.6) +
  # Simulated value
  geom_point(aes(x = thalf_sim_h, colour = thalf_flag),
             size = 4, shape = 18) +
  # Analytical value (open diamond for comparison)
  geom_point(aes(x = thalf_analytical_h),
             size = 3, shape = 5, colour = "grey40") +
  scale_colour_manual(
    values = c("PASS" = "#2E8B57", "BORDERLINE" = "#DAA520",
               "CONCERN" = "#C0392B", "N/A" = "grey60"),
    name = "t½ flag"
  ) +
  scale_x_log10(
    breaks = c(1, 3, 5, 10, 20, 50, 100, 200, 400),
    labels = scales::comma
  ) +
  labs(
    title    = "Model Qualification: Simulated t½ vs FDA Label Range",
    subtitle = "Blue band = label range  |  ◆ = simulated  |  ⬦ = analytical (0.693·Vd/CL)",
    x        = "Terminal half-life (h, log scale)",
    y        = NULL,
    caption  = paste0("Phase 3.5 qualification  |  seed = 2025  |  ",
                      format(Sys.Date(), "%Y-%m-%d"))
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position  = "bottom",
    plot.title       = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

plot_path <- file.path(OUT_DIR, "model_qualification_plot.png")
ggsave(plot_path, p_thalf, width = 10, height = 6, dpi = 300)
message("✓ Qualification plot saved: ", plot_path)

# ── 10. Flag summary ──────────────────────────────────────────────────────────
cat("\n── Overall flags ──\n")
qual_rows %>%
  select(Drug, thalf_flag, Tmax_flag, Css_flag, Overall) %>%
  print(n = Inf)

n_concern <- sum(qual_rows$Overall == "CONCERN — review", na.rm = TRUE)
if (n_concern == 0) {
  cat("\n✓ All drugs PASS or BORDERLINE. Safe to proceed to Phase 4.\n\n")
} else {
  cat(sprintf(
    "\n⚠ %d drug(s) flagged CONCERN. Review parameters before Phase 4.\n\n",
    n_concern
  ))
}