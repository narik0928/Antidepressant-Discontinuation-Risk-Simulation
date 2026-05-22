# Antidepressant Discontinuation Risk — PK/PD Simulation Study

**Author:** Nari Kim  
**Institution:** Kingston University London, BSc Pharmacology  
**Target Journal:** Translational and Clinical Pharmacology (TCP)  
**Status:** Phase 5 Complete — Manuscript in preparation

---

## Research Overview

This repository contains the computational pharmacokinetic (PK) simulation framework developed to quantify exposure-loss stress across selected SSRI and SNRI antidepressants under abrupt discontinuation scenarios.

### Research Question
> Can pharmacokinetic properties — including half-life, active metabolite accumulation, formulation, and CYP2D6 phenotype — explain differential exposure-loss stress across selected SSRIs and SNRIs under abrupt discontinuation?

### Key Findings
- Venlafaxine XR shows the fastest relative decline (t50% C0 = 7.2 h, stress rank 1)
- Fluoxetine shows the slowest decline (t50% C0 = 349 h) due to norfluoxetine buffering
- Paroxetine shows the largest CYP2D6 effect (η² = 0.197, large)
- Venlafaxine IR has the highest absolute AUC deficit (20.2 mg/L·h)

---

## Repository Structure

```
.
├── R/
│   ├── models/                    # Phase 3: Deterministic PK models
│   │   ├── sertraline_model.R
│   │   ├── fluoxetine_model.R
│   │   ├── paroxetine_model.R
│   │   └── venlafaxine_model.R
│   │
│   ├── qualification/             # Phase 3.5: Model qualification
│   │   └── model_qualification.R
│   │
│   ├── simulation/                # Phase 4: Monte Carlo simulation
│   │   ├── mc_params.R            # Parameter sampling (n=10,000)
│   │   ├── mc_run_single.R        # Per-patient PK simulation
│   │   └── mc_main.R              # Table 3 & 4 generation
│   │
│   ├── analysis/                  # Phase 5: Statistical analysis
│   │   ├── stability_check.R      # n=1,000/5,000/10,000 convergence
│   │   └── phase5_statistics.R    # Spearman, Kruskal-Wallis, Dunn test
│   │
│   └── figures/                   # Figure generation
│       ├── figure1_discontinuation.R
│       ├── figure2_venlafaxine_IR_XR.R
│       ├── figure3_mc_t50_violin.R
│       ├── figure4_mc_dCdt_ridge.R
│       └── figure5_mc_scatter.R
│
├── data/
│   └── pk_parameters/
│       └── pk_parameters.R        # PK parameters v4 (FDA label verified)
│
└── outputs/
    ├── figures/                   # Figure 1–5 (PNG, 300 dpi)
    ├── qualification/             # Model qualification table & plot
    ├── mc/                        # Monte Carlo results (CSV)
    └── tables/                    # Table 3–7 (CSV)
```

---

## Execution Order

```r
# Step 1. PK model qualification
source("R/qualification/model_qualification.R")

# Step 2. Monte Carlo parameter sampling
source("R/simulation/mc_params.R")          # N <- 10000

# Step 3. Monte Carlo simulation
source("R/simulation/mc_run_single.R")      # N <- 10000

# Step 4. Summary tables
source("R/simulation/mc_main.R")

# Step 5. Stability check
source("R/analysis/stability_check.R")

# Step 6. Statistical analysis
source("R/analysis/phase5_statistics.R")

# Step 7. Figures
source("R/figures/figure1_discontinuation.R")
source("R/figures/figure2_venlafaxine_IR_XR.R")
source("R/figures/figure3_mc_t50_violin.R")
source("R/figures/figure4_mc_dCdt_ridge.R")
source("R/figures/figure5_mc_scatter.R")
```

---

## Drug & PK Parameters (v4)

| Drug | Class | Model | t½ | Source |
|---|---|---|---|---|
| Sertraline | SSRI, IR | 1-cmt | 26 h | FDA label (Zoloft) |
| Fluoxetine +NFX | SSRI, IR | 1-cmt + metabolite | 153 h / 223 h | FDA label (Prozac); NZ Medicines |
| Paroxetine | SSRI, IR | 1-cmt | 21 h | FDA label (Paxil); DrugBank |
| Venlafaxine IR +ODV | SNRI, IR | 1-cmt + metabolite | 5 h / 11 h | FDA label (Effexor XR 2012) |
| Venlafaxine XR +ODV | SNRI, XR | 1-cmt + metabolite | 5 h / 11 h | FDA label (Effexor XR 2012) |

---

## Key PK Metrics

| Metric | Definition |
|---|---|
| C0 | Plasma concentration at time of discontinuation |
| t50% C0 | Time to reach 50% of C0 (linear interpolation) |
| t25% C0 | Time to reach 25% of C0 |
| \|Rel. dC/dt\| @24h | \|(C0 − C24h)\| / (24 × C0) — Css-normalised decline rate |
| AUC deficit 7d | C0 × 168h − area under post-discontinuation curve |

---

## Dependencies

```r
install.packages(c(
  "rxode2",     # ODE solver
  "tidyverse",  # Data manipulation
  "patchwork",  # Figure assembly
  "ggridges",   # Ridge plots
  "ggrepel",    # Label repulsion
  "rstatix"     # Statistical tests
))
```

---

## Reproducibility

- Random seed: `set.seed(2025)` applied in all stochastic steps
- Stability confirmed: CV < 5% at n=1,000; CV < 2% at n=5,000
- R version: 4.3.3 | rxode2: 5.0.2 | ggplot2: 4.0.3

---

## Notes

- All threshold lines use C0 (concentration at discontinuation), not mean Css
- Active metabolite contribution uses unweighted sum (base-case assumption)
- CYP2D6 variability applied to Fluoxetine, Paroxetine, Venlafaxine IR/XR
- Sertraline: CYP2D6 effect limited; single phenotype modelled
- Sensitivity analysis: removing t50% ≥ 1h filter changes estimates by ≤8.3%

---

*Last updated: 2025-05*