# ============================================================
# PK Parameters — Antidepressant Discontinuation Risk Study
# Author  : Nari Kim
# Version : v2 (검증 완료 2025-05)
# Changes : Sertraline CL/F 수정 (37→88 L/h), F=1.0 처리
#           Fluoxetine Vd 수정 (2450→2310 L), F=70%, NFX t½ 수정
# Sources : FDA label (Zoloft, Prozac, Effexor XR),
#           Monfort 2024 BJCP, Li 2022 Front Pharmacol,
#           Sagahon-Azua 2021 PRP, Troy 1994 J Clin Pharmacol
# Units   : Ka (1/h), Vd (L), CL (L/h), F (fraction), t½ (h)
# ============================================================

pk_params <- list(
  
  # ── Sertraline ─────────────────────────────────────────────
  # F 미확립 (IV 비교 연구 없음, Markowitz 2002: >44% 추정)
  # → F = 1.0 처리, CL/F = 88 L/h 사용 (apparent oral clearance)
  # CL/F 출처: Monfort 2024 BJCP: 1.09–1.41 L/h/kg → 1.26 × 70kg = 88 L/h
  sertraline = list(
    drug        = "Sertraline",
    class       = "SSRI",
    formulation = "IR",
    Ka          = 0.50,    # 1/h | Tmax 4.5–8.4h 역산 | FDA label
    Vd          = 1400,    # L   | 20 L/kg × 70kg | INCHEM/WHO, Monfort 2024
    CL          = 88,      # L/h | CL/F = 1.26 L/h/kg × 70kg | Monfort 2024 BJCP
    F           = 1.0,     # apparent (F 미확립 — CL/F에 통합)
    t_half      = 26,      # h   | FDA label (Zoloft)
    CV_CL       = 0.40,    # IIV | 40% (Monfort 2024)
    CV_Vd       = 0.35,
    model       = "1-compartment",
    F_note      = "절대 F 미확립. IV 비교 연구 없음 (Markowitz 2002: >44% 추정). CL/F 통합 처리.",
    source      = "FDA label (Zoloft); Monfort 2024 BJCP; Markowitz 2002 Clin Pharmacokinet"
  ),
  
  # ── Fluoxetine ─────────────────────────────────────────────
  # 비선형 PK (CYP2D6 포화) — 레퍼런스 약물로만 활용
  # Vd: 33 L/kg × 70kg = 2310 L (NZ Medicines: 33 L/kg 중앙값)
  # F: 60–80% 범위, 70% 중간값 사용
  # Norfluoxetine t½: FDA label 반복 투여 9.3일 = 223h
  fluoxetine = list(
    drug        = "Fluoxetine",
    class       = "SSRI",
    formulation = "IR",
    Ka          = 0.30,    # 1/h | Pauchard 2011 PopPK (고정값)
    Vd          = 2310,    # L   | 33 L/kg × 70kg | NZ Medicines monograph
    CL          = 10.5,    # L/h | 0.15 L/h/kg × 70kg | NZ Medicines monograph
    F           = 0.70,    # 60–80% 범위 중간값 | Wikipedia/DrugBank
    t_half      = 120,     # h   | 4–6일 chronic → 120h | Clin Pharmacokinet
    CV_CL       = 0.874,   # IIV | 87.4% | Sagahon-Azua 2021 PRP
    CV_Vd       = 0.50,
    model       = "2-compartment + metabolite",
    norfluoxetine = list(
      Ka_met    = 0.30,
      Vd_met    = 2310,
      CL_met    = 0.9,     # L/h | median | Sagahon-Azua 2021 PRP
      t_half    = 223,     # h   | 반복 투여 9.3일 × 24 = 223h | FDA label (Prozac)
      CV_CL_met = 0.564,   # IIV | 56.4% | Sagahon-Azua 2021
      Fm        = 0.72     # 대사 전환율
    ),
    source      = "NZ Medicines monograph; Pauchard 2011; Sagahon-Azua 2021 PRP; FDA label (Prozac)"
  ),
  
  # ── Paroxetine IR ─────────────────────────────────────────
  # 비선형 PK (CYP2D6 용량 포화): 용량 증가 시 CL 감소
  # CL/F = 21.2 L/h, V/F = 8850 L: Li 2022 PopPK (n = 전체 환자군)
  paroxetine_IR = list(
    drug        = "Paroxetine",
    class       = "SSRI",
    formulation = "IR",
    Ka          = 0.908,   # 1/h | Li 2022 Front Pharmacol (고정값)
    Vd          = 8850,    # L   | V/F | Li 2022 Front Pharmacol
    CL          = 21.2,    # L/h | CL/F | Li 2022 Front Pharmacol
    F           = 0.50,    # 약 50% (1st-pass 포화성 대사) | Kaye 1989
    t_half      = 21,      # h   | 평균값; 범위 7–65h (Kaye 1989)
    CV_CL       = 0.60,    # IIV | 높음 (비선형 PK)
    CV_Vd       = 0.50,
    model       = "1-compartment (nonlinear CYP2D6)",
    source      = "Li 2022 Front Pharmacol; Kaye 1989 Br J Clin Pharmacol"
  ),
  
  # ── Paroxetine SR ─────────────────────────────────────────
  # Li 2022: SR 제형 시 V/F IR 대비 66.6% 감소
  # CL은 IR과 동일 (제형이 CL에 영향 없음)
  paroxetine_SR = list(
    drug        = "Paroxetine",
    class       = "SSRI",
    formulation = "SR",
    Ka          = 0.200,   # 1/h | SR 흡수 지연 (추정값, Tmax 연장 기반)
    Vd          = 2955,    # L   | 8850 × (1 - 0.666) = 2955 L | Li 2022
    CL          = 21.2,    # L/h | IR과 동일 | Li 2022
    F           = 0.50,
    t_half      = 21,      # h
    CV_CL       = 0.60,
    CV_Vd       = 0.50,
    model       = "1-compartment (SR)",
    source      = "Li 2022 Front Pharmacol (V/F -66.6% for SR formulation)"
  ),
  
  # ── Venlafaxine IR ────────────────────────────────────────
  # FDA label (Effexor XR): CL = 1.3±0.6 L/h/kg, Vd = 7.5±3.7 L/kg
  # t½ VEN = 5±2h, ODV = 11±2h
  # F = 45% (절대 생체이용률, FDA label)
  venlafaxine_IR = list(
    drug        = "Venlafaxine",
    class       = "SNRI",
    formulation = "IR",
    Ka          = 1.50,    # 1/h | Tmax ~2h 기반 역산 | FDA label
    Vd          = 525,     # L   | 7.5 L/kg × 70kg | FDA label (Effexor XR)
    CL          = 91,      # L/h | 1.3 L/h/kg × 70kg | FDA label (Effexor XR)
    F           = 0.45,    # 절대 생체이용률 | FDA label (Effexor XR)
    t_half      = 5,       # h   | 5±2h | FDA label
    CV_CL       = 0.46,    # IIV | 0.6/1.3 = 46% (SD/mean) | FDA label
    CV_Vd       = 0.49,    # IIV | 3.7/7.5 = 49%
    model       = "1-compartment + metabolite (ODV)",
    ODV = list(
      Vd_odv    = 399,     # L   | 5.7 L/kg × 70kg | FDA label
      CL_odv    = 28,      # L/h | 0.4 L/h/kg × 70kg | FDA label
      t_half    = 11,      # h   | 11±2h | FDA label
      CV_CL_odv = 0.50,    # IIV | 0.2/0.4 = 50%
      Fm        = 0.80     # ODV 전환율 (80%)
    ),
    source      = "FDA label (Effexor XR 2012, accessdata.fda.gov)"
  ),
  
  # ── Venlafaxine XR ────────────────────────────────────────
  # FDA label: XR은 IR과 동일한 생체이용률(F=45%), 흡수 속도만 다름
  # Tmax XR = 5.5h (IR = 2h) → Ka 대폭 감소
  # CL, Vd, t½ 모두 IR과 동일
  venlafaxine_XR = list(
    drug        = "Venlafaxine",
    class       = "SNRI",
    formulation = "XR",
    Ka          = 0.25,    # 1/h | Tmax ~5.5h 기반 역산 | FDA label
    Vd          = 525,     # L   | IR과 동일 | FDA label
    CL          = 91,      # L/h | IR과 동일 | FDA label
    F           = 0.45,    # IR과 동일 | FDA label: "same extent of absorption"
    t_half      = 5,       # h   | IR과 동일
    CV_CL       = 0.46,
    CV_Vd       = 0.49,
    model       = "1-compartment + metabolite (ODV)",
    ODV = list(
      Vd_odv    = 399,
      CL_odv    = 28,
      t_half    = 11,
      CV_CL_odv = 0.50,
      Fm        = 0.80
    ),
    source      = "FDA label (Effexor XR 2012): slower absorption, same extent"
  )
)

# ── CYP2D6 genotype 보정 계수 ─────────────────────────────────
# 적용 약물: Paroxetine, Venlafaxine, Fluoxetine
# 출처: FDA label (Effexor XR), PharmGKB
cyp2d6 <- list(
  poor_metabolizer   = list(freq = 0.08, CL_factor = 0.30),
  normal_metabolizer = list(freq = 0.70, CL_factor = 1.00),
  rapid_metabolizer  = list(freq = 0.22, CL_factor = 1.60)
)

# ── 검증 요약 ─────────────────────────────────────────────────
validation_notes <- list(
  sertraline    = "CL 37→88 L/h (Monfort 2024). F 미확립 → F=1.0, CL/F 통합",
  fluoxetine    = "Vd 2450→2310 L (33 L/kg). F 72→70%. NFX t½ 360→223h (FDA label)",
  paroxetine    = "파라미터 유지. t½ 변동 매우 큼(7-65h) — IIV 60% 반영",
  venlafaxine   = "CL 84→91 L/h (FDA label 정확값). ODV Vd 525→399 L 수정"
)

cat("✅ PK 파라미터 v2 로드 완료 (검증 완료)\n")
cat("약물:", paste(names(pk_params), collapse = " | "), "\n")
cat("\n[검증 변경사항]\n")
for (drug in names(validation_notes)) {
  cat(sprintf("  %-20s: %s\n", drug, validation_notes[[drug]]))
}