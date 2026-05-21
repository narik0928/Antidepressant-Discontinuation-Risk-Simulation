# ============================================================
# PK Parameters — Antidepressant Discontinuation Risk Study
# Author  : Nari Kim
# Version : v4 (2025-05)
# Changes : Paroxetine 전면 수정
#   - Vd: 8850→609 L (8.7 L/kg × 70kg, DrugBank/medicine.com)
#   - CL: 21.2→20.1 L/h (0.693 × 609 / 21, FDA label t½ 역산)
#   - Ka: 0.908→0.580 /h (Tmax 5.2h 맞게 조정, Li 2022 Tmax 기준)
#   - paroxetine_SR 제거 (IR vs SR 중단 후 곡선 동일 → 비교 의미 없음)
#   Li 2022 V/F=8850은 비현실적 (정신과 환자 + 비선형 PK 왜곡) → 미적용
# ============================================================

pk_params <- list(
  
  # ── Sertraline ─────────────────────────────────────────────
  # CL = 0.693 × 1400 / 26 = 37.3 L/h (FDA label t½=26h 역산)
  # F 미확립 → F=1.0, CL/F 통합
  sertraline = list(
    drug        = "Sertraline",
    class       = "SSRI",
    formulation = "IR",
    Ka          = 0.50,    # 1/h | Tmax 4.5–8.4h 역산 | FDA label
    Vd          = 1400,    # L   | 20 L/kg × 70kg | Monfort 2024 BJCP
    CL          = 37.3,    # L/h | 0.693 × 1400 / 26 | FDA label t½ 역산
    F           = 1.0,     # apparent (F 미확립)
    t_half      = 26,      # h   | FDA label
    CV_CL       = 0.40,
    CV_Vd       = 0.35,
    model       = "1-compartment",
    source      = "FDA label (Zoloft); Monfort 2024 BJCP (Vd)"
  ),
  
  # ── Fluoxetine ─────────────────────────────────────────────
  fluoxetine = list(
    drug        = "Fluoxetine",
    class       = "SSRI",
    formulation = "IR",
    Ka          = 0.72,    # 1/h | Tmax 7.1h 역산 | FDA label 6-8h 기준
    Vd          = 2310,    # L   | 33 L/kg × 70kg | NZ Medicines
    CL          = 10.5,    # L/h | 0.15 L/h/kg × 70kg | NZ Medicines
    F           = 0.70,
    t_half      = 153,     # h   | 6.4일 (CL=10.5 역산값, 문헌범위 4-16일 내)
    CV_CL       = 0.874,
    CV_Vd       = 0.50,
    model       = "1-compartment + metabolite (단순화, terminal elimination 집중)",
    norfluoxetine = list(
      Ka_met    = 0.30,
      Vd_met    = 2310,
      CL_met    = 7.18,    # L/h | 0.693 × 2310 / 223 | FDA label t½=223h 역산
      t_half    = 223,     # h   | FDA label 반복투여 9.3일
      CV_CL_met = 0.564,
      Fm        = 0.72
    ),
    source      = "NZ Medicines (CL); FDA label (Prozac, Ka/t½_nfx); Pauchard 2011 (Tmax기준)"
  ),
  
  # ── Paroxetine IR ──────────────────────────────────────────
  # v4 수정: Vd 8850→609 L, CL 21.2→20.1 L/h, Ka 0.908→0.580
  # Li 2022 V/F=8850 미적용 이유:
  #   1) 정신과 환자(psychosis) 특화 집단 → 일반 성인 부적합
  #   2) 비선형 PK를 선형 모델로 fitting → V/F 과도 추정
  #   3) V/F=8850 + CL/F=21.2 → t½=289h: FDA label 21h와 불일치
  # → Vd: DrugBank 8.7 L/kg × 70kg = 609 L
  # → CL: FDA label t½=21h 역산: 0.693 × 609 / 21 = 20.1 L/h
  # → Ka: Tmax 5.2h (Li 2022) 맞게 역산: Ka=0.58/h
  # paroxetine_SR 제거: Vd, CL 동일 → 중단 후 소실 패턴 동일
  paroxetine = list(
    drug        = "Paroxetine",
    class       = "SSRI",
    formulation = "IR",
    Ka          = 0.580,   # 1/h | Tmax 5.2h 역산 | Li 2022 Tmax 기준
    Vd          = 609,     # L   | 8.7 L/kg × 70kg | DrugBank/medicine.com
    CL          = 20.1,    # L/h | 0.693 × 609 / 21 | FDA label t½ 역산
    F           = 0.50,    # 1st-pass 포화성 대사 | Kaye 1989
    t_half      = 21,      # h   | FDA label (범위 7–65h)
    CV_CL       = 0.60,    # IIV 높음 (비선형 CYP2D6)
    CV_Vd       = 0.35,
    model       = "1-compartment",
    CL_note     = "CL = FDA label t½ 역산. Li 2022 CL/F=21.2 미적용 (V/F=8850 왜곡)",
    source      = "FDA label (Paxil); DrugBank Vd=8.7 L/kg; Li 2022 (Ka/Tmax만)"
  ),
  
  # ── Venlafaxine IR ─────────────────────────────────────────
  venlafaxine_IR = list(
    drug        = "Venlafaxine",
    class       = "SNRI",
    formulation = "IR",
    Ka          = 1.50,    # 1/h | Tmax ~2h | FDA label
    Vd          = 525,     # L   | 7.5 L/kg × 70kg | FDA label
    CL          = 91,      # L/h | 1.3 L/h/kg × 70kg | FDA label
    F           = 0.45,
    t_half      = 5,       # h   | FDA label
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
    source      = "FDA label (Effexor XR 2012)"
  ),
  
  # ── Venlafaxine XR ─────────────────────────────────────────
  venlafaxine_XR = list(
    drug        = "Venlafaxine",
    class       = "SNRI",
    formulation = "XR",
    Ka          = 0.25,    # 1/h | Tmax ~5.5h | FDA label
    Vd          = 525,
    CL          = 91,
    F           = 0.45,
    t_half      = 5,
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
    source      = "FDA label (Effexor XR 2012)"
  )
)

# ── CYP2D6 genotype 보정 ───────────────────────────────────
cyp2d6 <- list(
  poor_metabolizer   = list(freq = 0.08, CL_factor = 0.30),
  normal_metabolizer = list(freq = 0.70, CL_factor = 1.00),
  rapid_metabolizer  = list(freq = 0.22, CL_factor = 1.60)
)

# ── 검증 요약 ──────────────────────────────────────────────
cat("✅ PK 파라미터 v4 로드 완료\n")
cat("약물:", paste(names(pk_params), collapse=" | "), "\n\n")
cat("[v4 변경사항]\n")
cat("  Paroxetine Vd: 8850→609 L (FDA label 기반)\n")
cat("  Paroxetine CL: 21.2→20.1 L/h (t½=21h 역산)\n")
cat("  Paroxetine Ka: 0.908→0.580 /h (Tmax 5.2h 맞게 조정)\n")
cat("  paroxetine_SR 제거: 중단 후 소실 패턴 IR과 동일\n")
cat("  Li 2022 V/F=8850 미적용: 정신과 환자 집단 + 비선형 PK 왜곡\n")