# 안드로겐 수용체(AR) 표적 PROTAC을 사례로 한 In-silico 예측 도구(분해 활성·투과성·삼원복합체 co-fold)의 비판적 평가

> *연구 질문: 현재 공개된 PROTAC 예측 도구는 학습에 없던 신규 화합물에 대해 "의미 있는 예측"이 가능한가?*

> 화학정보학 및 분자설계 기말 프로젝트 · 글로벌혁신신약학과 · 2026년 6월
> 본 보고서의 모든 정량 수치는 분석 스크립트(`scripts/`)와 그 출력(`outputs/`)을 단일 출처로 인용하며, 별도 하드코딩 값은 없다.

---

## 초록

PROTAC(proteolysis-targeting chimera) 신약 설계에는 분해 활성, 막 투과성, 삼원복합체 형성이라는 상충하는 다중 파라미터의 동시 최적화가 요구되며, 이는 인간의 직관만으로 다루기 어렵다. 이에 분해 활성 예측·투과성 예측·삼원복합체 구조 예측(co-fold)을 위한 여러 인공지능 도구가 공개되어 왔다. **본 보고서는 신약 설계 자체가 아니라, 이러한 공개 도구들이 실제로 "의미 있는 예측"을 할 만큼 성숙했는지를 비판적으로 검증하는 것을 목적으로 한다.** 임상적으로 검증되었고(GT20029 임상 2상) 공개 데이터가 풍부한 **안드로겐 수용체(AR; UniProt P10275)** 를 검증 사례(probe)로 삼아, 세 부류의 도구를 각 부류에 적합한 데이터에 대해 일관된 누수 통제(학습-미관측 전향 추론) 원칙 아래 평가하였다. 핵심 결론은 다음과 같다: **이 도구들은 자신이 학습한 분포(벤치마크) 안에서는 우수해 보이나, 학습에 포함되지 않은 신규 화합물에 대한 전향적(prospective) 예측에서는 성능이 크게 하락하며, 명목 순위조차 평가 설계에 따라 뒤집힐 만큼 불안정하다.** 따라서 현 단계의 도구는 단독 의사결정 근거가 아니라 실험 검증을 전제로 한 약한 사전확률(weak prior)로만 사용해야 한다.

---

## ① 타겟 소개와 연구 질문 (선정 이유)

### 1-1. 배경 — PROTAC 설계의 다중 파라미터 난제
안드로겐성 탈모(AGA)는 디하이드로테스토스테론(DHT)이 모낭의 **안드로겐 수용체(AR)** 와 결합해 모발을 소형화시키는 질환이다.¹ 기존 경구 치료제(피나스테리드 등)는 표적에 지속 결합하는 점유 구동(occupancy-driven) 방식이라 전신 부작용이 불가피하다.³ 이를 표적 단백질을 촉매적으로 **분해**하는 사건 구동(event-driven) PROTAC으로 대체하면 국소(topical) 적용에 유리하며, 실제 Kintor의 국소 AR-PROTAC **GT20029** 가 임상 2상에서 유효성과 전신 노출 최소화를 입증하였다.⁶ ⁸ 본 보고서의 직접적 동기는 2024년 *J. Med. Chem.* 논문(**jm4c02226**)¹¹의 AR-PROTAC 29개(A/B/C 시리즈)다.

그러나 새로운 국소 AR-PROTAC 설계는 화학적으로 난해하다. 충분한 **분해 효율**과 함께, 두피 각질층을 통과하되 전신으로는 빠져나가지 않는 **피부 체류성**을 동시에 만족해야 한다.¹⁰ 대부분의 PROTAC은 분자량 700~1200 Da·높은 TPSA로 리핀스키 5법칙을 벗어난 **bRo5** 공간에 존재해 PAMPA 등에서 저조한 투과성을 보인다.¹²

> **또한 PROTAC처럼 유연한 분자는 용매 환경에 따라 분자 내 수소 결합(IMHB)을 형성해 스스로 웅크리며 극성을 숨기는 '분자 카멜레온성(molecular chameleonicity)'을 띠는데, 이렇게 노출된 극성 표면적(EPSA)의 동적 궤적 등 수많은 변수가 복합 작용하여 투과성이 결정된다.¹⁵ 즉, 분자 구조의 미세한 조합 변화가 형태·투과성에 미치는 영향을 인간의 직관만으로 예측하는 것은 사실상 불가능하다.** 인간이 설계하기 까다로운 bRo5 화합물의 특성과, 낮은 용해도·막 내 갇힘으로 기존 체외(in vitro) 투과 분석법(PAMPA, Caco-2)이 빈번히 한계를 보이는 점을 고려할 때,¹⁴ 축적된 실험 데이터로 물성을 파악하고 인공지능으로 정량 예측하는 화학정보학적 접근의 의의는 점차 커지고 있다.

### 1-2. 연구 질문과 AR 선정 이유
이러한 기대에 비추어, 본 보고서는 다음을 묻는다 — **현재 공개된 PROTAC 예측 도구들이 실제로 신뢰할 만한 예측을 내놓는가?** 구체적으로 세 부류를 평가한다: (i) **분해 활성** 예측, (ii) **투과성** 예측, (iii) **삼원복합체(co-fold) 구조** 예측.

검증 사례로 **AR**을 택한 이유는 두 가지다. 첫째, AR은 GT20029·jm4c02226로 **임상·전임상 근거가 분명한 PROTAC 타겟**이다. 둘째, AR은 PROTAC-DB에 다수 수록된 **데이터가 풍부한 타겟**이어서, 도구의 입력 요건을 충족시키며 그 성능을 실측값과 대조하기에 적합하다. 즉 본 보고서에서 AR은 치료제 설계의 대상이라기보다, **도구의 성숙도를 시험하는 잘 특성화된 시금석**이다.

---

## ② 데이터·구조 확보 및 전처리

### 2-1. 데이터 자산 — 전체 규모와 AR 부분 ([표1])

**[표1] 데이터셋 전체 규모 vs AR(P10275)**

| 데이터셋 | 전체 행 | 전체 unique 화합물(ik14) | AR(P10275) 행 | AR unique ik14 |
|---|---|---|---|---|
| 구 PROTAC-DB (`ProtacDB_MTL_CLS.csv`) | 6,107 | 6,077 | 367 | 360 |
| 신 PROTAC-DB 3.0 (`protac.xlsx`)¹⁸ | 15,502 | 10,383 | 777 | 658 |
| case study (`jm4c02226`)¹¹ | 29 | 29 | 29 (전부 AR) | 29 |

**[표1b] 실험 측정값의 데이터 개수 (전체 / AR).** 구조에서 계산되는 값(분자량·XLogP3·TPSA·HBA/HBD·회전결합·고리수·InChIKey)은 제외하고, **실험으로만 얻는 값**의 비결측 개수를 집계하였다.

| 실험 측정값 | 데이터셋 | 전체 | AR |
|---|---|---|---|
| **DC50 (nM)** · **Dmax (%)** (분해) | 신DB | 3,592 / 2,729 | 167 / 139 |
| Percent degradation (%) | 신DB | 2,774 | 252 |
| IC50 (→Target / →E3) | 신DB | 1,831 / 453 | 9 / 0 |
| 세포 활성 IC50 / EC50 | 신DB | 4,405 / 388 | 191 / 22 |
| Caco-2 A2B / B2A / PAMPA Papp | 신DB | 91 / 69 / 32 | 19 / 19 / 10 |
| pDC50 / Dmax / Degrader_Class (실측 관측) | 구DB | 1,082 / 574 / 507 | 62 / 52 / 51 |
| Skin Retention / DC50 / Dmax | case study | 29 / 13 / 13 | (전부 AR) |

전체 규모는 크지만 **모델의 핵심 학습·평가 신호인 분해 라벨(DC50/Dmax)은 희소**하다. 신DB에서 DC50은 전체의 23%, AR의 21%에만 존재한다. "데이터가 많다"는 인상의 상당 부분은 계산 descriptor와 결합·세포 활성이며, 분해 라벨 자체는 모델 학습·검증에 빠듯하다 — 이는 ⑤에서 다룰 도구의 일반화 실패와 직결되는 구조적 제약이다.

### 2-2. 전처리 규칙
- **AR 필터는 `Uniprot == 'P10275'` 단일 기준.** `Target` 컬럼 정규식('AR'/'Androgen')은 AR-V7·T878A 등 변이체를 흡수해 1,193행까지 부풀려지므로 사용하지 않았다.
- **화합물 식별**은 두 DB의 precomputed InChIKey가 RDKit과 불일치하므로, 양쪽 모두 RDKit `MolToInchiKey` **앞 14자(골격키, ik14)** 로 재계산해 통일하였다(stereo 무관 골격 기준).
- **case study 파싱**: cp949 인코딩·다중헤더, '±'가 '÷'로 mojibake되어 `÷→±` 치환 후 선행 float만 취하였다. 3번째 칸 헤더는 `Skin Retention Rate(%)`(A·B·C 29/29 보유), DC50/Dmax는 B/C 13개만 보유(A 16개는 활성 미측정). 수치형 DC50은 5개(B3 8.97, B6 207.03, C1 103.87, C5 70.85, C6 199.5), 나머지는 `>1000` censored.
- **E3 리가아제 정체성**: 29개 SMILES를 직접 검사한 결과 글루타리미드(CRBN cue) 29/29, 하이드록시프롤린(VHL cue) 0/29 → 전 계열이 **CRBN(Q96SW2)** 기반이다. B/C의 클로로-시아노-페녹시-사이클로헥실 모티프는 VHL 리간드가 아니라 원논문 표기상 **enzalutamide analog 계열의 AR 길항제 워헤드**다.

### 2-3. 학습/검증 분할과 그 시간성
도구의 일반화를 정직하게 측정하려면 **학습에 포함되지 않은 화합물**로 평가해야 한다. 분할 기준은 골격키(ik14) 멤버십 — 즉 평가 도구가 학습 때 본 화합물인지 여부 — 로 정의하였다(엄밀한 캘린더 분할은 구DB에 연도 컬럼이 없어 불가). 다만 이 멤버십 분할이 시간성과 부합하는지를 사후 검증하였다: `protac.xlsx`의 Article DOI(unique 766개)를 **Crossref로 99% 게재연도 해소**한 뒤 비교한 결과, 학습-미관측(hold-out) 화합물의 게재연도 중앙값은 **2024년**(IQR 2023–2024)인 반면 기존 DB는 **2022년**(IQR 2021–2024)으로, hold-out이 체계적으로 더 최신이었다([그림 pubdate] `outputs/pubdate_kde.png`). 즉 멤버십 분할은 사실상 시간 분할과 대체로 일치한다(다만 분포가 겹쳐 하드 날짜 컷오프는 아니다).

**[표1c] AI 검증 분할 — 화합물·AR 수**

| 구분 | 화합물(unique) | AR(P10275) |
|---|---|---|
| **검증 hold-out** (STAN 미관측, AI 평가용) | 1,134 | 42 |
| **학습-관측 대조군** (leaked; 도구가 학습 때 본 화합물 → 암기 측정용 기준선) | 281 | 63 |
| 나머지 (PROTAC-DB 3.0 − hold-out) | 9,249 | 616 |
| (참고) PROTAC-DB 3.0 전체 | 10,383 | 658 |

**[표1d] 활성값(측정 항목)별 데이터 수 — hold-out vs 나머지** (unique 화합물 비결측 수, 전체 / AR)

| 활성값 | hold-out (전체 / AR) | 나머지 (전체 / AR) |
|---|---|---|
| **DC50 (nM)** | 1,134 / 42 | 1,590 / 94 |
| **Dmax (%)** | 1,134 / 42 | 921 / 71 |
| Percent degradation (%) | 99 / 1 | 2,079 / 190 |
| IC50 (→Target / →E3 결합) | 82 / 0 · 42 / 0 | 1,283 / 9 · 229 / 0 |
| 세포 활성 IC50 | 325 / 14 | 2,920 / 147 |
| Caco-2 A2B Papp | 26 / 0 | 36 / 7 |

hold-out 1,134개는 모델 입력 요건상 **DC50·Dmax를 둘 다 보유**한다(그래서 둘 다 1,134). "나머지"는 DC50(1,590) > Dmax(921)로, **Dmax 없이 DC50만 있어 hold-out 요건(둘 다)을 못 채운 화합물**이 많다. 분해 라벨(DC50/Dmax)은 결합·세포활성보다 희소하다(②-1).

### 2-4. 구조(3D) 확보의 한계
AR-LBD의 PROTAC 삼원복합체(AR+E3+PROTAC) 실험 구조는 부재하다. 이는 구조 기반 도구의 적용 가능성을 좌우하는 핵심 제약이며, ③·④-(다)·⑤에서 다룬다.

---

## ③ 적용한 방법과 선택 이유

세 부류의 도구를 평가하였으며, 모든 평가는 **사전학습 가중치 동결 후 학습-미관측 화합물에 대한 전향적 추론**을 원칙으로 한다.

| 부류 | 도구·방법 | 평가 설계 |
|---|---|---|
| 화학공간 | Morgan FP(r=2, 2048bit) + PCA, MW–TPSA overlay | AR PROTAC의 bRo5 위치 확인 |
| **(가) 분해 활성** | **PROTAC-STAN**(GNN+ESM) — 주 평가; **Ribes PROTAC-Degradation-Predictor**(포켓-프리 FP+XGBoost/MLP) — 교차검증; **DegradeMaster**(구조기반) — 참조 인용 | 전체 타겟 prospective 분할 + 도구 간 비교 |
| **(나) 투과성** | **Potts-Guy** QSPR(피부 투과) 직접 구현; **PROTAC-TS**(TabPFN, Caco-2 세포막 투과) | case study 29개 예측, 실측 잔류·계산 proxy와 대조 |
| **(다) 삼원복합체 co-fold** | **Boltz-2** | 리드 화합물(B3) AR-LBD+CRBN+PROTAC 삼원복합체 예측 |
| 신규 설계 | SAR 기반 수동 medchem 편집 | 도구를 weak prior로만 사용(순환논증 회피) |

**(가) 분해 활성 도구.** DeepPROTACs·DegradeMaster·AiPROTAC 등 구조 기반 예측기는 화합물마다 단백질 포켓과 3원 도킹 포즈를 필수로 요구하나 AR 삼원구조가 없어 직접 적용이 불가능하다. 따라서 SMILES·서열만으로 구동되는 **포켓-프리 도구(PROTAC-STAN, Ribes)** 를 주 평가 대상으로 삼고, 구조 기반 DegradeMaster는 저자 보고값으로만 참조하였다. 타겟 수준 degradability만 예측해 화합물 단위 랭킹이 불가능한 DegradoMap²⁴은 평가에서 제외하였다. PROTAC-STAN은 2-class 이진 분류기이므로 연속 DC50 회귀가 아닌 **활성확률 랭킹·순위상관**으로 평가하며, 양성 정의는 도구 자기 기준(`DC50<100 nM AND Dmax≥80%`)을 따른다(모델군마다 AND/OR 정의가 달라 통일 컷을 강요하지 않음). 추론 구현 세부(torch 호환 패치·class-1 확률 덤프 러너)는 부록에 정리한다.

**(나) 투과성 도구.** 피부 투과는 Potts-Guy QSPR `logKp(cm/s)=0.71·logP−0.0061·MW−6.3`¹⁹를 RDKit으로 직접 구현하였다(출처 불명의 "Mitragotri 식"은 임의 계수로 확인되어 철회). 세포막 투과는 **PROTAC-TS**의 Caco-2 예측 모델(PROTAC-DB Caco-2 약 90개로 학습한 TabPFN, count-Morgan 500차원)을 실제 구동해 29개를 예측하였다. 두 도구가 측정하는 **세포막 투과 ≠ 피부 투과 ≠ 피부 잔류**라는 도메인 구분을 데이터로 확인하는 것이 핵심 목적이다.

**(다) Co-fold 도구.** **Boltz-2**(AlphaFold3 계열 확산 모델)로 리드 B3을 포함해 **case study 29개 전부**의 삼원복합체를 예측하였다. 입력은 **AR-LBD(P10275 669–919, 251aa)** + **CRBN(Q96SW2, 442aa)** + 각 화합물 리간드이며(29개 모두 CRBN 기반), 전장 AR(920aa)은 대부분 비정형 N-말단이라 LBD로 한정하였다(워헤드가 LBD에 결합). 서열은 `ProtacDB_MTL_CLS`에서 취했고, 동일한 AR-LBD·CRBN MSA는 1회 계산 후 29개에 재사용하였다.

---

## ④ 결과

### 4-0. 화학공간 — AR PROTAC은 명백한 bRo5 ([그림1] `outputs/fig1_chemspace.png`)
세 데이터를 Morgan FP–PCA로 투영(PC1 20.8%, PC2 9.5%, 합 30.3%)하면, B/C 시리즈는 MW 725~792·TPSA 144~155로 Ro5 선(MW=500, TPSA=140)을 명백히 위반하고 bRo5 상한(MW≈1000, TPSA≈250; Doak 2014)²² 영역에 위치한다. A 시리즈(MW 466~597)는 상대적으로 작다. 전체 PROTAC-DB와 AR 부분의 인터랙티브 화학공간은 `chemspace_view1_AR.html`로 함께 제공한다.

### 4-1. (가) 분해 활성 도구 — 전체 타겟 전향적 검증 ([표2]·[표3]·[그림2])
AR 한정 평가는 표본이 작아 도구의 일반 성능을 말할 수 없으므로, 평가를 **전체 타겟으로 확대**하였다. PROTAC-DB 3.0에서 PROTAC-STAN이 학습 때 보지 못한 **prospective 화합물 1,134개(51개 타겟, 양성 388/음성 746)** 를 주 test로 삼고, 여러 in-distribution 기준과 대조하였다.

**[표2] PROTAC-STAN 분해예측 성능 — 평가 기반별** (`outputs/timesplit_metrics.csv`·`fig_testset_vs_holdout.png`; bootstrap 2,000회 CI, strict GT)

| 평가 기반 | n | AUROC [95% CI] | 성격 |
|---|---|---|---|
| 학습셋 재채점 | 1,296 | 0.971 | 암기 상한 |
| STAN 자체 test set | 207 | 0.910 [0.867, 0.949] | in-distribution (STAN 공개 벤치마크) |
| leaked (DB3.0 ∩ STAN학습) | 281 | 0.909 [0.873, 0.941] | in-distribution 대조 |
| **전향적 hold-out (전체)** | **1,134** | **0.574** [0.539, 0.607] | **OOD — 주 결과** |
| 전향적 hold-out (Ribes 공통 부분집합) | 786 | 0.550 [0.510, 0.588] | OOD, STAN↔Ribes 비교 기반 |
| 전향적 hold-out — AR(P10275) | 42 | 0.702 [0.524, 0.854] | AGA 타겟 |

- **in-distribution(자체 test 0.910·leaked 0.909·train 0.971)이 전향적 hold-out 0.574(MCC 0.098)로 붕괴**한다([그림2-1] `ts_roc.png`·[그림2-4] `fig_testset_vs_holdout.png`). in-distribution 기준을 STAN 공식 test로 잡든 leaked로 잡든 결론은 같다(둘 다 ≈0.91).
- STAN 성능은 **1,134(0.574)·786(0.550)에서 일관**되게 무작위 부근이다(786 = Ribes와 공통 예측 가능한 부분집합, 아래 [표3]).
- 타겟별로는 AUROC 중앙값 0.448·ERα(n=258) 0.276이고 양성-음성 쌍의 약 93%가 타겟 간(cross-target) 쌍이다([그림2-2] `ts_per_target.png`). AR은 0.702(n=42)지만 표본이 작다. (해석 ⑤-1)

**STAN vs Ribes — 공통 786개.** STAN은 서열(ESM)만으로 1,134 전부 예측하지만 **Ribes는 세포주·E3·타겟이 모두 자기 학습 어휘에 있어야** 해 786개만 예측 가능하다(신규 화합물의 ~31%는 적용 도메인 밖이라 예측 불가 — 그 자체가 한계). 공정 비교를 위해 두 도구가 **공통 예측한 786개**에 **두 활성 임계값을 양쪽 모두**에 적용하였다([표3], [그림2-3] `fig_multitool.png`·[그림2-6] `fig_threshold_x_model.png`).

**[표3] 동일 786개에서 STAN vs Ribes — 임계값·렌즈별 AUROC** (`outputs/multitool_metrics.csv`·`threshold_x_model.csv`)

| 모델 | strict GT(양성29%) | native GT(양성70%) | KRAS 제외(strict) | 타겟별 macro |
|---|---|---|---|---|
| PROTAC-STAN | 0.550 [0.51,0.59] | 0.559 [0.51,0.60] | 0.427 | 0.531 |
| Ribes-standard (random-CV) | **0.696** [0.65,0.74] | 0.491 [0.45,0.53] | 0.618 | 0.593 |
| Ribes-target (group-CV) | 0.685 [0.64,0.73] | **0.577** [0.54,0.62] | 0.670 | 0.501 |

*strict = `DC50<100nM & Dmax≥80%`(양성 29%), native = `DC50≤1000nM & Dmax≥60%`(양성 70%); 두 기준에서 라벨이 뒤바뀌는 화합물 410개(36%, 전부 DC50 100–1000nM·Dmax 60–80% 중간대).*

- **순위가 평가 설계에 따라 뒤집힌다**: strict에선 Ribes-standard(0.696)가 1위지만 native에선 0.491로 STAN보다 낮은 꼴찌가 된다. KRAS 단일클래스 블록(786 중 214개 거의 전부 음성) 제거 시 STAN 0.427, 타겟 구성을 통제한 타겟별 macro는 셋 다 0.50–0.59다.
- **STAN은 임계값에 둔감**(0.55~0.56)하나 Ribes-standard는 0.49~0.70으로 가장 크게 출렁인다. 두 도구의 일치도는 활성확률 Spearman 0.19·Cohen κ=0.13으로 낮다. **그럼에도 모든 칸이 무작위 부근**이라 "전향적 미성숙" 결론은 임계값·렌즈와 무관하다. (해석 ⑤-2)

**DegradeMaster 및 도구 공통 패턴.** 구조 기반 DegradeMaster는 화합물마다 3원 pocket이 필요해 본 세트에 적용할 수 없어 저자 보고값만 인용한다(supervised PROTAC-1K 0.854·semisupervised PROTAC-8K random 0.882)²⁸; 이 셋을 단일 held-out으로 공정 비교할 수도 없다(PROTAC-8K가 STAN 학습셋과 **100%**·Ribes와 **68%** 중첩). 각 도구의 *자체* 평가만 봐도 패턴은 같다([표3b]·[그림2-5] `fig_alltools_selftest.png`).

**[표3b] 도구별 자체 평가: random split vs 타겟 교차** (각 도구 자체 리포트 기준)

| 도구 | random / in-dist | 타겟 교차(target-split) | similarity split |
|---|---|---|---|
| PROTAC-STAN | 0.910 (자체 test) | 0.531 (per-target macro) | — |
| Ribes (XGBoost) | 0.906 | **0.585** | 0.897 |
| DegradeMaster | 0.854·0.882 | 측정 불가(구조 필요) | — |

→ Ribes는 similarity split(0.90)은 유지하나 **target-split(0.58)만 무너진다 — 새 골격엔 강하나 새 타겟엔 약함**(STAN per-target 붕괴와 동일 메커니즘). 도구·모델·split·임계값을 달리해도 결론은 견고하다: **in-distribution ≈0.85~0.91 → 타겟 교차/전향적 ≈0.45~0.59**. "벤치마크 성능"은 일관되게 실사용(신규 타겟·신규 골격) 성능을 과대평가한다.

### 4-2. (나) 투과성 도구 — 3중 도메인 구분의 실증 ([그림3]·[그림3b])
- **Potts-Guy(피부 투과)**: 29개 logKp는 −8.22 ~ −6.74로 압축되었고, MW>750이 11/29·TPSA>140이 13/29로 전체의 44.8%가 적용도메인(AD) 밖이다. MW항 지배로 모든 PROTAC이 균일하게 낮은 logKp로 압축되어 변별력을 잃는다(예상된 결과). 절대 Kp는 정량 신뢰가 불가하다. 예측 logKp와 실측 피부 잔류의 Spearman은 0.042로, 피부 잔류(depot)와 피부 통과(Kp)가 물리적으로 다른 endpoint임을 보여준다([그림3] `outputs/fig3_separation.png`).
- **PROTAC-TS(세포막 투과)**: PROTAC-DB Caco-2 약 90개로 TabPFN을 직접 학습(자체 LOOCV R²=0.78, `~/PROTAC-TS` make_model)해 29개를 예측하였다(누수 없음). 예측 세포막 투과는 실측 피부 잔류와 약한 양의 상관(Spearman +0.41), Potts-Guy 피부 logKp와는 무상관(−0.05)을 보였다([그림3b] `outputs/fig_protacts.png`).

→ 세 endpoint 간 상관(요약): Caco-2 예측 vs 피부 잔류 **+0.41**(약한 양), Potts-Guy logKp vs 피부 잔류 **0.042**, Caco-2 예측 vs Potts-Guy logKp **−0.05**. 두 예측 모두 좁은 범위로 압축됨(변별력 제한). 카멜레온성은 2D proxy(TPSA vs 회전결합)로만 근사했고 3D EPSA 정량(min-3D-PSA·ΔPSA·IMHB 등)은 거대 PROTAC의 형태 샘플링 비용 문제로 다음 단계로 둔다. (도메인 비환원성 해석은 ⑤-3)

### 4-3. (다) Co-fold 도구 — Boltz-2 삼원복합체 (B3 및 case study 전 29개) ([표5]·[그림4])
리드 화합물 **B3(TJA-107)** 을 포함해 **case study 29개 전부**의 AR-LBD+CRBN+PROTAC 삼원복합체를 Boltz-2로 예측하였다(**29/29 성공**; AR-LBD·CRBN의 MSA는 1회 계산 후 전 화합물에 재사용).

**[표5] Boltz-2 신뢰도 — 리드 B3 및 계열 평균** (`outputs/boltz_all29_confidence.csv`)

| 대상 | confidence | ipTM | ligand_ipTM | protein_ipTM | complex pLDDT |
|---|---|---|---|---|---|
| **B3** (리드) | 0.819 | 0.715 | **0.937** | 0.671 | 0.846 |
| A 계열(16) | 0.777 | 0.580 | 0.892 | **0.426** | 0.826 |
| B 계열(6) | 0.820 | 0.719 | 0.935 | 0.679 | 0.845 |
| C 계열(7) | 0.812 | 0.700 | 0.908 | 0.666 | 0.842 |

- **리간드 포즈 신뢰도(ligand_ipTM)는 29개 전부 0.84~0.96으로 일관되게 높다** — Boltz가 PROTAC 워헤드·E3 리간드를 포켓에 자신 있게 배치한다.
- 반면 삼원복합체의 핵심인 **단백질-단백질 계면(AR-LBD:CRBN) 신뢰도(protein_ipTM)는 중간이고 편차가 크다**(A 0.43 vs B/C 0.67~0.68). 특히 enzalutamide-analog 워헤드가 아닌 A 계열이 체계적으로 낮다([그림4] `outputs/fig_boltz_all29.png`).
- **예측 신뢰도는 실측 분해능과 무의미한 상관**이다(수치형 DC50 5개에서 Spearman(ipTM, −DC50)=+0.30, p=0.62). 최강 활성 B3(8.97 nM)와 중간 활성들의 ipTM이 비슷해, 그럴듯한 co-fold가 곧 활성 분해제를 뜻하지 않는다. (해석은 ⑤-4)

### 4-4. case study(AGA) 적용과 신규 분자 가설
case study 29개에 PROTAC-STAN을 적용하면 leakage-free 활성 앵커(C5 70.85 nM, C1 103.87 nM, C6 199.5 nM)가 모두 상위에 랭크된다. 다만 누수 화합물 B3·B5의 확률은 오히려 낮았는데(0.25, 0.16), B5는 실측 비활성이라 정답 예측이고 활성인데 낮은 사례는 B3 한 건뿐이어서, 누수가 점수를 단순 부풀리지 않았음을 확인할 수 있다(소표본이라 정성 관찰). 한편 도구들의 한계가 드러난 만큼 신규 분자 3종(D1~D3; 부모 C1/C5/C6)은 **모델 점수가 아니라 SAR 규칙**(CRBN 글루타리미드 + enzalutamide-analog 워헤드 약리단 보존)에 근거해 제안하고, STAN 점수는 weak prior로만, Tanimoto 최근접 거리(0.63~0.70)와 "실험 검증 필요"를 동반해 제시한다(`outputs/table4.csv`).

---

## ⑤ 결과에 대한 비판적 고찰 (가장 중요) — 도구는 충분히 성숙했는가?

본 보고서의 핵심 질문에 대한 답은 **"아직 아니다"** 이며, 도구 부류별 판정을 먼저 요약하면 다음과 같다([표6]).

**[표6] 도구 부류별 성숙도 판정 요약**

| 도구 부류 | 보고/in-distribution | 본 검증(전향적·AR 사례) | 판정 |
|---|---|---|---|
| **분해 활성** | AUROC 0.85~0.91 | pooled 0.57·**타겟 내 ≈0.45(무작위 부근)**, 도구 간 순위 불안정·κ=0.13 | **미성숙** — 화합물 단위 변별 거의 불가 |
| **투과성** | (소분자/펩타이드 검증) | Potts-Guy 44.8% AD 밖·압축, Caco-2 예측도 압축; 두 지표 무상관 | **제한적** — 변별력 낮고 도메인 비환원 |
| **Co-fold** | (AF3급 일반 성능) | ligand_ipTM 0.937 / **protein_ipTM 0.671**(단일 시드) | **부분적** — 리간드 배치는 ○, 단백질 계면 불확실 |

세부 근거는 다음과 같다.

**5-1. 분해 예측 — 벤치마크 성능과 실사용 성능의 괴리.** PROTAC-STAN은 자신이 학습한 화합물에서는 AUROC 0.909로 거의 완벽하지만, 학습에 포함되지 않은 화합물에서는 0.574(MCC 0.098)로 무작위에 근접하게 붕괴한다. 즉 **보고되는 높은 성능의 상당 부분은 일반화가 아니라 암기(memorization)** 다. 더 결정적으로 pooled 0.574조차 부풀려진 값이다 — AUROC가 비교하는 양성-음성 쌍의 약 93%가 서로 다른 타겟 간(cross-target) 쌍이기 때문이다. 타겟 구성을 통제해 타겟별로 보면 AUROC 중앙값 0.448·n-가중 0.45로 무작위 부근이고, 데이터가 가장 많은 ERα(n=258)에서는 오히려 0.276이다. 다시 말해 모델은 "어느 타겟이 분해가 잘 되는 편인가"(타겟 수준 prior)는 부분적으로 반영하지만, 실무에서 정작 필요한 "한 타겟 안에서 어떤 화합물이 더 활성인가"는 거의 구분하지 못한다. 이 한계는 단일 타겟·소표본 평가로는 결코 드러나지 않으며, "전체 타겟 + 누수 대조 + 타겟별 분해"라는 평가 설계를 통해서만 드러난다 — 이것이 본 검증의 방법론적 핵심이다.

**5-2. 순위의 불안정성과 도구 간 불일치.** 방법론이 전혀 다른 두 포켓-프리 도구(STAN: GNN+ESM, Ribes: FP+트리/MLP)를 동일 held-out에 적용했을 때, "어느 도구가 우수한가"라는 결론 자체가 **평가 설계에 따라 뒤집힌다**: pooled-STAN라벨에서는 Ribes(0.696)>STAN(0.550)이지만, Ribes를 그 자신이 학습한 라벨로 평가하면 0.491(무작위)로 떨어지고, 단일클래스 KRAS 블록을 제거하면 STAN이 0.427(무작위 이하)로 떨어지며, 타겟 구성을 통제한 macro에서는 셋 다 0.50~0.59로 수렴한다. 더욱이 두 도구의 예측은 서로 크게 불일치한다(활성확률 Spearman 0.19, Cohen κ=0.13 — 우연 수준을 겨우 넘는다). 따라서 **어느 단일 도구의 단일 숫자도 신약 설계의 독립 근거로 신뢰하기 어렵다.** 한편 논문들이 보고하는 0.85~0.91 AUROC(DegradeMaster PROTAC-1K 0.854·PROTAC-8K 0.882 포함)는 대부분 train/test가 chemotype을 공유하는 in-distribution split의 값으로, STAN의 학습-관측 0.909와 같은 낙관 영역이다. 게다가 이 도구들을 동일 held-out으로 공정 비교하는 것조차 불가능하다 — 구조 기반 도구가 구조를 가진 유일한 셋(PROTAC-8K)이 STAN 학습셋과 100%·Ribes와 68% 겹쳐, **분야 벤치마크 자체에 학습-평가 중첩(contamination)이 내재**하기 때문이다.

**5-3. 투과 예측 — "투과성"은 단일 측정으로 환원되지 않는다.** Potts-Guy는 소분자 회귀식이라 bRo5 PROTAC의 44.8%가 적용도메인 밖이고, MW항 지배로 모든 화합물을 균일하게 낮은 logKp로 압축해 변별력을 잃는다. PROTAC-TS는 PROTAC으로 학습되어 도메인 안쪽이지만 예측값이 좁게 압축되어 역시 변별력이 제한된다. 세 endpoint 간 상관을 보면, **세포막 투과(Caco-2 예측)와 피부 잔류는 약한 양의 상관(+0.41)으로 일부 공통 의존(친유성)을 공유**하지만, **두 in-silico 투과 지표(Caco-2 예측 vs Potts-Guy logKp)는 서로 무관(−0.05)** 하고 피부 통과(Kp) vs 잔류도 무상관(0.042)이다. 즉 세 endpoint가 동일량으로 환원되지 않으며, 특히 두 계산 투과 지표가 서로를 예측하지 못한다는 점이 "permeability를 하나의 숫자로 다루는 것"의 위험을 보여준다. 더욱이 어떤 투과 모델도 모낭(follicular) 경로·유한용량·제형 효과를 반영하지 못하는데, 이는 두피 국소제의 핵심 메커니즘이다.

**5-4. Co-fold — 그럴듯한 리간드 포즈, 불확실한 계면, 활성과 무관.** case study 29개 전부의 삼원복합체를 예측한 결과, 리간드 포즈 신뢰도는 일관되게 높았으나(ligand_ipTM 0.84~0.96) 삼원복합체 형성의 핵심인 단백질-단백질 계면(protein_ipTM)은 중간 수준이고 계열별 편차가 컸다(A 계열 0.43 vs B/C 0.67~0.68). 더 중요하게, **Boltz 신뢰도는 실측 분해 활성과 무의미한 상관**(ipTM vs −DC50 Spearman +0.30, n=5, p=0.62)이어서, 신뢰도가 높은 co-fold가 곧 활성 분해제를 의미하지 않는다. 모두 단일 시드 예측이므로 pLDDT·ipTM만으로 결합·활성을 단정해서는 안 된다.

**5-5. 공통 근본 원인 — 데이터 희소성.** 위 한계의 상당 부분은 분해 라벨의 희소성(②)에서 비롯된다. AR 관측 양성은 25개, 신DB DC50도 전체의 23%뿐이다. 데이터가 빈약한 영역에서는 복잡한 모델일수록 일반화보다 암기로 기우는데, 단순한 Ribes(FP+트리)가 prospective에서 무거운 GNN인 STAN과 비등하거나 앞섰던 것이 그 방증이다. 즉 현 단계의 병목은 "더 큰 모델"이 아니라 "더 많고 깨끗하며 누수 없는 데이터"다.

**5-6. 방법론적 메타-교훈.** 본 보고서의 가장 일반적인 결론은 **도구가 보고한 숫자를 그대로 믿어서는 안 된다**는 것이다. 누수 통제·타겟 통제·라벨 정합·다중 도구 교차검증 중 하나라도 빠지면, 화려한 벤치마크 점수가 실사용 성능을 크게 과대평가한다. 현 단계의 PROTAC 예측 도구(분해·투과·co-fold)는 의미 있는 *단독* 예측을 하기에는 아직 이르며, 실험 검증을 전제로 한 **약한 사전확률(weak prior)** 로, 가급적 **여러 도구의 합의(consensus)** 형태로만 사용하는 것이 정당하다.

---

## ⑥ 한계와 다음 단계

**한계.**
- 분해 도구의 전향적 성능이 약하고(미관측 AUROC≈0.57, 타겟 내 변별 ≈무작위) 순위가 불안정하다.
- 도구 적용 범위가 제한적이다 — STAN은 ESM 사전(283 Uniprot)이 커버하는 타겟·E3로, 구조 기반 도구는 삼원구조 부재로 제약된다.
- 분할은 골격 멤버십 기준이며(게재연도와 대체로 일치하나 하드 컷오프는 아님), 엄밀한 temporal generalization을 주장하지 않는다.
- co-fold는 단일 시드·1 샘플 예측이다.

**다음 단계.**
1. Boltz-2를 다중 시드·`diffusion_samples` 앙상블로 재실행해 삼원복합체 신뢰구간을 정량하고, 예측 구조를 구조 기반 분해 예측기 입력으로 활용한다.
2. PROTAC-STAN을 scaffold/temporal split로 재학습·재보정하고, **타겟 내 랭킹**을 주 평가축으로 삼아 일반화를 높인다.
3. 3D EPSA·ΔPSA·IMHB 기반 카멜레온성 정량(min-3D-PSA, 이중 유전율 형태 샘플링)으로 투과 proxy를 보강한다.
4. 신규 3종(D1~D3)의 합성·HDPC AR 분해 assay·피부 잔류 측정으로 weak prior를 실험 검증한다.
5. 모낭 표적 전달(제형·나노입자)을 고려한 전신/국소 분리 모델링.

---

## 부록 — 재현 정보
- **환경**: 메인 분석 `~/anaconda3/envs/protac`(Python 3.10, torch 2.8.0+cu128, RDKit 2026.03.1). 투과(PROTAC-TS) `protac_ts`(TabPFN). co-fold(Boltz-2) `protac_boltz`. 도구 교차검증(Ribes) `protac_ribes`.
- **스크립트(`scripts/`)**: `phase0_freeze.py`, `fig1_chemspace.py`, `run_stan_inference.py`/`run_stan_ts.py`(STAN 추론), `timesplit_build.py`/`timesplit_eval.py`(전향적 평가), `multitool_eval.py`(STAN vs Ribes), `phase3_skin.py`(Potts-Guy), `protacts_predict.py`(PROTAC-TS Caco-2), `phase5_design.py`(신규 설계), `resolve_doi_years.py`/`build_pubdate_kde.py`(게재연도), `build_chemspace_html.py`(인터랙티브 화학공간), `stan_patch.py`, `stan_testset_build.py`(STAN 자체 train/test 평가), `boltz_all29_build.py`/`boltz_all29_collect.py`(co-fold 29개 생성·집계). Boltz 입력 `~/PROTAC_MTL_v5/boltz_{B3_ternary,all29}/`. Ribes 자체 test 메트릭은 `~/PROTAC-Degradation-Predictor/reports/`.
- **산출물(`outputs/`)**: `table1~5.csv`, `timesplit_metrics.csv`·`timesplit_per_target.csv`·`multitool_metrics.csv`, 그림 `fig1_chemspace.png`·`ts_roc.png`·`ts_per_target.png`·`fig_multitool.png`·`fig3_separation.png`·`fig_protacts.png`·`pubdate_kde.png`·`fig_boltz_all29.png`·`fig_testset_vs_holdout.png`·`fig_alltools_selftest.png`·`fig_threshold_x_model.png`, `boltz_all29_confidence.csv`·`threshold_x_model.csv`, 인터랙티브 `chemspace_view1_AR.html`·`chemspace_view2_split.html`, Boltz 구조 29개 `~/PROTAC_MTL_v5/boltz_all29/out/.../predictions/<ID>/<ID>_model_0.pdb`(+ B3 별도).

---

## 참고 자료

1. Androgenetic Alopecia — StatPearls (NIH). https://www.ncbi.nlm.nih.gov/books/NBK430924/
3. The diagnosis and treatment of androgenetic alopecia — Forum Dermatologicum. https://journals.viamedica.pl/forum_dermatologicum/article/view/95391/76623
6. GT20029 Topical Option for AGA — Bauman Medical. https://www.baumanmedical.com/gt20029-hair-loss-treatment-phase-2-results/
8. Kintor Pharma AGA Trial Reaches Primary Endpoint — Clival. https://clival.com/news/kintor-pharma-announces-its-trial-for-androgenetic-alopecia-reaches-primary-endpoint
10. Expanding the Scope of PROTACs: Topical Delivery — *J. Med. Chem.* 2025 (PMC12670424). https://pmc.ncbi.nlm.nih.gov/articles/PMC12670424/
11. **Discovery of a Novel Non-invasive AR PROTAC Degrader for the Topical Treatment of AGA.** *J. Med. Chem.* **2024, 67(24), 22218–22244** (online 2024-12-06; PubMed 39641607). DOI 10.1021/acs.jmedchem.4c02226
12. From Ro5 to bRo5: Refinement of Physicochemical Properties — PMC9511483. https://pmc.ncbi.nlm.nih.gov/articles/PMC9511483/
14. Systematic Investigation of the Permeability of AR PROTACs — PMC7429968. https://pmc.ncbi.nlm.nih.gov/articles/PMC7429968/
15. A rational control of molecular properties in the bRo5 chemical space (Garcia-Jimenez thesis). https://iris.unito.it/retrieve/1ba8b290-4ce5-448f-8fff-833f9a0d754c/Thesis_Diego_Garcia_Jimenez_final_version.pdf
18. **Ge J. et al. PROTAC-DB 3.0.** *Nucleic Acids Research* **2025;53(D1):D1510–D1515.** DOI 10.1093/nar/gkae768
19. **Potts RO, Guy RH.** Predicting skin permeability. *Pharm. Res.* **1992, 9(5), 663–669.**
22. **Doak BC, Over B, Giordanetto F, Kihlberg J.** Oral druggable space beyond the rule of 5. *Chem. Biol.* **2014, 21(9), 1115–1142.**
23. **PROTAC-STAN** — GNN+ESM 기반 PROTAC 분해 분류 모델. https://github.com/PROTACs/PROTAC-STAN
24. **DegradoMap** — 타겟 수준 degradability 예측(화합물 랭킹 불가). https://github.com/bryanc5864/DegradoMap
25. 구조 기반 분해 예측기(AR 삼원구조 부재로 직접 적용 불가): DeepPROTACs (https://github.com/fenglei104/DeepPROTACs), AiPROTAC (https://github.com/LiZhang30/AiPROTAC).
26. **Boltz-2** (co-folding, 삼원복합체 예측). https://github.com/jwohlwend/boltz
27. **Ribes S. et al.** Modeling PROTAC degradation activity with machine learning. *Artificial Intelligence in the Life Sciences* 2024 — 포켓-프리 예측기 PROTAC-Degradation-Predictor. https://github.com/ribesstefano/PROTAC-Degradation-Predictor
28. **Liu J. et al.** Accurate PROTAC targeted degradation prediction with DegradeMaster. *Bioinformatics* 2025, 41(Suppl 1):i342 — 보고값: supervised PROTAC-1K AUROC 0.854, semisupervised PROTAC-8K(random split) 0.882(모두 in-distribution). https://github.com/ABILiLab/DegradeMaster
29. **PROTAC-TS** — ChemTSv2 기반 PROTAC 링커 설계 + 세포막 투과(Caco-2, TabPFN) 예측. *JACS Au* 2025, DOI 10.1021/jacsau.6c00033. https://github.com/ycu-iil/PROTAC-TS

> ※ 문헌 확인 결과 "Mitragotri 식(logKp = logP − 0.0082·MW − 4.84)"은 출처가 불명확한 임의 계수로 판단되어 사용하지 않았다.
