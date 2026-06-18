# AR PROTAC의 안드로겐성 탈모(AGA) 국소 치료제 가능성 평가 및 AI 기반 신규 분자 설계

> 화학정보학·분자설계 기말 프로젝트 · 글로벌혁신신약학과 · 2026-06
> 모든 수치는 동결 스크립트(`scripts/phase0_freeze.py`) 출력을 단일 출처로 인용한다.

---

## ① 타겟 소개와 선정 이유

**안드로겐 수용체(Androgen Receptor, AR; UniProt P10275)** 를 타겟으로 선정하였다. AR은 핵수용체 전사인자로, 디하이드로테스토스테론(DHT)과 결합해 모낭(dermal papilla) 세포에서 모발 성장기를 단축시키는 안드로겐성 탈모(AGA)의 핵심 분자다. 기존 경구 AR 길항제·5α-환원효소 억제제(피나스테리드 등)는 전신 노출에 따른 성기능 부작용으로 장기 사용이 제한된다.

**PROTAC(proteolysis-targeting chimera)** 은 표적 단백질을 분해(degradation)하는 이중기능 분자로, 단순 점유(occupancy)가 아닌 촉매적 분해를 통해 낮은 농도에서도 작용한다. **국소(topical) AR-PROTAC** 은 "두피 모낭에서는 AR을 분해하되 전신 노출은 최소화"라는 분리(separation)를 달성할 수 있다면 이상적인 AGA 치료 전략이 된다. 본 프로젝트의 직접적 동기는 2024년 *J. Med. Chem.* 논문(**jm4c02226**, "Discovery of a Novel Non-invasive AR PROTAC Degrader for the Topical Treatment of AGA")이며, 이 논문의 29개 화합물(A/B/C 시리즈)을 검증 케이스로 삼는다.

타겟 선정 정당화의 보조 근거로 DegradoMap(타겟 수준 degradability 예측) 결과를 정성 인용하되, 이 도구는 PROTAC 구조를 입력받지 않아 **화합물 단위 랭킹에는 사용할 수 없으며**(target-unseen AUROC ~0.60, gradient boosting 대비 p=0.556) 신규 추론은 수행하지 않았다.

---

## ② 데이터/구조 확보 및 전처리

### 데이터 자산 ([표1])

**[표1] 전체 규모 vs AR(P10275)** — AR 별도 열

| 데이터셋 | 전체 행 | 전체 unique ik14 | AR 행 | AR unique ik14 |
|---|---|---|---|---|
| 구 PROTAC-DB (`ProtacDB_MTL_CLS.csv`) | 6,107 | 6,077 | 367 | 360 |
| 신 PROTAC-DB 3.0 (`protac.xlsx`) | 15,502 | 10,383 | 777 | 658 |
| case study (`jm4c02226`) | 29 | 29 | 29(전부 AR) | 29 |

(보조: 신DB AR new-only 골격 298[구조분할 test]; case study A16/B6/C7, **누수 B3·B5**.)

**[표1b] 실험 측정값(계산 descriptor 제외) 데이터 개수 — 전체 / AR.** MW·XLogP3·TPSA·HBA/HBD·RotBond·InChIKey 등 계산값은 제외.

| 측정값(실험) | 데이터셋 | 전체 | AR |
|---|---|---|---|
| **DC50 (nM)** | 신DB | 3,592 | 167 |
| **Dmax (%)** | 신DB | 2,729 | 139 |
| Percent degradation | 신DB | 2,774 | 252 |
| IC50(→Target) / IC50(→E3) | 신DB | 1,831 / 453 | 9 / 0 |
| Kd(→Target) / Kd(→E3) | 신DB | 501 / 287 | 0 / 0 |
| 세포 IC50 / EC50 | 신DB | 4,405 / 388 | 191 / 22 |
| Caco-2 A2B / B2A / PAMPA Papp | 신DB | 91 / 69 / 32 | 19 / 19 / 10 |
| pDC50 / Dmax / Degrader_Class (mask 관측) | 구DB | 1,082 / 574 / 507 | 62 / 52 / 51 |
| Skin Retention / DC50 / Dmax | case study | 29 / 13 / 13 | (전부 AR) |

> **시사점**: 전체 규모는 크나 **핵심 분해 라벨(DC50/Dmax)은 희소** — 신DB DC50은 전체의 23%·AR의 21%뿐. "데이터 많음"의 상당 부분은 계산 descriptor·결합/세포 활성이고, 분해 라벨 자체는 학습·검증에 빠듯하다(⑤-2·⑤-6).

- **AR 필터는 `Uniprot=='P10275'` 단일 기준**으로 강제하였다. `Target` 컬럼 정규식('AR'/'Androgen') 매칭은 AR-V7·T878A 등 변이체를 substring으로 흡수해 1193행까지 부풀려지므로 **사용 금지**하였다.
- **InChIKey 매칭**은 precomputed 컬럼이 RDKit과 불일치하므로, 양 DB 모두 RDKit `MolToInchiKey` **앞 14자(골격키)** 로 재계산해 통일하였다.
- **case study 파싱**: `cp949` 인코딩, `header=[0,1]` 다중헤더. '±'가 cp949에서 '÷'로 mojibake되어 있어 `÷→±` 치환 후 선행 float만 취하였다. 3번째 칸 헤더는 정확히 `Skin Retention Rate(%)`로 확인(A·B·C 29/29 전부 보유). DC50/Dmax는 **B/C 13개만** 보유(A16개는 활성 데이터 없음). 수치형(비검열) DC50은 5개(B3 8.97, B6 207.03, C1 103.87, C5 70.85, C6 199.5)이며 나머지는 `>1000` censored.
- **구조적 시간분할**: 구DB 골격 360개는 신DB 658개에 완전 포함(old⊂new), 신DB에만 존재하는 **298개 골격**이 prospective test 후보다. 구DB에는 DOI/연도 컬럼이 없어 엄밀한 연도분할은 주장하지 않는다. DOI를 Crossref로 해소(99%)하면 hold-out 게재연도 중앙값 2024 vs 기존 2022로, 멤버십 분할이 시간분할과 대체로 일치한다.

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

hold-out 1,134개는 모델 입력 요건상 **DC50·Dmax를 둘 다 보유**(그래서 둘 다 1,134). "나머지"는 DC50(1,590)>Dmax(921)로, **Dmax 없이 DC50만 있어 hold-out 요건(둘 다)을 못 채운 화합물**이 많다. 분해 라벨(DC50/Dmax)은 결합·세포활성보다 희소하다.

### 구조(3D) 확보의 한계
AR-LBD의 PROTAC 삼원복합체(AR+E3+PROTAC) 실험구조는 부재하다. 이는 ⑤·⑥에서 핵심 제약으로 다룬다.

---

## ③ 적용한 방법과 선택 이유

| 단계 | 방법 | 선택 이유 |
|---|---|---|
| 화학공간 | Morgan FP(r=2, 2048bit) + PCA, MW–TPSA overlay | bRo5 영역의 시각적 위치 확인 |
| 분해 예측 | **PROTAC-STAN**(GNN+ESM, 이진 분류) — 전체 타겟 prospective(학습셋 미관측) 검증(주) + AGA 29개 적용 | SMILES+서열만으로 가동; 학습-미관측 화합물에 prospective 평가로 도구 자체를 검증 |
| 피부 투과 | Potts-Guy QSPR(직접 구현) | 투명·재현 가능한 전신투과 proxy |
| 신규 설계 | SAR 기반 수동 medchem 편집 | 모델 점수 의존(순환논증) 회피 |

**분해예측 모델 선택의 핵심 제약.** DeepPROTACs·DegradeMaster·AiPROTAC 등 구조기반 예측기는 모두 **단백질 포켓(.mol2/PDB) + 3원 도킹 포즈**를 필수 입력으로 요구한다(코드 직접 확인). AR 삼원구조가 없어 2일 내 적용이 불가능하므로, 구조 입력이 필요 없는 **PROTAC-STAN을 주 분해예측 모델**로 채택하였다. "구조기반 예측기 3종(DeepPROTACs·DegradeMaster·AiPROTAC)은 AR 3원구조 부재로 적용 불가, 포켓-프리 도구(STAN·Ribes)만 가동"이라는 점 자체가 ⑤의 정직한 발견이다.

**PROTAC-STAN 적용 세부.**
- 사전학습 가중치(`protac-stan.pt`)를 **동결**하고 prospective 추론만 수행(AR 관측 양성 절대수 25개로 재학습 비현실적).
- 입력 `custom.csv`(29행)의 `E3 ligase Uniprot`는 **전 29행 CRBN(Q96SW2)** 로 지정하였다. 29/29 화합물이 글루타리미드(piperidine-2,6-dione, CRBN cue)를 보유하고 하이드록시프롤린(VHL cue)은 0/29이기 때문이다. B/C의 클로로-시아노-페녹시-사이클로헥실 모티프는 VHL 리간드가 아니라 **AR 길항제 워헤드(원논문 표기상 enzalutamide analog 계열)**다.
- 9개 물성 컬럼은 RDKit 재계산(XLogP3←Crippen MolLogP, TPSA←RDKit TPSA 등)으로 STAN의 `columns` 리스트와 1자 단위 일치시켰다.
- **재현성 패치**: torch 2.8에서 `torch.load`의 `weights_only=True` 기본값 때문에 데이터 파이프라인이 `UnpicklingError`로 실패한다. 원본 리포를 수정하지 않고 monkeypatch 래퍼(`scripts/stan_patch.py`)로 `weights_only=False`를 강제하여 해결하였다.
- 원본 `inference.py`는 argmax 라벨만 반환하므로, `F.log_softmax` 출력의 class-1 성분을 `exp`하여 **활성확률**을 덤프하도록 별도 러너(`scripts/run_stan_inference.py`)를 작성하였다(랭킹 가능).

**평가 임계값**: 모델군마다 양성 라벨 정의가 달라(STAN/DeepPROTACs=양성 AND vs DegradeMaster/DegradoMap=양성 OR) 통일 컷 단일 채점을 금지하고, STAN은 자기 기준(`DC50<100nM AND Dmax≥80%`)으로 평가하였다.

---

## ④ 결과

### [그림1] 화학공간 — AR PROTAC은 명백한 bRo5
구·신 DB와 case study를 Morgan FP–PCA로 투영(PC1 20.8%, PC2 9.5%, 합 30.3%). MW–TPSA 평면에서 B/C 시리즈는 **MW 725~792, TPSA 144~155** 로 Ro5 선(MW=500, TPSA=140)을 명백히 위반하고 bRo5 상한(MW≈1000, TPSA≈250, Doak 2014) 영역에 위치한다. A 시리즈(MW 466~597)는 상대적으로 작다.

### [표5]·[그림TS] 전체 타겟 prospective(학습셋 미관측) 도구 검증 (★ 주 결과)
"P10275만으로 도구 성능을 판단하지 말라"는 요구에 따라, 동결 PROTAC-STAN을 최신 PROTAC-DB 3.0의 **학습-미관측(prospective) 화합물 1,134개(51개 타겟, 양성 388/음성 746)** 에 적용해 일반 분해예측 성능을 검증하고, 학습-관측 281개를 memorization 대조군으로 두었다(STAN 자기 라벨 `DC50<100 AND Dmax≥80`; 분할=STAN 학습 골격 ik14 1,164개 대비 미관측 여부[엄밀 시간분할은 데이터 제약으로 불가→prospective 멤버십 분할로 대체], ESM 283 Uniprot 커버 타겟·E3 한정, compound-level dedup; 스크립트 `timesplit_build.py`/`run_stan_ts.py`/`timesplit_eval.py`, 로그 `outputs/timesplit_log.txt`).

| 코호트 | n | AUROC [95% CI] | AUPR | F1 | MCC |
|---|---|---|---|---|---|
| **prospective (미관측)** | 1,134 | **0.574** [0.539, 0.607] | 0.379 | 0.473 | **0.098** |
| leaked (관측, 대조) | 281 | **0.909** [0.873, 0.941] | 0.895 | 0.858 | 0.687 |
| prospective–AR | 42 | 0.702 [0.524, 0.854] | 0.700 | 0.731 | 0.347 |

- **Memorization gap**: 학습-관측 AUROC 0.909 → 미관측 0.574(MCC 0.098, 거의 무작위)로 붕괴([그림TS-1] `outputs/ts_roc.png`).
- **pooled 0.574조차 부풀려짐**: 타겟별(양·음 각≥5인 14개, 화합물 65% 커버) AUROC 중앙값 **0.448**, >0.5는 5/14, 최대 타겟 ERα(n=258)는 0.276 → 타겟 내 화합물 변별은 ≈무작위(n-가중 타겟내 AUROC ≈0.40). AUROC 비교 쌍의 ≈93%가 타겟간이라 pooled는 **타겟간 분리**를 반영할 뿐이다([그림TS-2] `outputs/ts_per_target.png`). Spearman(prob vs −DC50)=0.348(미관측) vs 0.614(관측).
- **결론(도구 검증):** STAN은 학습 분포 안에선 잘 맞으나 **새 PROTAC의 prospective 분해예측력은 약하다**(⑤-2). 화합물 랭킹 근거로 직접 쓰기 어렵다.

**▸ 다중 도구 교차검증 (STAN vs Ribes, +DegradeMaster 참조)** ([표6] `outputs/multitool_metrics.csv`·`multitool_pertarget.csv`, [그림TS-3] `outputs/fig_multitool.png`). 포켓-프리 예측기 **Ribes**(FP+XGB/MLP)를 동일 **mutually-held-out 786개**(STAN·Ribes 학습셋 교집합 0; Ribes 입력 커버리지가 정의한 부분집합이라 클래스·타겟 비균형)에 적용:

| 도구 | pooled (strict GT) | Ribes-native GT | KRAS 제외 | 타겟별 macro |
|---|---|---|---|---|
| PROTAC-STAN | 0.550 | 0.559 | **0.427** | 0.531 |
| Ribes std (random-CV) | **0.696** | **0.491** | 0.618 | 0.593 |
| Ribes tgt (group-CV) | 0.685 | 0.577 | 0.670 | 0.501 |

- **순위가 분석 렌즈마다 뒤집힌다**(라벨 정의·KRAS 단일클래스 블록[786 중 214개 거의 전부 음성]·분할): 공정한 타겟별 macro는 셋 다 **0.50–0.59(무작위 부근)**, Ribes-std는 자기 학습 라벨에서 0.491, STAN은 KRAS 제외 시 0.427로 추락 → 두 도구 모두 신규 화합물의 타겟 내 변별을 거의 못 함. "어느 모델이 낫다"는 라벨·구성 의존이라 신뢰 불가.
- **라벨·구성 무관하게 견고**: 두 도구 크게 불일치(proba ρ=0.19, **κ=0.13**; 0.5컷 호출일치 56%는 보조).
- 구조기반 **DegradeMaster**는 우리 세트 적용 불가(삼원 pocket 필요) → 저자 보고값 인용: supervised PROTAC-1K AUROC **0.854**(+6.9%), semisupervised PROTAC-8K(random split) **0.882**(+11.76%) — chemotype 공유 in-distribution 수치(STAN 학습-관측 0.909와 동류).
- **메시지:** 리더보드(0.85–0.91)는 대부분 in-distribution이고, held-out 신규 화합물에선 macro AUROC 0.50–0.59로 추락하며 **순위마저 불안정**하고 도구 간 합의도 약하다 → 단일 벤치마크 숫자 단독 신뢰 불가, 앙상블·실험검증 필수.
- **3-way 동일 held-out 셋이 불가능한 이유:** DegradeMaster 구조가 있는 유일 셋 PROTAC-8K(1,164골격)는 **STAN 학습셋과 100% 동일**·Ribes와 68% 중첩(STAN held-out 0개)이고, STAN·Ribes에 held-out인 786개는 3D구조가 없어 DegradeMaster 미실행 → "구조 구비 ∩ 세 도구 미관측" 양립 불가(벤치마크 contamination). 786개에 삼원구조 생성(Boltz/도킹) 없이는 공정 3-way 단일셋 없음.

### [그림2]·[표2] case study(AGA) — AR PROTAC 29개 활성확률 랭킹
- leakage-free 활성 앵커가 **모두 상위 랭크**: C5(70.85nM)=1.00, C1(103.87nM)=0.99, C6(199.5nM, 리드)=0.80.
- **누수 화합물 B3(8.97nM, 최강)·B5(실측 비활성: DC50>1000, Dmax 37%)의 확률**(0.25, 0.16). B5의 낮은 점수는 비활성에 대한 정답 예측이고, 활성인데도 낮은 사례는 B3 1건뿐이다 → 학습셋 누수가 점수를 단순 부풀리지 않았다(④의 핵심 관찰, ⑤에서 해석).
- jm4c B/C 13개에 STAN 자기 기준(`DC50<100 AND Dmax≥80`)을 적용하면 양성이 **0개**(최대 Dmax 76%)라 AUROC/AUPR/F1 산출이 불가능하다. 따라서 **주지표는 활성확률 vs 실측의 Spearman 순위상관**으로 한다.

| 지표 | rho | 95% CI(bootstrap 2000) | n |
|---|---|---|---|
| 확률 vs −DC50 (수치형, 누수 포함) | −0.10 | [−1.00, 1.00] | 5 |
| 확률 vs −DC50 (**leakage-free**) | **+0.80** | [−1.00, 1.00] | 4 |
| 확률 vs Dmax (B/C 전체) | +0.37 | [−0.24, 0.81] | 13 |

표본이 극히 작아(n=4~13) 신뢰구간이 매우 넓다. 결과는 **정량 검증이 아니라 정성적 경향/가설**로 해석한다.

### [그림3]·[표3] 피부 투과(전신누출 proxy) vs 피부 잔류
- Potts-Guy logKp = 0.71·logP − 0.0061·MW − 6.3 (RDKit 직접 구현). 29개 logKp는 **−8.22 ~ −6.74** 로 압축, B/C std=0.41에 불과 → MW항 지배로 **변별력 거의 상실**.
- **적용도메인(AD) 외삽**: MW>750이 11/29, TPSA>140이 13/29, AD 밖(MW>750 또는 TPSA>140)이 **13/29(44.8%)**. B/C 13개는 전부 TPSA>140 기준으로 AD 밖(MW>750은 11개, C4·C5는 MW≤750).
- logKp와 실측 retention의 Spearman = **0.042**. 단 logKp가 거의 상수로 압축(B/C std=0.41, 변별력 상실)되어 이 값은 정보를 거의 담지 못하므로, **"피부 잔류(depot) ≠ 피부 통과(Kp, 전신흡수)"** 가설과 일관될 뿐 그것을 입증하지는 않는다(⑤-1·⑤-7).
- 흥미롭게도 활성 B/C(적색 두꺼운 테두리)는 retention이 낮고, A 시리즈가 "이상 영역(낮은 Kp·높은 retention)"에 더 가깝다 — 활성과 잔류가 동시에 최적화되지 않음을 시사.
- **▸ PROTAC-TS 세포막 투과(Caco-2) 예측으로 "3중 도메인 구분" 실증** ([그림3c] `outputs/fig_protacts.png`, `protacts_predict.py`): PROTAC-DB Caco-2 89개로 학습한 TabPFN(LOOCV R²=0.78)으로 29개 예측 → 예측 세포막투과 vs 피부잔류 Spearman **+0.41**(약한 양), vs Potts-Guy logKp **−0.05**(무상관), vs STAN활성 −0.10. **세포막투과≠피부투과≠피부잔류가 데이터로 확인**. 누수 0. (29개는 PROTAC-TS 자체 AD 필터상 도메인 안쪽 — AD 외삽이 아니라 예측값 압축 한계; −0.05는 logKp도 압축돼 독립 단정보다 "서로 예측 못 함"으로 해석.)

### [표4] 신규 분자 3종 (가설 생성)
leakage-free 활성 앵커의 약리단(CRBN 글루타리미드 + AR 길항제 워헤드(enzalutamide analog 계열))을 **보존**하고 말단만 보수적으로 변형:

| ID | 부모 | MW | logP | TPSA | AD 최근접 Tanimoto | STAN weak prior | 설계 근거(SAR) |
|---|---|---|---|---|---|---|---|
| D1 | C1 | 770 | 4.80 | 147 | 0.70 | 0.998 | CRBN benzamide ortho-F (대사 차단·logD 조정) |
| D2 | C5 | 751 | 5.49 | 144 | 0.65 | 0.995 | phenyl-glutarimide ortho-Me (링커 회전 제약) |
| D3 | C6 | 783 | 6.10 | 144 | 0.63 | 0.835 | AR 워헤드 인접 F (시아노-할로-아릴 워헤드 다중치환 모사) |

STAN 점수는 weak prior일 뿐이며 실험 검증이 필수다(⑤·⑥).

---

## ⑤ 결과에 대한 비판적 고찰 (가장 중요)

도구가 내놓은 숫자를 어디까지 믿을 수 있는가가 본 프로젝트의 핵심이다.

**5-1. endpoint 선택 오류의 발견과 수정.** 초기 가설은 "높은 logKp = 좋은 국소약"이었으나 이는 개념 오류다. 국소 AGA 약의 목표는 **높은 피부 잔류 + 낮은 전신 투과**이며, logKp(통과)와 retention(잔류)은 음의 관계이거나 무상관일 수 있다. 실제 Spearman은 0.042였다. 다만 logKp가 거의 상수로 압축(B/C std=0.41, 변별력 상실)되어 이 상관계수 자체가 정보를 거의 담지 못하므로, 0.042는 depot≠flux 가설과 일관될 뿐 입증하지는 않는다(실증 근거는 C6 skin-to-plasma >1,700배 등 독립 증거에 둔다). endpoint를 통과(Kp)가 아닌 잔류로 바로잡은 방향은 옳으며, logKp를 "전신 누출 proxy(낮을수록 유리)"로 재정의하였다.

**5-2. [핵심] 전체 타겟 prospective(학습셋 미관측) 검증 — STAN의 memorization gap.** 본 프로젝트 도구 검증의 핵심 결과다. 동일 모델·지표에서 학습-관측 AUROC 0.909 vs 미관측 0.574(MCC 0.098) — 실제 사용(미관측 신규 화합물)에서 신뢰할 값은 0.574 쪽으로 거의 무작위다. 더구나 pooled 0.574는 타겟별 중앙값 0.448(n-가중 타겟내 ≈0.40)을 가린다(ERα n=258 0.276): AUROC 비교 쌍의 ≈93%가 타겟간이라 pooled는 타겟간 분리를 반영할 뿐 타겟 내 화합물 변별은 못 한다(단 모델이 타겟 양성률을 *명시적으로* 학습한 것은 아니다 — 타겟 평균확률 vs 양성률 상관 p≈0.10, "타겟 양성률" 오라클 0.839 ≫ 0.574). 따라서 STAN 점수를 화합물 랭킹 근거로 직접 쓰면 위험하며, 신규 설계(5-11)에서 weak prior로만 쓴 결정이 정당화된다. 이 한계는 단일 타겟·소표본으로는 드러나지 않았을 것 — "전체 타겟+leakage 대조+타겟별 분해" 설계가 도구의 실체를 드러냈다.

**5-3. PROTAC-STAN은 분류기다.** 모델은 2-class 이진 분류기(`config.toml class=2`)로 **연속 DC50(nM)을 재현할 수 없다.** 따라서 RMSE/R²/Pearson은 원천적으로 산출 불가하며, 활성확률 랭킹과 Spearman 순위상관만 의미를 가진다.

**5-4. 코호트 의존 z-score.** `data.py`가 입력 코호트 내부에서 per-column z-score를 계산하므로 동일 SMILES도 코호트가 바뀌면 점수가 변한다. 이를 절차로 못박아 **항상 동일한 29행으로 1회만 추론**하고 누수 표시는 사후 마스킹으로만 하였다. 신규 설계 3종은 3행만으로는 일부 물성의 std=0 → NaN이 발생하여, 29개 참조 코호트에 합친 **32행 코호트**에서 추론하였음을 명시한다.

**5-5. 데이터 누수와 그 역설.** case study 29개 중 B3·B5가 STAN 학습셋 및 구DB와 골격키 일치(누수)다. 하필 B3는 최강(8.97nM) 화합물이라 이전 분석의 핵심 앵커였으나 **memorization 위험**으로 주 근거에서 제외하였다. 누수된 B3·B5의 STAN 확률은 각각 0.25, 0.16으로 낮았다 — 다만 B5는 실측 비활성(DC50>1000, Dmax 37%)이라 낮은 확률은 정답 예측이고, 활성(8.97nM)인데도 낮은 사례는 B3 1건(n=1)뿐이다. 따라서 누수가 점수를 부풀리지 않았다는 정도만 말할 수 있을 뿐 활성 학습 실패를 일반화할 수는 없다. leakage-free 앵커(C5/C1/C6)만으로 정성 일치를 보고한다.

**5-6. 단일 E3(CRBN) 시리즈.** 29/29가 CRBN 기반이므로 STAN의 E3 채널은 상수다. 즉 랭킹은 사실상 PROTAC 그래프 + 물성 z-score만으로 결정된다. ESM E3 임베딩의 변별 기여는 본 데이터에서 검증 불가하다.

**5-7. Potts-Guy의 AD 외삽.** Potts-Guy는 Flynn 데이터셋(대부분 MW<500 수준의 소분자) 회귀식으로, 전체 29개 중 44.8%(=B/C 13개 전부)가 AD 밖이다(B/C는 TPSA>140 기준 13/13). MW항 지배로 모든 PROTAC이 균일하게 낮은 logKp로 압축되어 변별력을 잃는다(예상된 결과). 절대 Kp 값은 정량 신뢰 불가, 상대 랭킹·정성 용도로만 사용하였다. 또한 모낭(follicular) 경로·유한용량·제형 효과를 무시하는데, 이는 두피 국소제의 핵심 메커니즘이다.

**5-8. 3중 도메인 구분(실측 확인).** 세포막 수동투과(Caco-2/PAMPA) ≠ 피부 투과 ≠ 피부 잔류 — **PROTAC-TS Caco-2 예측으로 실증**(④[그림3c]): 예측 세포막투과는 피부잔류와 +0.41, Potts-Guy logKp와 −0.05(무상관). 카멜레온성은 TPSA(2D, conformer 무관 고정값)로만 근사했고, "IMHB로 TPSA 감소"는 범주 오류다(변하는 것은 3D EPSA). 거대 PROTAC의 3D conformer 샘플링은 비용·실패율 문제로 NICE-TO-HAVE로 강등하였다.

**5-9. 표본 크기(AR/case study 한정).** retention n=29, DC50+Dmax n=13, 수치형 DC50 n=5에 불과해 모든 상관의 bootstrap CI가 [−1,1]에 걸친다. 정량 결론이 아닌 가설 생성으로 톤다운한다.

**5-10. 재현성.** torch 2.8 비호환을 5줄 monkeypatch로 해결한 사실, argmax→확률 덤프 패치 사실을 명시한다.

**5-11. 신규 설계의 순환논증 위험.** Phase 2/3에서 한계가 드러난 STAN으로 자기 설계물을 "검증"하면 순환논증이다. 따라서 설계 근거를 모델 점수가 아니라 SAR 규칙(약리단 보존 + 보수적 말단 변형)에 두고, STAN 점수는 weak prior로만, AD 거리(Tanimoto 0.63~0.70)를 동반해 제시하며, 결론은 "실험 검증 필요"다.

---

## ⑥ 한계와 다음 단계

**한계.**
- **PROTAC-STAN prospective 분해예측력 약함**(주 발견): 미관측 AUROC ≈ 0.57·타겟 내 변별 ≈ 무작위(타겟별 중앙 0.448, n-가중 0.40) vs 학습-관측 0.909 — memorization gap. 화합물 랭킹 신뢰 제한(④[표5]·5-2). 분할은 ik14 멤버십(학습셋 미포함) 기준이나, Crossref로 DOI→게재연도 99% 해소 후 검증 시 **hold-out 중앙값 2024 vs 기존 2022**로 사실상 시간 분할과 대체로 일치(분포 겹침=하드 컷오프 아님; [그림 pubdate] `outputs/pubdate_kde.png`).
- 도구 검증 범위는 ESM 사전이 커버하는 **283 Uniprot → 평가셋 51개 타겟·E3 5종(CRBN/VHL 중심)**으로 제한.
- AR 삼원복합체 구조 부재로 구조기반 분해예측기 3종(DeepPROTACs/DegradeMaster/AiPROTAC) 적용 불가; 포켓-프리(STAN·Ribes)만 적용 가능.
- AR 관측 양성 절대수 25개로 AR-only 재학습 비현실적.
- jm4c B/C에 STAN 자기 기준 적용 시 양성 0개 → 분류지표 산출 불가, 순위상관만.
- 피부 잔류의 실험 조건(도포 후 시점·ex vivo/in vivo 모델)에 대한 정규화가 일관되지 않을 수 있다.

**다음 단계.**
1. Boltz/AlphaFold3로 AR-LBD+CRBN(±DDB1)+PROTAC 삼원복합체를 예측(리드 C6 1케이스부터)하여 구조기반 예측기를 활성화.
2. 포켓-프리 리간드 기반 예측기(PROTAC-Degradation-Predictor, Ribes 2024)를 추가해 STAN과 "모델 간 불일치"를 정량.
3. ETKDG 다중 conformer로 3D EPSA·ΔPSA·IMHB를 계산해 카멜레온성을 정량(앵커 5~6개부터).
4. 신규 3종 D1~D3의 합성·HDPC AR 분해 assay·피부 잔류 측정으로 weak prior를 실험 검증.
5. 모낭 표적 전달(나노입자·제형)을 고려한 전신/국소 분리 모델링.

---

### 부록 — 재현 정보
- 환경: `~/anaconda3/envs/protac` (Python 3.10, torch 2.8.0+cu128, RDKit 2026.03.1).
- 스크립트(`scripts/`): `phase0_freeze.py`, `stan_patch.py`, `run_stan_inference.py`, `build_custom_csv.py`, `phase2_analysis.py`, `fig1_chemspace.py`, `phase3_skin.py`, `phase5_design.py`.
- 산출물(`outputs/`): `table1~4.csv`, `fig1_chemspace.png`, `fig2_ranking.png`, `fig3_separation.png`, `fig3b_chameleon.png`, 각 phase 로그.
- 참고문헌: Potts & Guy (1992); Flynn (1990); Doak et al. (2014); jm4c02226 (*J. Med. Chem.* 2024); PROTAC-STAN; DegradoMap. (※ SETUP의 "Mitragotri 식"은 출처불명 임의 계수로 확인되어 사용 철회.)
