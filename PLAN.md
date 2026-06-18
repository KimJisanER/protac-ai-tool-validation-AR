# 실행계획 (PLAN.md) — AR PROTAC의 AGA 치료제 가능성 평가 + AI 신규 분자 설계

> 대학원 화학정보학·분자설계 기말 프로젝트. 마감 2026-06-19(금) 자정. 오늘 2026-06-17. 남은 시간 약 2일.
> 보고서 10쪽 이내, 한국어. 본 계획은 데이터/코드/환경을 직접 실행해 검증한 사실과 2건의 적대적 검토(과학적 정확성·실현가능성)에서 제기된 critical·major 지적을 전부 반영해 이전 PLAN(=`PLAN_v1_previous.md`)을 대체한다.
> **기조: "결과가 안 나와도 좋다. 실패한 시도와 원인 분석도 훌륭한 보고서"라는 교수 지침(과제 docx 원문에서 ⑤ '결과에 대한 비판적 고찰'이 명시적으로 '(가장 중요)'로 표기됨)을 적극 활용한다. 도구가 내놓는 숫자를 어디까지 믿을 수 있는지 따지는 것(채점 ⑤)이 가장 중요하다.**

---

## 0. 검증된 사실 및 이전 계획 대비 수정사항

아래 표의 **[CRITICAL]** 항목은 그대로 두면 분석의 타당성 자체가 무너지는 것이므로 반드시 반영한다. 모든 수치는 §0-END의 "동결 스크립트"가 단일 출처이며, 보고서·표·그림에 하드코딩 금지(스크립트 출력만 인용).

| # | 이전 가정 | 코드/데이터로 검증한 결과 | 조치 |
|---|---|---|---|
| **C1** | jm4c02226 3번째 칸 = 피부 잔류율, "예측 Kp vs 실측 retention 양의 상관"이 핵심 분석 | **[검증완료]** 헤더 라벨은 정확히 `Skin Retention Rate（%）`이며 **A·B·C 전 29행 모두** 값 보유(파싱 확인). 값 분포 A 0.68~9.12 / B 0.58~1.41 / C 0.81~2.48 — 모두 한 자릿수 %, 동일 스케일. 단, **A(16개)는 retention만 있고 DC50/Dmax 전무**, B/C(13개)만 DC50/Dmax 보유. 3번째 칸은 단일 의미(Skin Retention)로 확정 | A vs B/C는 척도가 같으므로 retention축은 공유 가능. 다만 **A는 활성(DC50/Dmax) 데이터가 없어 "활성-잔류 2D 플롯"에 활성색으로 올릴 수 없음** → Phase 3의 활성 연동 분석은 B/C(13개) 중심, A는 retention 분포 참고용으로만. SI 캡션 1회 교차확인을 §1-3 체크리스트로 유지(가점) |
| **C2** | logKp가 높을수록 좋은 국소약 | **[CRITICAL]** 개념 오류. 국소 탈모약 목표 = **높은 피부 잔류 + 낮은 전신 투과**(C6 skin-to-plasma >1700배가 이 분리의 실증). Kp(통과)와 retention(잔류)은 음의 관계이거나 무상관일 수 있음 | logKp를 **"전신 누출(undesired flux) proxy, 낮을수록 국소제로 유리"**로 명명. "낮은 상관"을 실패로 적지 않고 "애초에 다른 물리량"으로 프레이밍(§1, Phase 3) |
| **C3** | 분해모델로 DC50 회귀(RMSE/R²/Pearson) 평가 | **[CRITICAL]** PROTAC-STAN은 **2-class 이진 분류기**(`config.toml class=2`, `model.py` h_net out=2, `inference.py` `torch.max(...,dim=1)` argmax). 연속 DC50 출력 불가 → RMSE/R²/Pearson 원천 산출 불가 | Phase 2 지표를 **이진분류(AUROC, AUPR, F1, MCC)** + 활성확률 랭킹으로 전면 교체. argmax만으로는 랭킹 불가 → **class-1 확률 덤프 패치 필수**. "왜 회귀가 아닌 분류인가(모델 설계 제약)"를 ③에 명시 |
| **C4** | Train=구DB AR / Test=신DB AR로 AR-only 재학습 후 time-split | **[CRITICAL]** 사전학습 가중치 제공(`protac-stan.pt`). **AR 관측 양성 절대수가 25개에 불과**(구DB P10275 367행 중 `Degrader_Class_mask==1`은 51행, 그중 양성 25 / 음성 26) → 멀티엔티티 GNN 재학습 비현실적 | **AR-only 재학습 포기.** 사전학습 가중치 **동결(frozen)** 후 unseen 셋에서 **prospective 추론 평가**로 재설계. 재학습 불가 사유 = "AR 관측 양성 절대수(25) 부족"(불균형이 아님, C14 참조)을 ⑥에 명시 |
| **C5** | 2024 신규 논문이니 학습셋에 없음(누수 無) | **[CRITICAL/실측]** InChIKey 골격(앞14자) 매칭: case study 29개 중 **정확히 2개가 PROTAC-STAN 학습셋(`train_compound_smiles.csv`, 968 unique ik14) 및 구DB와 일치 = 누수.** 누수 정체는 **B3(QJZUVOHJKPMREF), B5(IFCKDGNCNVHZTM)**. **하필 B3가 이전 plan의 핵심 concordance 앵커(8.97nM)** → memorization 위험 | 누수 2건 정체를 [표2]/§C5에 명시. **B3는 leaked → STAN concordance 주 근거에서 제외**(참고용 마킹). leakage-free 주 앵커는 **C6(199.5nM), C1(103.87nM), C5(70.85nM)**(모두 비누수, 실측 확인). 누수포함/제외 **이중 보고** |
| **C6** | DegradeMaster/DeepPROTACs는 체크포인트 있으니 예측만 | **[CRITICAL]** 두 모델 모두 **단백질 포켓 .mol2 + 3원 도킹 포즈**(또는 DGL 포켓그래프) 입력 요구. AR 포켓·29개 3원 도킹은 2일 내 비현실적. AiPROTAC는 학습 .pth 부재 + py3.7/dgl 환경 별도 | **PROTAC-STAN을 단독 PRIMARY 분해모델로 채택**(SMILES+물성만으로 가동). 포켓 의존 3모델은 "AR 3원구조 부재로 적용 불가"를 ⑤ 핵심 한계로 서술 |
| **C7** | DegradoMap을 분해예측 후보로 나열 | **[검증]** DegradoMap은 **AlphaFold 구조 + E3 정체성만** 입력, PROTAC 구조 불입력. 단일 타겟 AR 내 29개 화합물 랭킹 **불가**(target-unseen AUROC ~0.60, GB 대비 p=0.556) | 역할 재정의: **①타겟 선정 정당화 전용**(per-compound 점수 금지). **기존 `results/*.json` 정성 인용으로만 확정, 신규 추론은 시도 안 함**(NICE에서도 환경 비호환 리스크 큼) |
| **C8** | SETUP의 Mitragotri 식 `logKp=logP-0.0082MW-4.84` | **[검증]** 출처불명 임의 계수. 실제 Mitragotri(2002-03 SPT) 모델 아님. Potts-Guy와 함수형 동일(공선) → 독립검증 가치 없음 | "Mitragotri 모델" 표기 **철회**. Potts-Guy 1종을 정식 인용. 두 번째 비교군 필요 시 **출처 있는 QSPR**(예: Magnusson MW-only) 또는 "단순 휴리스틱"으로 명시 강등 |
| **C9** | Potts-Guy를 PROTAC 29개에 적용해 정량 Kp 보고 | **[검증]** Potts-Guy는 Flynn 데이터셋(주로 MW<750 소분자) 회귀. B/C MW 725~792·TPSA 144~155는 **적용도메인(AD) 밖 외삽**. −0.0061·MW 항이 모든 PROTAC을 균일하게 "투과 안 됨"으로 압축 → 변별력 상실 | Kp는 **정량 신뢰 불가, 상대 랭킹/정성 용도**로만. AD 밖 비율을 표로 정량화하고 "예상 실패=결과"로 프레이밍(Phase 3, ⑤) |
| **C10** | 카멜레온성 = TPSA vs RotBonds 산점도, "IMHB로 TPSA 감소" | **[검증]** TPSA는 원자기여 합으로 conformer와 무관하게 **고정**, IMHB로 안 변함 → **"IMHB로 TPSA 감소"는 범주 오류.** 변하는 것은 3D EPSA/3D PSA | 카멜레온성 proxy 업그레이드: ETKDG 다중 conformer → **min3D_PSA, ΔPSA(=TPSA−min3D_PSA), IMHB 개수(거리/각도), Fsp3, PBF/globularity, Rg**. 설계전략을 "conformer별 3D EPSA 저하(TPSA 동일해도)"로 정정. **단, 비용 문제로 NICE-TO-HAVE 강등**(C17) |
| **C11** | AR PROTAC 30개 검증셋 | **[검증]** cp949·header=[0,1] 파싱 시 **29행**(A16/B6/C7). SMILES 29/29 정상 파싱. **A=16, B=6, C=7.** retention 29/29 보유, DC50/Dmax는 B/C **13행만**. 수치형 DC50는 **5개**(B3 8.97[누수], B6 207.03, C1 103.87, C5 70.85, C6 199.5), 나머지 '>1000' censored. A계열 DC50 전무 | "30개"를 "**29개(A16/B6/C7), retention 29, DC50+Dmax 13(B+C), 수치형 DC50 5개, '>1000' censored 다수**"로 정정. 모든 상관/지표에 **bootstrap 95% CI + leave-one-out 민감도** 의무화. 결론을 "정성 경향/가설"로 격하 |
| **C12** | AR PROTAC = protac.xlsx 777행 (Target 기준) | **[검증]** **`Uniprot=='P10275'` = 777행(658 unique ik14)** 사용. Target 정규식('AR'/'Androgen')은 substring 매칭으로 1193까지 부풀려짐(AR-V7, T878A 등 변이체 혼입) → **사용 금지**. 구DB(P10275)=367행(360 unique). 신DB AR DC50 78.5% 결손, Dmax 82.1% 결손 | AR 필터는 **`Uniprot=='P10275'` 단일 기준**으로 명문화. 변이체(AR-V7 등)는 부속 컬럼으로 분리 보고 |
| **C13** | 시간분할은 DOI 연도 기준 | **[검증]** 구DB(`ProtacDB_MTL_CLS.csv`)는 **DOI/연도 컬럼 없음**. 신DB AR에만 있는 **298 unique(ik14)**가 구DB에 부재(=test 후보; new 658 − old 360 overlap = **298**, old⊂new). 깨끗한 **구조적 분할** 가능 | 1차 기준을 **구조적 분할(old=train후보 / new-only 298 AR=test)** + 모델별 학습 InChIKey 누수 제거로. DOI 연도는 보조 annotation으로만(엄밀 연도분할 주장 안 함, DOI 커버리지 낮음). **이전 plan의 '297' → 검증값 '298'로 통일**(MolToInchiKey 앞14자 골격키 기준) |
| **C14** | Degrader_Class 254/6107(4.2%) 심한 불균형 → class weight·AUPR·재학습 포기 근거 | **[CRITICAL 정정]** `*_mask` 적용 시 **관측 라벨(`Degrader_Class_mask==1`)은 507행, 그중 양성 254 / 음성 253 ≈ 50:50 균형**(실측). 4.2%는 라벨 미관측(0-imputed 5600행)을 분모에 넣은 수치로 mask를 쓰는 순간 **정의상 모순**. 실제 "심한 불균형"은 사실이 아님 | C14를 **"관측 라벨 507행 ≈ 50:50(균형)"**으로 정정. 4.2%는 **"라벨 결손율(미관측 5600/6107 = 91.7%)"**로 명칭 변경해 별도 보고. **class weight 필요성은 mask-적용 후 실제 비율(균형)로 재판단 → 불필요**. AUPR은 정보용으로 병기하되 "불균형 대응"이라는 근거는 삭제. 재학습 포기 사유는 C4로(AR 양성 절대수 25) 일원화 |
| **C15** | protac.xlsx 기술자/InChIKey를 RDKit 재계산 | **[검증]** xlsx에 MW, XLogP3, TPSA, HBA/HBD, RotBond, Ring Count precomputed 존재. 단 **precomputed 'InChI Key' 컬럼은 RDKit과 불일치 → cross-DB 매칭엔 RDKit 재계산(MolToInchiKey 앞14자)을 양측 동일 적용** | Phase 1 기술자는 xlsx precomputed 1차 사용(시간절약). **InChIKey 매칭은 양 DB 모두 RDKit 재계산.** logP 출처(XLogP3 vs RDKit MolLogP)는 Potts-Guy 입력에서 통일하고 명시. 신DB는 **89열**(이전 plan '92열' 정정) |
| **C16** | bRo5 = MW>700 & TPSA>140 | **[검증]** 임의 컷, 표준 정의 아님 | Ro5 선(MW=500, TPSA=140)과 bRo5 상한(MW≈1000, TPSA≈250, Doak 2014) 복수 참조선으로 overlay. 임의 컷은 보조선·근거 인용 명시 |
| **C17** *(신규)* | PROTAC-STAN custom/ 파이프라인 "실행 입증됨" | **[CRITICAL/실측]** torch 2.8 환경에서 `data.py:181`의 processed 캐시 로드(`protac.pt`=torch_geometric Data/slices)가 `weights_only=True` 기본값 때문에 **`UnpicklingError`로 즉시 실패**(데모 custom 데이터로 재현 확인). 모델 가중치 `protac-stan.pt`는 plain state_dict라 로드되지만, **데이터 파이프라인은 패치 없이는 0% 가동.** monkeypatch로 `torch.load(..., weights_only=False)` 적용 시 데모 10행 추론 정상 출력 확인 | **Day1 최우선으로 STAN 스모크 테스트.** `torch.load` 호출부(`data.py:181`, `data_loader.py:47/49/52`, `inference.py:69`)에 **`weights_only=False`(또는 `torch.serialization.add_safe_globals`) 패치 명문화.** 위험표의 'STAN 실패=낮음'을 "환경 비호환으로 확실 실패하나 5줄 패치로 해결"로 정정. 패치 사실을 ②/⑤(재현성·환경 의존)에 기록 |
| **C18** *(신규)* | 입력 파일 = `data/custom/raw.csv`(3열: Uniprot, E3, Smiles) | **[CRITICAL/실측]** `data.py:179`는 `{name}.csv`를 읽음 → `--name custom`이면 실제 읽는 파일은 **`custom.csv`**(`raw.csv` 아님). 또한 `process()`는 `columns` 9개 precomputed 물성(`Molecular Weight, Exact Mass, XLogP3, Heavy Atom Count, Ring Count, Hydrogen Bond Acceptor Count, Hydrogen Bond Donor Count, Rotatable Bond Count, Topological Polar Surface Area`)을 필수로 요구하고 `columns[1:]`에 z-score 적용 → 3열만 있으면 **KeyError로 실패** | 입력 파일을 **`data/custom/custom.csv`**(29행)로 작성, 컬럼 = `Uniprot, E3 ligase Uniprot, Smiles` + **위 9개 RDKit 물성**(컬럼명을 `data.py`의 `columns` 리스트와 1자 단위로 정확히 일치). XLogP3←RDKit `MolLogP`, TPSA←RDKit `TPSA`, 나머지 대응 RDKit 함수. **`raw.csv` 언급 전부 삭제.** 재실행 전 `data/custom/processed/custom` 삭제(캐시 재생성) |
| **C19** *(신규)* | A=CRBN(Q96SW2), B/C=VHL(P40337)로 E3 매핑 | **[CRITICAL/실측]** 29개 SMILES를 RDKit substructure로 직접 검사: **글루타리미드(piperidine-2,6-dione, `O=C1NC(=O)CCC1`) 29/29 보유**, 하이드록시프롤린 **0/29**, 티아졸(VHL cue) **0/29**. 즉 **전 계열 CRBN 기반.** B/C의 클로로-시아노-페녹시-사이클로헥실 모티프는 VHL 리간드가 아니라 **darolutamide/ARN-509류 AR 워헤드**. `custom.csv`의 'E3 ligase Uniprot'는 ESM 임베딩(`esm_s_map.pkl`) 키로 **직접** 쓰이므로, 잘못된 VHL 라벨은 B/C 13개에 엉뚱한 E3 채널을 주입해 모든 하류 확률·랭킹·[그림2]/[표2]를 오염 | **`custom.csv`의 'E3 ligase Uniprot'를 29행 전부 CRBN=Q96SW2로 지정.** 이상적으로는 코드 자동판정(글루타리미드 SMARTS 매치→CRBN, 하이드록시프롤린 매치→VHL)으로 작성. A vs B/C 구분은 **E3가 아니라 워헤드(B/C=darolutamide류 AR 길항제 vs A=다른 워헤드)·링커**로 재서술. **가점: 단일 E3(CRBN) 시리즈라 STAN의 E3 채널이 상수 → 랭킹이 사실상 PROTAC 그래프+물성 z-score만으로 결정됨**을 ⑤에 명시 |
| **C20** *(신규)* | 분해모델 평가를 단일 `DC50<100nM` 컷으로 통일 채점 | **[CRITICAL/실측]** 모델마다 양성 라벨 정의가 다름(코드·논문 직접 확인): **PROTAC-STAN = `DC50<100nM AND Dmax≥80%`**(코드 `clean.ipynb` cell16-17, Dmax는 **inclusive `≥`** + 단일값 완화·%-백필), DeepPROTACs/AiPROTAC = `DC50<100 AND Dmax>80%`(strict AND, 관례 기원=DeepPROTACs 2022), **DegradeMaster/DegradoMap = `DC50<100 OR Dmax≥80%`(양성=OR, 음성만 AND)**, 우리 GT(`Degrader_Class` mask==1 507행) = `pDC50>7 AND Dmax>0.80`. **AND 위치가 모델군마다 반대** → 같은 숫자라도 양성집합 상이. 또한 **jm4c B/C 13개에 STAN AND 기준 적용 시 양성 0개**(최대 Dmax=76%), `DC50<100` 단독이어도 양성 2개(B3[누수]/C5)뿐 → 이진분류 지표(AUROC/AUPR/F1) **산출 불가** | 각 모델은 **자기 라벨 정의로 평가**(통일 컷 금지). 주지표를 **class-1 확률 vs 실측 DC50/Dmax Spearman 순위상관**으로, 분류지표는 양성≥2 컷에서만 보조. 사용자 가설 `Dmax>80`은 STAN 코드 기준 `Dmax≥80`로 부등호 정정 |

> **§0-END 동결 스크립트:** Phase 0에서 위 모든 카운트(777/658, 367/360, 298, 507/254/253, 91.7%, AR 양성 25, 누수 B3/B5, 신DB 89열, 29=A16/B6/C7, DC50 13행/수치형 5개, jm4c AND-양성 0개·DC50<100단독 양성 2개[B3·C5])를 1개 스크립트(`scripts/phase0_freeze.py`)가 출력하게 하고, [표1]·본문은 이 출력만 인용. **하드코딩·기억 인용 금지.**

> **사용자 확정 결정 (2026-06-17):**
> 1. **STAN `torch.load` 패치 = monkeypatch 래퍼**(`scripts/stan_patch.py`) — 원본 리포 미수정, 재현성·git 청결 유지(§Phase 2 절차 1).
> 2. **두 번째 분해모델 = DegradeMaster NICE 슬롯 추진** — 학습데이터는 Zenodo(rec. 14728925)/GitHub에서 확보 가능하나, **코드 확인 결과 우리 29개 적용은 3원(ternary) 구조 입력(AR/CRBN 포켓 PDB + 포지셔닝 mol2 + PROTAC 3D)이 필수 → AR 3원구조 부재로 미적용**. [최종: 미적용 문서화만] 설치·데모 생략, '구조기반 모델 3종 미적용 vs STAN 단독 가동'을 ⑤ 공통 한계로 서술(§Phase 2 SECONDARY).
> 3. **Phase 5 신규설계 E3 = CRBN 확정**(VHL 탐색 안 함; case study 29/29가 CRBN/글루타리미드라 일관).
> 4. **case study 논문(jm4c02226) PDF 보유** → SI에서 Skin Retention **측정 조건(도포 후 시점·ex vivo/in vivo 피부 모델·정규화 기준)**을 확인해 ②/⑤에 기록(C1 가점). PDF를 세션에 공유하면 직접 파싱.

---

## 1. 핵심 개념 정정 (보고서 ③·⑤에 그대로 사용)

### 1-1. Skin Retention(피부 잔류) ≠ 투과계수 Kp(피부 통과)
- **logKp (Potts-Guy, cm/h)**: 정상상태·무한용량·수용액 조건에서 각질층(SC)을 **통과해 수용체상/혈류로 빠져나가는** 투과계수 = **전신 흡수 플럭스** 지표.
- **Skin Retention Rate(%)**: 도포량 중 피부 조직 **안에 잔류하는** 분율 = **depot 양**.
- **국소 탈모약(AGA) 설계 목표 = 높은 피부 잔류 + 낮은 전신 투과.** C6의 skin-to-plasma >1700배가 바로 이 "분리"를 의미.
- 따라서 잔류가 높은 화합물이 오히려 Kp가 낮을 수 있어 **둘은 음의 관계이거나 무상관일 수 있음.** 이전 plan의 "높은 Kp = 좋은 약" 암묵 가정은 이 적응증에서 정반대.
- **재정의된 분석(Phase 3):** logKp를 "전신 누출 proxy(낮을수록 유리)"로 명명하고, 실측 retention과 직접 상관시키지 않는다. 대신 **2D 산점도(x=예측 logKp 전신투과축, y=실측 retention 피부저류축)**를 그려 **"이상적 사분면(낮은 Kp·높은 retention)"**을 시각화한다.
- ⑤ 가점 포인트: "**평가지표(endpoint) 선택 오류를 발견·수정한 과정 자체**"를 비판적 고찰로 서술.

### 1-2. 세포막 투과 ≠ 피부 투과 ≠ 피부 잔류 (3중 구분)
- PROTAC-TS(Caco-2/PAMPA류)는 **세포막 수동투과**를 예측 — 인지질 단일막/세포 단층, 짧은 거리. 피부(다층 SC, mm 두께, depot)와 물리적으로 다른 장벽.
- "PROTAC-TS 예측 vs 피부 잔류 직접 상관"은 **두 단계 도메인 점프**. 대신 "세포막 투과·피부 잔류가 **같은 분자기술자(logP, TPSA, HBD)**에 의해 같은 방향으로 움직이는가"를 기술자 매개로 **간접** 비교.

### 1-3. case study 3번째 칸 — SI 교차확인 체크리스트(확인 완료, 보고서 ② 기록용)
원논문: **J. Med. Chem. 2024, jm4c02226** "Discovery of a Novel Non-invasive AR PROTAC Degrader for the Topical Treatment of AGA".
- [x] 3번째 칸 헤더 = `Skin Retention Rate（%）`, **A·B·C 29/29 전부 값 보유**(파싱 확인). 단일 의미로 확정.
- [x] 4번째 칸 = `AR Degradation in HDPCs / DC50 (nM)`, 5번째 = `Dmax (%)` → **A1~A16은 DC50/Dmax 없음(측정 안 됨), B/C 13개만 보유.**
- [x] '±'가 cp949에서 '÷'로 mojibake(예 `8.45 ÷ 0.16` = 8.45 ± 0.16). 파싱 시 선행 float만 취함.
- [x] C6만 PK 컬럼(T1/2, Tmax, Cmax, AUC) 보유 = 리드 화합물. C6 skin-to-plasma 분리는 본문 인용.
- **결론(분기 (a) 확정):** col-3는 전 계열 Skin Retention. Phase 3 2D 분리 시각화는 **활성색이 필요하므로 B/C(13개) 중심**, A(16개)는 활성 데이터 부재로 retention 분포 참고용. SI 캡션 확인 과정을 ②에 1문장 기록(가점).

---

## 2. 데이터 자산 및 전처리 규칙 (Phase 0)

### 파일 (모두 `/home/kimjisan95/ar_protac_project/`)
- `ProtacDB_MTL_CLS.csv` — 구 PROTAC-DB. **6107행×25열.** `Protac_SMILES`, `Uniprot`, `pDC50`, `Dmax`, `Degrader_Class`, **`*_mask` 컬럼**(`pDC50_mask`, `Dmax_mask`, `Degrader_Class_mask` 등). **DOI/연도 없음.**
- `protac.xlsx` — 신 PROTAC-DB 3.0. **15502행×89열.** `Smiles`, `Uniprot`, `DC50 (nM)`, `Dmax (%)`, `Article DOI`, precomputed 물성(MW/XLogP3/TPSA/HBA/HBD/RotBond/Ring). precomputed `InChI Key`는 cross-DB 매칭에 신뢰 불가.
- `jm4c02226_si_002.csv` — case study(encoding **cp949**, header **[0,1]**). **29 데이터 행(A16/B6/C7).** 3번째 칸=Skin Retention, 4/5번째=DC50/Dmax(B/C만).

### 전처리 규칙(코드로 강제, `scripts/phase0_freeze.py`)
1. **AR 필터 = `Uniprot=='P10275'` 단일 기준.** Target 정규식 금지. → 구DB 367행, 신DB 777행.
2. **InChIKey 매칭은 양 DB 모두 RDKit `MolToInchiKey` 앞 14자(골격키)로 재계산.** ±1 변동 방지 위해 골격키 단일 기준으로 동결.
3. **DC50/Dmax 파싱:** 문자열 정제 — `>`, `<=`, `>=` 부등호 분리, `/` 슬래시 다중값 분리, `N.D.`→NaN, 콤마/공백 정규화. censored 처리:
   - 분류(Active/Inactive): `'>1000'`/`'>30000'`→**Inactive**, `'<=50'`→**Active**(정보 보존).
   - 회귀(쓸 경우): censored **제외**하고 명시(점추정 금지).
4. **구DB 라벨은 `*_mask` 컬럼으로 실측 행만 선택**(pDC50/Dmax 0-imputed 제거). `Degrader_Class`는 `Degrader_Class_mask==1`만 → 507행(양성 254/음성 253, 균형). class weight 불필요(C14).
5. **case study:** `pd.read_csv(..., encoding='cp949', header=[0,1])`. '÷'→'±' 인식, 선행 float만 파싱. A/B/C 계열 **분리 보존**. A는 활성 데이터 없음 명시.
6. **시간분할(구조적):** train 후보 = 구DB(old) / test = 신DB에만 있는 **AR 298 unique(ik14)**(골격키). DOI 연도는 보조 annotation으로만 — 엄밀 연도분할 주장 안 함.

### Phase 0 산출물
- `data/processed/ar_old.csv`, `ar_new.csv`, `case_study_clean.csv`(A/B/C 라벨 + 활성유무 플래그), `leakage_report.csv`(모델별 학습셋 ∩ 평가셋 InChIKey + 누수 화합물 ID).
- **[표1]** 데이터 무결성: 행수, AR(P10275) 수, unique ik14, DC50/Dmax 결손율, 구조분할 크기, 누수 화합물 수/ID. **§0-END 동결 스크립트 출력 단일 인용.**

---

## 3. Phase별 실행계획

> 우선순위: **MUST-HAVE**(보고서 골격) / **NICE-TO-HAVE**(시간 여유 시). 상세 우선순위·일정은 §6.

### Phase 1 — Chemical Space & 물성 (MUST-HAVE)
- 데이터: 구DB AR, 신DB AR, case study 29개 3집합.
- 기술자: **xlsx precomputed 물성 1차 사용**(MW, XLogP3, TPSA, HBA/HBD, RotBond, Ring). case study는 RDKit 재계산. Morgan FP(r=2, 2048bit).
- 분석:
  - PCA + t-SNE(또는 UMAP)로 3집합 + case study A/B/C 계열을 색/마커 구분 시각화.
  - **bRo5 overlay**: Ro5 선(MW=500, TPSA=140) + bRo5 상한(MW≈1000, TPSA≈250, **Doak 2014 인용**). AR PROTAC이 Ro5를 명백히 위반(MW≫500)함을 보임.
- 실행 env: `protac` (`~/anaconda3/envs/protac/bin/python`).
- **⑤ talking point:** (1) Morgan FP+PCA 분산설명률 낮으면 t-SNE 거리 비해석. (2) logP 출처(XLogP3 vs RDKit MolLogP) 차이가 하류 Kp에 전파 → logP ±1 민감도 1줄.
- 산출물: **[그림1]** chemical space, **[표1]** 데이터셋·필터·분할·누수 요약.

### Phase 2 — 분해 예측: prospective 평가 (MUST-HAVE)

**모델 선택(현실화):**
- **PRIMARY = PROTAC-STAN** (`~/PROTAC-STAN`, env `protac`). 유일하게 **SMILES+물성만으로 29개 AR을 랭킹 가능.** 가중치 `saved_models/protac-stan.pt` 존재, ESM 임베딩(`data/custom/esm_s_map.pkl`)에 AR(P10275)/CRBN(Q96SW2) 포함. **단, torch 2.8 비호환으로 패치 전에는 미가동(C17).**
- **SECONDARY = DegradeMaster (NICE-TO-HAVE, ⑤ '모델 간 불일치' 논증용)** [사용자 결정]. 양성=OR 라벨이라 STAN(양성=AND)과 대비가 ⑤ 논증으로 가치 큼.
  - **[feasibility gate 판정완료 — `prepare_data.py:GraphData.process()` 코드 확인]** 화합물마다 **① target_pocket PDB(AR) + 포지셔닝된 warhead mol2, ② ligase_pocket PDB(CRBN) + E3 ligand mol2, ③ PROTAC 3D mol2(EGNN은 좌표 필수), ④ 사전계산 `features/{protac,target,e3}_feature.npy`**를 요구하고, PyMOL `cmd`로 '리간드 X Å 이내 잔기' 포켓을 잘라냄 → **본질적으로 3원(ternary) 구조 입력 필요.** 또한 `prepare_data.py`는 모듈 임포트 시 `from pymol import cmd`·`from openbabel import pybel` 필요(현 `protac` env 미설치).
  - **결론 [사용자 결정: 미적용 문서화만]:** **PROTAC-8K(Zenodo 14728925) '데이터 확보'와 무관하게, 우리 29개 적용은 AR 3원구조 부재로 2일 내 비현실적**(PROTAC별 도킹 포즈·AR/CRBN 포켓·3D conformer 생성 필요). → **설치·데모 없이** "**구조기반 분해예측기(DegradeMaster/DeepPROTACs/AiPROTAC)는 AR 3원구조 부재로 우리 29개에 미적용, PROTAC-STAN만 SMILES+서열로 가동**"을 **⑤의 정직한 공통 한계로 서술**. pymol/openbabel 설치·`case_study.py` 데모는 시간 예산상 생략(⑤ 논증에 데모 불필요).
- **차순위(동일 구조-입력 장벽, 모두 포켓 필요로 확정):** DeepPROTACs(웹서버 `bailab`, mol2 포켓 4개 필요)·**AiPROTAC(`TargetGraphs_8A.bin`/`LigaseGraphs_8A.bin` = 리간드 8Å 포켓 그래프 → 구조기반 확정**; README가 PyMOL 포켓 추출 명시; py3.7/dgl env, `.pth` 미배포). 셋 다 막히면 **공통 장벽("구조기반 분해예측기 3종은 AR 3원구조 없이는 적용 불가, PROTAC-STAN만 SMILES+서열로 가동")** 자체가 ⑤의 강한 정직한 발견.
- **포켓-프리 대안 후보(웹 검증 대기 — API 과부하로 미완):** 리간드 기반 예측기 **PROTAC-Degradation-Predictor (Ribes et al., 2024, JCIM/arXiv)** — POI 리간드·E3 리간드·링커 Morgan FP + 셀라인/E3 정체성 → XGBoost/MLP로 active/inactive(DC50/Dmax 기반). **포켓 불필요**, STAN(서열 트랜스포머)과 방법이 달라 ⑤ '모델 간 불일치' 논증에 이상적. **GitHub/가중치 공개 여부를 웹으로 확인 후 채택 판정**(adopt 후보). 참고: PrePROTAC(ESM 기반)은 **타겟-수준** 분해가능성(화합물 랭킹 불가, DegradoMap과 동류) → ①타겟 정당화용만.
- **DegradoMap** = 화합물 랭킹 아님. ①타겟 정당화 전용(Phase 1.5).

**PROTAC-STAN 실행 절차(C17/C18/C19 반영):**
1. **STAN 스모크 테스트(최우선, 전처리 전).** **[결정] monkeypatch 래퍼 `scripts/stan_patch.py`**: 원본 미수정. 래퍼에서 `import torch; _orig=torch.load; torch.load=lambda *a,**k:{**k,'weights_only':False} and _orig(*a,**{**k,'weights_only':False})`(또는 `functools.partial`)로 `torch.load` 기본값을 덮어쓴 뒤 `inference.py`를 import/호출 — `data.py:181`/`data_loader.py:47/49/52`/`inference.py:69`가 모두 패치된 `torch.load`를 사용하게 됨. 데모 `custom.csv`로 1회 추론 정상 출력 확인.
2. **입력 파일 = `data/custom/custom.csv`**(29행, **`raw.csv` 아님**). 컬럼 = `Uniprot=P10275`, `E3 ligase Uniprot=Q96SW2`(**전 29행 CRBN, C19**), `Smiles`(29개) + **9개 RDKit 물성**(`Molecular Weight, Exact Mass, XLogP3(=MolLogP), Heavy Atom Count, Ring Count, Hydrogen Bond Acceptor Count, Hydrogen Bond Donor Count, Rotatable Bond Count, Topological Polar Surface Area` — `data.py` `columns`와 1자 일치).
3. 재실행 전 **`data/custom/processed/custom` 삭제**(process()가 캐시).
4. `inference.py`를 **class-1 확률 출력**으로 패치(`exp(F.log_softmax)`의 class-1 성분 덤프) — argmax(`torch.max(...,dim=1)`)만으로는 랭킹 불가.
5. 실행: `cd ~/PROTAC-STAN && ~/anaconda3/envs/protac/bin/python inference.py --root data/custom --name custom`.

**평가 프레이밍(엄밀):**
- 사전학습 가중치 **동결**, unseen 셋에서 **prospective 추론**만. **AR-only 재학습 안 함**(사유=AR 관측 양성 25개, C4).
- **누수 점검 의무:** `data/PROTAC-fine/train_compound_smiles.csv`의 ik14 집합(968) 추출 → 평가셋 교집합 제거 → **leakage-free 부분집합 + 전체** 둘 다 보고. **누수 = B3, B5**(명시).
- **표준화 코호트 고정(C19/⑤):** `data.py:209`가 입력 코호트 내부에서 per-column z-score를 계산하므로, 동일 SMILES도 추론 코호트(29 전체 vs leakage-free 부분집합 vs N)가 바뀌면 mean/std가 바뀌어 점수·랭킹이 달라짐 → **부분집합과 전체의 STAN 확률은 직접 비교 불가.** 절차로 못박음: **항상 동일한 29행 `custom.csv`로 1회만 추론**하고, 누수 표시는 **사후 마스킹으로만**(부분집합 재추론 금지). 보고서에 이 사실 명시.
- **지표(소표본·양성부재 대응):** 주지표 = **class-1 활성확률 vs 실측 DC50/Dmax의 Spearman 순위상관**. 이진분류 지표(AUROC/AUPR/F1/MCC)는 **양성 클래스가 ≥2개 확보되는 컷에서만 보조 산출**(양성 수 명기). **RMSE/R² 안 씀(모델이 분류기).** AUPR은 정보용(불균형 대응 근거 아님).
- **라벨/임계값 (각 모델 자기 정의 준수 — 통일 컷 단일 채점 금지, C20):**
  - **PROTAC-STAN(평가 대상)** 양성 = `DC50<100nM AND Dmax≥80%` — 코드 `clean.ipynb`(cell 16-17)에서 DC50 strict `<100`, **Dmax는 inclusive `≥80`**(+ 단일값 완화·%-분해 백필). 사용자 가설 `Dmax>80`은 부등호만 `≥`로 정정.
  - 참고(타 모델): DeepPROTACs/AiPROTAC = `DC50<100 AND Dmax>80%`(strict AND; 이 AND 관례의 기원=DeepPROTACs 2022). **DegradeMaster/DegradoMap = `DC50<100 OR Dmax≥80%`(양성=OR, 음성만 AND)** → 같은 숫자라도 AND 위치가 반대라 양성집합이 다름.
  - 우리 GT(`Degrader_Class`, mask==1 507행) = `pDC50>7(=DC50<100) AND Dmax>0.80`(strict AND). Dmax 부등호(STAN `≥` vs GT/DeepPROTACs `>`)는 경계행에서만 충돌 → [표2] 각주 표기.
  - **[CRITICAL] jm4c B/C 13개: STAN AND 기준(`DC50<100 AND Dmax≥80`) 양성 = 0개(최대 Dmax=76%) → AUROC/AUPR/F1 산출 불가.** `DC50<100` 단독이어도 양성 2개(**B3[누수]/C5**)뿐 → 분류지표는 조건부, 주지표는 순위상관.
  - censored `>1000`=Inactive 일관. Dmax 스케일 통일(구DB fraction 0~1.12 vs jm4c % → 0.80 vs 80%). 민감도 컷(50/100/1000nM)은 각각 양성수 명기.
- **소표본 처리:** 검증셋 실태표(retention n=29, DC50+Dmax n=13, 수치형 DC50 n=5, censored 다수)를 정직히 보고. 모든 지표에 **bootstrap 95% CI(2000 resample)**. 결론은 **"정성 경향성/가설 생성"**.
- **STAN concordance 검증:** **leakage-free 주 앵커 = C6(199.5nM), C1(103.87nM), C5(70.85nM)**. **B3(8.97nM)·B5는 leaked → "memorization 가능, 참고용"으로만 분리 표기**(검증 근거 아님). small N → RMSE 아님, 정성 일치만.
- **⑤ talking point:** (1) STAN은 binary → **DC50 nM 재현 불가**, 활성확률 랭킹뿐. (2) **코호트 의존 z-score** → 부분집합/전체 비교 불가·점수 변동(고정 절차 명시). (3) **단일 E3(CRBN) 시리즈 → E3 채널 상수, 랭킹은 PROTAC 그래프+물성만으로 결정**(C19). (4) 사전학습 분포 밖 bRo5/AR 외삽. (5) **leakage 2건(B3, B5)**과 n 작음. (6) torch 비호환 패치(재현성, C17). (7) 모델군마다 **AND 위치 반대**(STAN/DeepPROTACs=양성 AND vs DegradeMaster/DegradoMap=양성 OR) → 자기 기준 평가; STAN 자기 기준 적용 시 jm4c 양성 0개라 분류지표는 조건부·주지표는 순위상관(C20).
- 산출물: **[그림2]** STAN 활성확률 랭킹 vs 실측(B/C DC50·Dmax; **C6/C1 강조 + B3/B5 leaked 마킹**) + leakage-free 표시, **[표2]** 29개 물성+STAN확률+실측+누수플래그 요약.

### Phase 1.5 — DegradoMap: 타겟 druggability 정성 점검 (NICE-TO-HAVE, ①에 사용)
- 역할: **AR이 PROTAC 분해 친화 타겟인가** 정성 근거. **per-compound 점수 절대 금지.**
- 실행: **`~/DegradoMap/results/*.json`(ablation/baseline/e3_evaluation/case_studies) 정성 인용으로만 확정.** **신규 추론은 시도 안 함**(타겟 정당화 ①은 기존 결과로 충분; STAN과 동일하게 torch/환경 비호환 가능성 큼). 굳이 한다면 Day2 NICE 슬롯에서 `weights_only` 패치 가능성 먼저 점검.
- **⑤ talking point:** target-unseen **AUROC ~0.60(6-seed), gradient boosting 대비 p=0.556** → 약한 신뢰. PROTAC 구조 미입력이라 29개 변별 불가.
- env: `protac`.

### Phase 3 — 피부 투과/잔류 (MUST-HAVE, 재설계)

**3-1. 전신투과 proxy (Potts-Guy 1종) [MUST]:**
- 식: `logKp(cm/h) = 0.71·logP − 0.0061·MW − 6.3` (Potts-Guy 1992, Flynn 데이터셋). RDKit MolLogP + MolWt로 직접 구현(투명·재현 가능).
- **AD 외삽 정량화:** B/C 13개 전부 MW 725~792로 적합 상한(~750) 경계 밖임을 표로. "MW항 지배로 모든 PROTAC이 균일하게 낮은 Kp로 압축 → 변별력 0"을 **예상 결과로 미리 명시**.
- **절대 Kp 값 신뢰구간 없이 제시 금지.** 상대 랭킹/정성 용도로만.

**3-2. 비교 QSPR(선택):** "Mitragotri 모델" 표기 **철회**. 필요 시 **출처 있는** Magnusson MW-only 또는 "단순 휴리스틱"으로 명시 강등. SETUP 임의식 사용 금지.

**3-3. 2D 분리 시각화(§1-1) [MUST]:** x=예측 logKp(전신투과), y=실측 retention(피부저류). **B/C 13개는 활성색**, A 16개는 retention 분포로 별도 표시. 이상적 사분면 시각화. **직접 상관계수를 주 결과로 쓰지 않음**; 보더라도 Spearman + bootstrap CI + LOO.

**3-4. 카멜레온성(NICE-TO-HAVE로 강등, §C10·C17):**
- **MUST는 저비용 2D 카멜레온(TPSA vs RotBond, 2D proxy 명시)만.** 3D는 시간 여유 시.
- 3D 수행 시: ETKDG **conformer 분자당 10~20개 상한**, **앵커 5~6개(C6/C1/C5/B6 등) 또는 B/C 13개로 범위 축소** → min3D_PSA, ΔPSA(=TPSA−min3D_PSA), IMHB 개수, Fsp3, PBF, Rg. MW 700~1100·RotBond 20+ 거대 PROTAC은 ETKDG 임베딩 실패·분자당 수십초~수분이므로 29개 전수 신규 파이프라인을 반나절에 완성 불가.
- 설계전략 표현: "TPSA 감소"가 아니라 **"IMHB 배치로 conformer별 3D EPSA 저하(TPSA 동일해도)"**.
- **ETKDG 임베딩 실패율 자체를 ⑤ 소재(거대분자 conformer 한계)로 프레이밍.**

**3-5. PROTAC-TS Caco-2(NICE-TO-HAVE):** `protac_ts` env. xlsx Caco-2 Papp 부분집합 → `make_feature.py`(Morgan-count) + `make_model.py`(TabPFN) → 29개 예측. **§1-2 3중 구분 명시**, leakage 점검.

- **⑤ talking points:** (1) Potts-Guy AD 외삽 → 정량 신뢰 불가·검증 부재. (2) 모낭(follicular) 경로·유한용량·제형 효과 무시(두피 탈모제 핵심 메커니즘 누락) → ⑥. (3) 세포투과≠피부투과≠피부잔류. (4) 3D 기술자는 conformer 샘플링·force field 의존 → 정성적. (5) n=13(또는 더 작음) → 모든 상관은 가설.
- 산출물: **[그림3]** 2D 분리 시각화(retention vs logKp, A/B/C 구분)**(기본 확정)**, [그림3 부록] 카멜레온성(2D, 여유 시 3D), **[표3]** AD 외삽 정량(MW>750·TPSA>140 화합물 수).

### Phase 4 — Boltz 3원 복합체 (NICE-TO-HAVE, 명시적 OPTIONAL, 런타임 리스크)
- env: `protac_boltz` (boltz 2.2.1).
- **AR 서열:** SETUP의 polyQ/polyP 단편은 오류. **PDB 2AM9 SEQRES에서 실제 존재하는 잔기를 코드로 추출**해 AR-LBD 서열로 사용(잔기 범위 하드코딩 금지; UniProt P10275 LBD 코어는 출처마다 ±수 잔기 차이가 있음). 전장 사용 안 함. 2AM9는 참조 구조.
- E3 파트너: **전 화합물 CRBN(C19)이므로 CRBN(±DDB1)** 단순화. 풀 복합체 비용↑ → **1 케이스만(C6, leakage-free 리드)**.
- **리스크:** AR-LBD(~250aa)+E3+PROTAC(원자 수백~1000+) 공동접힘은 GPU 수십분~시간, 수렴 보장 없음. **보고서 핵심 그림을 Boltz에 의존시키지 말 것.**
- **시간 게이트:** **모든 MUST 완료 + 보고서 초안 작성 후에만 1런.**
- **⑤ talking point:** 단일 시드·pLDDT/ipTM만으로 결합 단정 금지. 미수렴/실패도 결과로 서술.
- 산출물: **[그림4(선택)]** Boltz 포즈 + pLDDT/ipTM.

### Phase 5 — 신규 분자 설계: hypothesis generation (MUST-HAVE, 순환논증 경고)
- **재프레이밍:** "in-silico 검증"이 아니라 **"가설 생성(hypothesis generation)"**. Phase 2/3에서 한계가 드러난 모델로 자기 설계물을 "검증"하면 **순환논증** → 채점 ⑤와 정면 충돌.
- 설계 근거는 **모델 점수가 아니라 Phase 1 SAR/물성 규칙**(case study 29개는 전부 **CRBN/글루타리미드 + darolutamide류 워헤드** 기반(C19), 적정 MW/TPSA, IMHB 가능 배치)에 둔다. **VHL/하이드록시프롤린은 case study 근거가 아님**(0/29). 만약 'VHL 선호'를 주장하려면 **별도 데이터(구DB P10275 내 VHL 성공사례)**를 제시.
- 신규 SMILES **3개**(여유 시 5개) 생성(수동 medchem 편집 또는 PROTAC-TS ChemTSv2 링커설계).
- 모델 출력은 **"약한 사전확률(weak prior)"**로만 해석하고 반드시:
  - (a) **applicability domain 점검**: 신규 분자가 학습 화학공간 내부인지 Tanimoto 최근접 거리.
  - (b) 모델 불확실성/앙상블 분산 동반.
  - (c) **"실험 검증 필요"**를 결론으로 명시.
- 산출물: **[표4]** 신규 3개 SMILES + 설계근거(SAR) + STAN 사전확률 + AD 거리.

---

## 4. ⑤ 비판적 고찰 talking points 마스터 (보고서 1.5쪽 핵심)

- **데이터:** AR 필터 P10275(Target 정규식이면 1193 부풀림). 구DB `Degrader_Class`는 mask 적용 시 **507행 50:50 균형**(4.2%는 결손율 91.7%의 오해석). case study 3번째 칸=Skin Retention(SI 확인). censored DC50 처리.
- **Phase 1:** Morgan FP+PCA 분산설명률 낮으면 t-SNE 거리 비해석. logP 출처 차이 전파.
- **Phase 2:** (1) STAN은 분류기 → DC50 nM 재현 불가, 활성확률 랭킹뿐. (2) **코호트 의존 z-score** → 부분집합/전체 직접 비교 불가, 고정 절차 명시. (3) **단일 E3(CRBN) → E3 채널 상수, 랭킹=PROTAC 그래프+물성**(C19). (4) 사전학습 분포 밖 bRo5/AR 외삽. (5) **leakage 2건(B3, B5) 정량**, 특히 이전 핵심 앵커 B3가 누수. (6) n 작음 → bootstrap CI. (7) torch 비호환 5줄 패치(재현성). (8) 포켓 의존 3모델 적용 불가(AR 3원구조 부재) = 정직한 공통 한계. (9) 모델별 라벨 정의 상이(AND 위치 반대) → 자기 기준 평가, STAN AND 기준 jm4c 양성 0개 → 분류지표 대신 순위상관(C20).
- **Phase 1.5:** DegradoMap AUROC ~0.60, p=0.556, 화합물 변별 불가. 신규 추론 안 함(기존 결과 인용).
- **Phase 3:** (1) Potts-Guy AD 외삽 → 정량 신뢰 불가. (2) Skin Retention≠Kp(endpoint 오류 수정 과정). (3) 모낭 경로/유한용량 누락. (4) 세포투과≠피부투과≠잔류. (5) 카멜레온성은 conformer 의존 정성 proxy + 거대분자 ETKDG 실패율.
- **Phase 4:** 단일 시드·pLDDT만으로 결합 단정 금지. 미수렴도 결과.
- **Phase 5:** 순환논증 회피 — 모델 점수는 weak prior, 실험 검증 필요. 설계 근거는 CRBN/글루타리미드 SAR(VHL 아님).

---

## 5. 보고서 그림/표 manifest (10쪽 제약, ④는 6~7개 / ⑤에 지면 양보)

모든 산출물은 **생성 스크립트 절대경로 → 출력 경로**를 고정한다(모든 경로 prefix `/home/kimjisan95/ar_protac_project/`).

| 산출물 | 내용/인코딩 사양 | 스크립트 → 출력 | Phase | 위치(보고서) | 폴백 |
|---|---|---|---|---|---|
| **[표1]** | 행=[구DB AR, 신DB AR, case study]×열=[행수, P10275수, unique ik14, DC50결손%, Dmax결손%, 구조분할 크기, 누수수/ID] | `scripts/phase0_freeze.py` → `outputs/table1.csv` | Phase 0/2 | ② 데이터 (1쪽) | 필수 |
| **[그림1]** | Chemical space PCA(또는 t-SNE). x/y=PC1/PC2, 색=데이터셋(3), 마커=A/B/C, Ro5(MW=500,TPSA=140)·bRo5(MW=1000,TPSA=250) 점선 overlay | `scripts/fig1_chemspace.py` → `outputs/fig1.png` | Phase 1 | ④ (0.5쪽) | 필수 |
| **[표2]** | case study 29개: ID, 계열, MW, TPSA, RotBond, STAN class-1 확률, 실측 DC50/Dmax, **leaked 플래그(B3/B5)** | `scripts/phase2_stan.py` → `outputs/table2.csv` | Phase 2 | ④ (0.5쪽) | 필수 |
| **[그림2]** | STAN 확률 랭킹 vs 실측. x=화합물(확률순), y=class-1 확률, 색=실측 활성(DC50<100 등), **C6/C1 강조 라벨 + B3/B5 leaked 마킹(다른 마커)** | `scripts/fig2_ranking.py` → `outputs/fig2.png` | Phase 2 | ④ (0.5쪽) | 필수 |
| **[그림3]** | **2D 분리(기본 확정).** x=예측 logKp, y=실측 retention(%), 색=계열(A/B/C), B/C는 활성 테두리. 이상 사분면 음영 | `scripts/fig3_separation.py` → `outputs/fig3.png` | Phase 3 | ④ (0.5쪽) | 필수(2D 분리 고정) |
| [그림3-부록] | 카멜레온성: 2D(TPSA vs RotBond) 기본, 3D(ΔPSA/Fsp3, 활성색)는 여유 시 | `scripts/fig3b_chameleon.py` → `outputs/fig3b.png` | Phase 3 | 부록 | 2D만이라도 |
| **[표3]** | Potts-Guy AD 외삽: 화합물별 logKp, MW>750 플래그, TPSA>140 플래그, AD 밖 비율 | `scripts/phase3_pottsguy.py` → `outputs/table3.csv` | Phase 3 | ④ 또는 ⑤ (0.3쪽) | 필수 |
| **[표4]** | 신규 3개 SMILES + 설계근거(SAR) + STAN 사전확률 + Tanimoto AD 거리 | `scripts/phase5_design.py` → `outputs/table4.csv` | Phase 5 | ④ (0.5쪽) | 필수 |
| **[그림4]** *(선택)* | Boltz 포즈 + pLDDT/ipTM | `scripts/fig4_boltz.py` → `outputs/fig4.png` | Phase 4 | ④(선택) | 생략 가능 |

**그림 동결 체크포인트:** Day2 집필 시작 전, 위 필수 산출물(표1·2·3·4, 그림1·2·3)이 `outputs/`에 존재하는지 확인 후 포맷 동결. 보고서 6섹션 ↔ 채점기준 1:1 매핑, **⑤를 가장 두껍게(1.5쪽)**. 나머지는 부록/생략.

---

## 6. 우선순위 및 일정 (남은 2일)

### MUST-HAVE (이것만으로 완결된 보고서 성립)
Phase 0(데이터 동결·AR 필터·누수표) → Phase 1(chemical space + bRo5) → Phase 2(PROTAC-STAN 추론 + leakage 표 + bootstrap CI) → Phase 3(Potts-Guy + AD 외삽 + 2D 분리 + **저비용 2D 카멜레온**) → Phase 5(신규 3개) → 보고서 집필.

### NICE-TO-HAVE (시간 여유 시만)
3D 카멜레온성(앵커 5~6개), PROTAC-TS Caco-2, DegradoMap 신규 추론, Boltz 1런, 두 번째 분해모델.

### Day 1 (6/17 화) — 모델 파이프라인·데이터
- **오전(최우선):** **PROTAC-STAN 스모크 테스트 + torch.load 5곳 `weights_only=False` 패치(C17)** — 데모 `custom.csv`로 추론이 도는지 먼저 확인(이게 막히면 전체 막힘). 병행: Phase 0 전처리 코드(`scripts/phase0_freeze.py`: AR 필터 P10275, 양DB RDKit ik14 재계산, DC50/Dmax 파싱, mask, leakage 표). 산출물 CSV 동결.
- **오후:** `custom.csv` 29행 작성(**전 행 CRBN=Q96SW2, 9개 RDKit 물성, C18/C19**) + `inference.py` class-1 확률 패치 → 추론 실행(29행 1회만, 코호트 고정). Phase 1 chemical space + bRo5([그림1]).
- **저녁:** Phase 2 leakage 점검(B3/B5 확인) + bootstrap CI 지표 + [그림2]/[표2]. 누수 사후 마스킹.

### Day 2 (6/18 수) — 투과성·설계·보고서
- **오전:** Phase 3 MUST(Potts-Guy 직접 구현, AD 외삽 표, 2D 분리 시각화, **저비용 2D 카멜레온**) → [그림3]/[표3]. (3D 카멜레온성은 시간 남으면 앵커 5~6개로만.)
- **이른 오후:** Phase 5 신규 3개 설계(CRBN SAR 근거) + STAN 사전확률 + AD 거리 → [표4].
- **오후 절반 이후~저녁: 그림 동결 체크포인트 → 보고서 집필 못박기**(6섹션, ⑤ 최우선). manifest 그림/표 삽입.

### 6/19 (목/금 자정 전) — 예비/제출
- 누락 보강, 문헌 인용 정리(Potts-Guy 1992/Flynn, Doak 2014, Mitragotri 철회 명시, jm4c02226, PROTAC-STAN/DegradoMap 논문), (게이트 통과 시) Boltz 1런, 최종 검토·제출.

---

## 7. 위험요소 표

| 위험 | 가능성 | 영향 | 대응(실패를 결과로 전환) |
|---|---|---|---|
| **PROTAC-STAN torch 2.8 비호환(C17)** | **확실(패치 전)** | 높음 | **Day1 오전 1순위.** `torch.load` 5곳 `weights_only=False` 패치(검증 완료). 패치 사실을 ②/⑤(재현성)에 기록 |
| **E3 라벨 오류로 입력 오염(C19)** | (방지됨) | 높음 | **전 29행 CRBN=Q96SW2** 고정(글루타리미드 29/29 실측). VHL 매핑 금지 |
| **입력 파일명/스키마 오류(C18)** | (방지됨) | 높음 | `custom.csv`(raw.csv 아님) + 9개 물성 컬럼명 1자 일치 |
| 포켓 의존 모델 적용 불가 | 확실 | 중 | 시도 안 함. "AR 3원구조 부재로 구조기반 분해예측기 적용 불가"를 ⑤ 핵심 발견으로 |
| Potts-Guy가 PROTAC을 균일 저Kp로 압축(변별력 0) | 높음 | 중 | **예상 결과로 미리 명시.** "소분자 QSPR의 bRo5 외삽 무력화"가 ⑤ 소재 |
| retention vs Kp 상관 낮음/음 | 높음 | 낮음 | **이것이 정답.** "애초에 다른 물리량(잔류≠통과)"으로 프레이밍 |
| 누수(B3/B5)로 concordance 무효 | 확실 | 중 | leakage-free 앵커(C6/C1/C5) 사용, B3/B5는 누수 사례로 ⑤ 활용 |
| 3D 카멜레온 ETKDG 비용 초과 | 중 | 낮음 | NICE-TO-HAVE 강등. MUST는 2D 카멜레온. ETKDG 실패율도 ⑤ 소재 |
| n 작아 지표 불안정 | 확실 | 중 | bootstrap CI + LOO 의무. "정량 검증 아닌 가설 생성"으로 톤다운 |
| Boltz 미수렴/런타임 초과 | 중 | 낮음(optional) | MUST+초안 후에만 1런, 핵심 그림 비의존. 미수렴도 ⑤ 소재 |
| 시간 부족 | 중 | 높음 | NICE 전부 포기 가능. Day2 오후 절반은 집필 고정. 그림 동결 체크포인트 |

---

## 8. 실행 환경 빠른참조
- 메인 분석: `~/anaconda3/envs/protac/bin/python` (py3.10, **torch 2.8.0+cu128**, rdkit 2026.3.1, torch-geometric, sklearn, umap, ETKDG). → Phase 0/1/2/3/5.
- PROTAC-TS: `protac_ts` (py3.11, chemtsv2, tabpfn). → Phase 3-5(선택).
- Boltz: `protac_boltz` (py3.10, boltz 2.2.1). → Phase 4(선택).
- 주요 경로:
  - `~/PROTAC-STAN`: `inference.py`(argmax→class-1 확률 패치), `data.py`(L179 `{name}.csv`, L181 torch.load 패치, L206/209 z-score, `columns` 9물성), `data_loader.py`(L47/49/52 torch.load 패치), `saved_models/protac-stan.pt`, **`data/custom/custom.csv`(작성 대상)**, `data/PROTAC-fine/train_compound_smiles.csv`(누수 점검 968 ik14).
  - `~/DegradoMap/results/*.json`(정성 인용), `~/PROTAC-TS`(make_feature.py/make_model.py).
- AiPROTAC/DeepPROTACs/DegradeMaster: 입력(포켓/도킹/.pth/py3.7) 미충족 → 적용 불가, ⑤ 한계 인용용.
