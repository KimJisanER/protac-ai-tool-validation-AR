# AR PROTAC 탈모 치료제 프로젝트 실행 계획

**마감: 2026-06-19 (금) 자정**  
**남은 시간: 약 2일**  
**사용 conda 환경: `protac` (rdkit, pandas, torch, sklearn 등) / `protac_boltz` (boltz)**  
**Python 경로: `~/anaconda3/envs/protac/bin/python`**

---

## 데이터 구조 정리

| 파일 | 설명 | 규모 |
|---|---|---|
| `ProtacDB_MTL_CLS.csv` | **구버전** PROTAC-DB (시간 분할 기준 Train set) | 6,107개 |
| `protac.xlsx` | **신버전** PROTAC-DB 3.0 (전체 데이터, DOI 컬럼 포함) | 15,502개 |
| `jm4c02226_si_002.csv` | case study 논문 (AGA 치료 AR PROTAC 30개, 실험값 有) | 30개 |

---

## Phase 0: 데이터 전처리 및 Time Split 구성 ✦ 우선순위 최고

### 0-1. 데이터 로딩 및 AR PROTAC 필터링
- `protac.xlsx`에서 `Target == 'AR'` (또는 Uniprot `P10275`) 항목 추출
- `ProtacDB_MTL_CLS.csv`에서도 AR 타겟 항목 추출
- AR PROTAC 개수 확인 (기대치: 수십~수백 개)

### 0-2. Time Split 구성
- `protac.xlsx`의 `DOI` 컬럼으로 논문 출판 연도 수집
  - `requests` / `crossref API` 또는 `habanero` 라이브러리 활용
  - 연도 기준 정렬 → 구버전 ProtacDB에 없는 신규 항목 = **Test set**
- `jm4c02226_si_002.csv` 30개 중 `ProtacDB_MTL_CLS.csv`에 없는 화합물 추출 (SMILES 기준 InChIKey 매칭)
  - 이 화합물들 = **skin retention 실측값이 있는 핵심 검증 셋**

### 0-3. 데이터 정제
- SMILES 유효성 검증 (RDKit `MolFromSmiles`)
- `DC50` 값 전처리: `'>1000'`, `'> 1000'` 등 → `NaN` 또는 상한값 처리
- 피처 컬럼 정규화

---

## Phase 1: Chemical Space 탐색

### 1-1. 분자 기술자 계산 (RDKit)
다음 기술자를 전체 AR PROTAC에 대해 계산:
- MW (분자량), TPSA, cLogP, HBD, HBA, RotBonds, ArRings
- 링커 구조 유형 (PEG / alkyl / aromatic 분류)
- E3 리가아제 리간드 유형 (VHL / CRBN / MDM2 등 분류)

### 1-2. Chemical Space 시각화
- **PCA** + **t-SNE** (Morgan Fingerprint 기반)
- 3개 집합을 다른 색으로 표시:
  1. 전체 PROTAC-DB
  2. AR 타겟 PROTAC (필터링)
  3. jm4c02226 case study 30개 (마커 강조)
- bRo5 경계선 (MW > 700, TPSA > 140) overlay

### 1-3. ADMET 프로파일링
- **SwissADME** 또는 **ADMETlab 3.0** (웹 서버 활용 가능)
- 또는 RDKit + sklearn으로 기본 ADMET 기술자 계산
- jm4c02226 30개의 물성 분포를 분해 효율별로 비교

---

## Phase 2: 분해 효율 AI 모델 평가 ✦ 핵심

### 2-1. 모델 선택 (시간 고려, 2~3개 선택)
우선순위 순:

| 모델 | 레포 | 특징 | 난이도 |
|---|---|---|---|
| **DeepPROTACs** | `fenglei104/DeepPROTACs` | GNN 기반, 구조 정보 활용 | 중 |
| **DegradeMaster** | `ABILiLab/DegradeMaster` | 최신, multi-task | 중 |
| **PROTAC-STAN** | `PROTACs/PROTAC-STAN` | Transformer 기반 | 중~고 |
| AiPROTAC | `LiZhang30/AiPROTAC` | 초기 모델 | 낮 |
| DegradoMap | `bryanc5864/DegradoMap` | 시각화 중심 | 낮 |

### 2-2. Time Split 평가 방법
- **Train**: ProtacDB_MTL_CLS.csv 내 AR PROTAC (구버전 데이터)
- **Test**: protac.xlsx 중 ProtacDB_MTL_CLS에 없는 신규 AR PROTAC
- **검증셋**: jm4c02226_si_002.csv 30개 (실험값 비교)

평가 지표:
- DC50 예측: RMSE, R², Pearson correlation
- Active/Inactive 분류: AUROC, AUPR (DC50 < 100nM = Active 등)

### 2-3. 비판적 고찰 준비
- 모델별 학습 데이터 포함 여부 확인 (data leakage 점검)
- bRo5 화합물에서 GNN/Transformer의 한계 분석
- 예측 실패 케이스 분석 (>1000nM 예측 vs 실측 불일치)

---

## Phase 3: 피부 투과성 예측 ✦ 핵심 차별점

### 3-1. Potts-Guy 다중선형회귀 모델 (직접 구현)
```
log(Kp) = 0.71 × log P - 0.0061 × MW - 6.3
```
- RDKit의 `MolLogP`, `MolWt`로 즉시 계산 가능
- 30개 jm4c02226 화합물에 적용
- 예측 Kp vs 실측 Skin Retention Rate 상관관계 분석

### 3-2. Mitragotri 결정론적 모델 (직접 구현)
- 분자 크기(반지름), 지질친화도 기반 모델
- TPSA 보정 항 추가 고려

### 3-3. PROTAC-TS (ycu-iil/PROTAC-TS)
- 세포막 투과성 ML 예측
- 레포 설치 후 `protac` conda 환경에서 실행
- 예측값 vs 실측 skin retention 비교

### 3-4. 분자 카멜레온성(Chameleonicity) 분석
- TPSA vs 분자 유연성(RotBonds) 산점도
- IMHB(Intramolecular H-bond) 가능성 지표
- 30개 화합물을 피부 잔류율 기준 색상 구분

---

## Phase 4: Boltz 공동접힘 예측 (선택, 시간 여유 시)

### 4-1. 대상
- **C6** (가장 좋은 PK 데이터, Skin-to-plasma ratio 우수)
- **B3** (DC50 = 8.97 nM, Dmax = 76%)

### 4-2. 실행
- `protac_boltz` 환경 사용 (`~/anaconda3/envs/protac_boltz/bin/python`)
- AR (PDB: 2AM9 또는 AlphaFold) + E3 리가아제(VHL/CRBN) + PROTAC SMILES → 3원 복합체 예측
- 결합 포즈 시각화 (PyMol 또는 py3Dmol)

### 4-3. 평가
- 예측 구조의 pLDDT 점수 확인
- 링커가 허용하는 구형화(Globularity) 정도 분석
- 도킹 점수와 DC50 상관관계

---

## Phase 5: 신규 분자 설계

### 5-1. SAR 인사이트 도출
Phase 1~3 결과를 통합하여:
- DC50 < 100nM **AND** Skin Retention > 1% → 우수 화합물 패턴 추출
- 링커 길이/강성, E3 리간드 유형, clogP 범위 정리

### 5-2. 신규 SMILES 제안
- 기존 best 화합물(C6 계열) 구조 기반
- 링커 변형 (PEG → rigid aromatic, alkyl 길이 조절)
- 피부 투과성 향상을 위한 TPSA 감소 전략 (IMHB 형성 유도)
- 3~5개 후보 SMILES 제시

### 5-3. 후보 인실리코 검증
- Phase 3 모델로 예측 피부 투과성 계산
- Phase 2 모델로 예측 DC50 계산
- SwissADME / ADMETlab으로 ADMET 프로파일 확인

---

## 보고서 구조 (권장 10페이지 이내)

| 장 | 내용 | 분량 |
|---|---|---|
| 1. 타겟 소개 및 선정 이유 | AR, AGA, PROTAC 기술, GT20029 임상 데이터 | ~1.5p |
| 2. 데이터/구조 확보 | 3개 데이터셋 설명, time split 방법론 | ~1p |
| 3. 방법 및 선택 이유 | Chemical space, AI 모델, permeability 모델 선택 근거 | ~1.5p |
| 4. 결과 | 표·그림 중심 (PCA, AUROC 그래프, 상관관계 플롯) | ~3p |
| 5. **비판적 고찰** | 각 모델의 한계, bRo5 한계, 예측 신뢰도 | ~1.5p |
| 6. 한계 및 다음 단계 | 실험 부재, 모델 일반화 문제, 향후 계획 | ~0.5p |

---

## 일정 (2일)

### Day 1 (2026-06-17, 오늘)
- [ ] Phase 0: 데이터 전처리, AR 필터링, Time Split 구성
- [ ] Phase 1: Chemical space 탐색, PCA/t-SNE, 기술자 계산
- [ ] Phase 2 준비: 모델 1~2개 설치 및 테스트 실행

### Day 2 (2026-06-18)
- [ ] Phase 2: AI 모델 평가, time split 성능 측정
- [ ] Phase 3: Potts-Guy, Mitragotri 직접 구현, PROTAC-TS 실행, skin retention 검증
- [ ] Phase 4: Boltz (시간 여유 시)
- [ ] Phase 5: 신규 분자 제안
- [ ] 보고서 작성

---

## 주요 위험 요소 및 대응

| 위험 | 대응 |
|---|---|
| AI 모델 설치 실패 | DegradeMaster / DeepPROTACs 우선, 나머지 건너뜀 |
| jm4c02226 화합물이 전부 ProtacDB에 없음 | 예상된 상황 (신규 논문). 30개 전체를 검증셋으로 활용 |
| AR PROTAC 수 부족으로 학습 어려움 | 보고서에서 data scarcity를 한계로 기술, 전체 PROTAC 대상으로 모델 평가 |
| Boltz 실행 시간 초과 | Phase 4 생략, Phase 1~3+5 집중 |
| Skin retention vs permeability 상관 낮음 | 이것 자체가 분석 결과이자 비판적 고찰 소재 |
