# AR PROTAC 프로젝트 환경 설정 및 파이프라인 요약

**마감: 2026-06-19 (금) 자정**  
**작성일: 2026-06-17**

---

## 1. 프로젝트 개요

**주제:** 안드로겐 수용체(AR) 표적 PROTAC의 탈모 치료제 가능성 평가 및 AI 기반 신규 분자 설계  
**과목:** 화학정보학 및 분자설계 (2026학년도 1학기, 대학원 글로벌혁신신약학과)  
**보고서 분량:** 권장 10쪽 이내

### 과제 요구사항 (화학정보학_기말_프로젝트_안내.docx)

| 보고서 섹션 | 내용 |
|---|---|
| ① 타겟 소개와 선정 이유 | AR 단백질, AGA 질환, PROTAC 기술, GT20029 임상 데이터 |
| ② 데이터/구조 확보 | 3개 데이터셋 출처 및 전처리 |
| ③ 방법과 선택 이유 | Chemical space, AI 분해 모델, 피부 투과 모델 선택 근거 |
| ④ 결과 | 표·그림 중심 (PCA, AUROC, 상관관계 플롯) |
| ⑤ 비판적 고찰 (가장 중요) | 모델 한계, bRo5 문제, 예측 신뢰도 discussion |
| ⑥ 한계와 다음 단계 | 부족한 점, 향후 계획 |

---

## 2. 데이터셋 현황

| 파일 | 내용 | 규모 | AR 항목 |
|---|---|---|---|
| `ProtacDB_MTL_CLS.csv` | 구버전 PROTAC-DB (Train set) | 6,107개 | 별도 필터링 필요 |
| `protac.xlsx` | 신버전 PROTAC-DB 3.0 (전체) | **15,502개** | **777개** (`Target=='AR'` 또는 `Uniprot=='P10275'`) |
| `jm4c02226_si_002.csv` | case study 논문 AR PROTAC | 30개 | 30개 (전부) |

**핵심 컬럼:**
- `protac.xlsx`: `Compound ID`, `Uniprot`, `Target`, `E3 ligase`, `Smiles`, `DC50 (nM)`, `Dmax (%)`, `Article DOI`
- `ProtacDB_MTL_CLS.csv`: `Protac_SMILES`, `Uniprot`, `Target`, `E3 ligase`, `pDC50`, `Dmax`, `Degrader_Class`
- `jm4c02226_si_002.csv` (encoding=`cp949`): `Compounds`, `Smiles`, `Skin Retention Rate(%)`, `DC50 (nM)`, `Dmax (%)`, `T1/2`, `Tmax`, `Cmax`, `AUC0-t`

---

## 3. conda 환경

| 환경 이름 | Python 경로 | 용도 |
|---|---|---|
| `protac` | `~/anaconda3/envs/protac/bin/python` | 메인 분석 (Phase 0~3, 5) |
| `protac_boltz` | `~/anaconda3/envs/protac_boltz/bin/python` | Boltz 공동접힘 예측 (Phase 4) |

---

## 4. 설치된 패키지 현황 (`protac` 환경)

### 핵심 패키지 (설치 완료)

| 패키지 | 버전 | 용도 |
|---|---|---|
| `rdkit` | 2026.3.1 | SMILES 처리, 분자 기술자 계산 (MW, TPSA, cLogP 등) |
| `pandas` | 2.3.3 | 데이터 처리 |
| `numpy` | 2.2.6 | 수치 계산 |
| `torch` | 2.8.0+cu128 | 딥러닝 모델 |
| `torch-geometric` | 2.7.0 | GNN (DegradeMaster, DeepPROTACs) |
| `scikit-learn` | 1.7.2 | ML 모델, PCA, t-SNE |
| `scipy` | 1.15.3 | 통계 분석 |
| `matplotlib` | 3.10.8 | 시각화 |
| `seaborn` | 0.13.2 | 시각화 |
| `openpyxl` | 3.1.5 | `.xlsx` 파일 읽기 |
| `requests` | 2.32.5 | HTTP 요청 |
| `networkx` | 3.4.2 | 그래프 처리 (DegradeMaster) |
| `umap-learn` | 0.5.11 | UMAP 차원 축소 (선택) |
| `tqdm` | 4.67.1 | 진행 표시 |

### 이번 세션에서 추가 설치

| 패키지 | 버전 | 설치 명령 | 용도 |
|---|---|---|---|
| `habanero` | 2.4.0 | `pip install habanero` | CrossRef API로 DOI → 출판 연도 조회 |
| `python-docx` | 1.2.0 | `pip install python-docx` | `.docx` 파일 파싱 |

### `protac_boltz` 환경 주요 패키지

| 패키지 | 버전 |
|---|---|
| `boltz` | 2.2.1 |
| `rdkit` | 2025.9.5 |
| `torch` | 2.8.0+cu128 |
| `pytorch-lightning` | 2.5.0 |

---

## 5. 클론된 AI 모델 레포지토리

### 5-1. DeepPROTACs (`~/DeepPROTACs/`)
- **출처:** `github.com/fenglei104/DeepPROTACs`
- **논문:** Li et al., Nature Communications 2022
- **방식:** GNN 기반 분해 활성 예측 (binary: Active/Inactive)
- **입력:** 리가아제 포켓, 타겟 포켓, 리간드 `.mol2` + 링커 `.smi`
- **주의:** mol2 포켓 파일 필요 → 웹서버 사용 권장 (bailab.siais.shanghaitech.edu.cn)
- **실행 환경:** 별도 `DeepPROTACs` conda env 필요 (Python 3.7, CUDA 11.1)

### 5-2. DegradeMaster (`~/DegradeMaster/`)
- **출처:** `github.com/ABILiLab/DegradeMaster`
- **논문:** Liu et al., bioRxiv 2025 (ISMB/ECCB'25)
- **방식:** EGNN 기반 multi-task 분해 예측
- **데이터:** PROTAC-8K (Zenodo: 14728925) — **미다운로드, `processed/` 비어 있음**
- **체크포인트:** `~/DegradeMaster/checkpoint/1000/`, `2000/` 존재
- **실행 명령:** `python main.py --config config/config.yml` (mode: Test)
- **실행 환경:** `protac` 환경에서 실행 가능 (networkx, torch-geometric 충족)

### 5-3. AiPROTAC (`~/AiPROTAC/`)
- **출처:** `github.com/LiZhang30/AiPROTAC`
- **특징:** AR degrader 설계에 직접 적용된 도구 (GT19 리드 화합물 도출), PROTAC-ZL 자체 데이터셋 포함
- **방식:** GNN + 지도/비지도 학습 혼합 (PROTAC-DB 2.0 기반)
- **데이터:** `data/` 폴더에 전처리 완료 데이터 포함되어 있음
- **실행 환경:** Python 3.7, DGL, CUDA 11.7 (별도 env 필요)
- **난이도:** 낮음 (초기 모델, AR 특화)

### 5-4. PROTAC-STAN (`~/PROTAC-STAN/`)
- **출처:** `github.com/PROTACs/PROTAC-STAN`
- **논문:** Advanced Science 2025, Ternary Attention Framework
- **방식:** Transformer 기반 3원 복합체 상호작용 모델링 (원자-분자-물성 계층 표현)
- **데이터:** `data/PROTAC-fine` 폴더에 포함 (PROTAC-DB 2.0 기반 정제 데이터)
- **특징:** Colab 데모 제공 (`demo.ipynb`)
- **실행 환경:** `protac` 환경에서 시도 가능 (requirements 확인 필요)
- **난이도:** 중~고

### 5-5. DegradoMap (`~/DegradoMap/`)
- **출처:** `github.com/bryanc5864/DegradoMap`
- **논문:** ACM BCB 2026 accepted
- **방식:** AlphaFold 구조 + E3 ligase 정체성만으로 분해 가능성 예측 (PROTAC 구조 불필요)
- **특징:** PROTAC 설계 전 단계에서 타겟 druggability 스크리닝 용도
- **성능:** Target-unseen AUROC 0.603 (6-seed mean), E3-unseen 0.811
- **난이도:** 낮음 (시각화 중심, 빠른 실행)

### 5-6. PROTAC-TS (`~/PROTAC-TS/`)
- **출처:** `github.com/ycu-iil/PROTAC-TS`
- **논문:** JACS Au 2026
- **방식:** ML 기반 세포막 투과성 예측 + 링커 설계 (ChemTSv2 기반)
- **주의:** `tabpfn`, `medchem`, `chemtsv2` 별도 설치 필요 (Python 3.11 환경)
- **투과성 예측만 사용 시:** `make_feature.py` + `make_model.py` → `protac` 환경에서 부분 실행 가능

---

## 6. 각 Phase별 실행 계획 및 명령어

### Phase 0: 데이터 전처리 (`protac` 환경)

```python
# AR PROTAC 추출
import pandas as pd
df = pd.read_excel('protac.xlsx')
ar_df = df[df['Uniprot'] == 'P10275']  # 777개

# DOI → 출판 연도 조회
from habanero import Crossref
cr = Crossref()
def get_year(doi):
    try:
        result = cr.works(ids=doi)
        return result['message']['published']['date-parts'][0][0]
    except:
        return None

# jm4c02226 CSV 읽기 (인코딩 주의)
case_df = pd.read_csv('jm4c02226_si_002.csv', encoding='cp949', header=[0,1])

# SMILES 유효성 검증
from rdkit import Chem
def valid_smiles(s):
    return Chem.MolFromSmiles(s) is not None
```

### Phase 1: Chemical Space 탐색 (`protac` 환경)

```python
from rdkit.Chem import Descriptors, rdMolDescriptors
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from rdkit.Chem import AllChem

# Morgan fingerprint 계산
def get_fp(mol, radius=2, nbits=2048):
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius, nbits)

# 분자 기술자 계산
def calc_descriptors(mol):
    return {
        'MW': Descriptors.MolWt(mol),
        'TPSA': rdMolDescriptors.CalcTPSA(mol),
        'cLogP': Descriptors.MolLogP(mol),
        'HBD': rdMolDescriptors.CalcNumHBD(mol),
        'HBA': rdMolDescriptors.CalcNumHBA(mol),
        'RotBonds': rdMolDescriptors.CalcNumRotatableBonds(mol),
    }
```

### Phase 2: AI 분해 예측 (`protac` 환경)

**DegradeMaster case study 실행:**
```bash
cd ~/DegradeMaster
~/anaconda3/envs/protac/bin/python case_study.py
```

**주의사항:**
- PROTAC-8K 전체 데이터는 Zenodo에서 별도 다운로드 필요
- 현재 체크포인트(1000, 2000 epoch)는 존재하므로 예측만 가능

### Phase 3: 피부 투과성 예측 (`protac` 환경)

#### Potts-Guy 모델 (직접 구현)
```python
from rdkit.Chem import Descriptors
import numpy as np

def potts_guy(mol):
    logP = Descriptors.MolLogP(mol)
    MW = Descriptors.MolWt(mol)
    logKp = 0.71 * logP - 0.0061 * MW - 6.3  # cm/h
    Kp = 10 ** logKp
    return logKp, Kp
```

#### Mitragotri 결정론적 모델 (직접 구현)
```python
def mitragotri(mol):
    # 분자 반지름 추정: r ≈ (3MW/4πρNA)^(1/3), ρ≈1 g/cm³
    MW = Descriptors.MolWt(mol)
    logP = Descriptors.MolLogP(mol)
    r = (3 * MW / (4 * 3.14159 * 6.022e23)) ** (1/3) * 1e7  # Angstrom
    # Mitragotri: log(Kp) = log(D_lip) - log(D_aq) - δ_lip²/RT * (logP의 함수)
    # 간략화: log(Kp) ≈ logP - f(MW)
    logKp = logP - 0.0082 * MW - 4.84
    return logKp
```

### Phase 4: Boltz 공동접힘 (`protac_boltz` 환경)

```bash
conda activate protac_boltz
boltz predict input.yaml --out_dir output/
```

**입력 형식 (YAML):**
```yaml
sequences:
  - protein:
      id: AR
      sequence: "MEVQQGLPYGPGAQHPYQPQHPYQPQHPYPQPQHPYPQPQHPYPQPQHQ..."
  - protein:
      id: VHL
      sequence: "..."
  - ligand:
      id: PROTAC
      smiles: "C1CC(=O)N..."
```

**대상 화합물:** C6 (최적 PK), B3 (DC50=8.97nM, Dmax=76%)

---

## 7. 핵심 위험 요소 및 대응 전략

| 위험 | 현황 | 대응 |
|---|---|---|
| DegradeMaster PROTAC-8K 데이터 미보유 | `processed/` 비어 있음 | case_study만 실행, time-split은 자체 구현 |
| DeepPROTACs mol2 파일 필요 | 포켓 추출 도구(PyMOL) 필요 | 웹서버 활용 (bailab.siais.shanghaitech.edu.cn) |
| PROTAC-TS 별도 환경 필요 | tabpfn, chemtsv2 미설치 | Potts-Guy + Mitragotri로 대체 |
| jm4c02226 CSV 인코딩 | cp949 (UTF-8 아님) | `encoding='cp949'` 명시 필요 |
| AR PROTAC 수량 | protac.xlsx에 777개 존재 | 충분, time split 가능 |

---

## 8. Time Split 전략

```
Train set: ProtacDB_MTL_CLS.csv 내 AR PROTAC (구버전, ~2022년 이전)
Test set:  protac.xlsx 中 ProtacDB_MTL_CLS에 없는 신규 AR PROTAC
           (InChIKey 매칭으로 구별, DOI 연도로 정렬)
검증셋:    jm4c02226_si_002.csv 30개 (실측 Skin Retention Rate 보유)
```

**InChIKey 기반 매칭 코드:**
```python
from rdkit.Chem import MolToInchiKey, MolFromSmiles

def smiles_to_inchikey(smi):
    mol = MolFromSmiles(smi)
    return MolToInchiKey(mol) if mol else None

old_keys = set(df_old['Protac_SMILES'].apply(smiles_to_inchikey).dropna())
new_only = df_new[~df_new['Smiles'].apply(smiles_to_inchikey).isin(old_keys)]
```

---

## 9. 보고서 실행 순서 요약 (Day 1 → Day 2)

### Day 1 (오늘, 2026-06-17)
1. `~/ar_protac_project/` 에서 Phase 0 스크립트 실행 (AR 필터링, time split, DOI 연도 조회)
2. Phase 1 Chemical space 탐색 (PCA/t-SNE/UMAP, 기술자 계산, bRo5 overlay)
3. DegradeMaster `case_study.py` 테스트 실행

### Day 2 (2026-06-18)
1. Phase 2 AI 모델 time-split 평가 (DC50 예측 RMSE, AUROC)
2. Phase 3 Potts-Guy + Mitragotri 구현 및 30개 case study 검증
3. Phase 4 Boltz (시간 여유 시, C6 화합물)
4. Phase 5 신규 분자 SMILES 제안 (3~5개)
5. 보고서 작성

---

## 10. 참고 자원

| 도구 | 접근 방법 |
|---|---|
| DeepPROTACs 웹서버 | bailab.siais.shanghaitech.edu.cn/services/deepprotacs/ |
| SwissADME | swissadme.ch |
| ADMETlab 3.0 | admetlab3.scbdd.com |
| PROTAC-DB 3.0 | cadd.zju.edu.cn/protacdb/ |
| DegradeMaster 데이터 | zenodo.org/records/14728925 |
