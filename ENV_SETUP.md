# 모델 환경 설정 가이드

> **작성일:** 2026-06-17  
> **전략:** 디스크 최소화 — 기존 `protac` env에 최대한 통합, 부득이한 경우만 별도 env 생성

---

## 환경 구성 요약

| conda 환경 | Python | 용도 | 모델 |
|---|---|---|---|
| `protac` (기존) | 3.10 | 메인 분석 환경 | DeepPROTACs, DegradeMaster, AiPROTAC, PROTAC-STAN, DegradoMap |
| `protac_ts` (신규) | 3.11 | PROTAC-TS 전용 | PROTAC-TS |
| `protac_boltz` (기존) | 3.10 | Boltz 공동접힘 | Boltz |

> chemtsv2는 Python ≥3.11 필수 → `protac_ts` 별도 생성 불가피  
> 나머지 5개 모델은 모두 `protac` env에서 실행 가능

---

## 1. `protac` 환경 — 5개 모델 통합

### Python 경로
```
~/anaconda3/envs/protac/bin/python
```

### 핵심 패키지 현황

| 패키지 | 버전 | 설치 방법 |
|---|---|---|
| rdkit | 2026.3.1 | 기존 |
| torch | 2.8.0+cu128 | 기존 |
| torch-geometric | 2.7.0 | 기존 |
| torch-scatter | 2.1.2 | 기존 |
| torch-sparse | 0.6.18 | 기존 |
| scikit-learn | 1.7.2 | 기존 |
| pandas | 2.3.3 | 기존 |
| numpy | 2.2.6 | 기존 |
| networkx | 3.4.2 | 기존 |
| umap-learn | 0.5.11 | 기존 |
| wandb | 0.27.2 | 기존 |
| toml | 0.10.2 | 기존 |
| e3nn | 0.6.0 | 기존 |
| fair-esm | 2.0.0 | 기존 |
| habanero | 2.4.0 | `pip install habanero` |
| tabpfn | 2.0.3 | `pip install tabpfn==2.0.3` |
| medchem | 2.0.5 | `pip install medchem==2.0.5` |
| dgl | 1.1.3 (CPU) | `pip install dgl==1.1.3` |

### DGL 호환성 패치

DGL 1.1.3 설치 후 torchdata 최신 버전과의 호환성을 위해 stub 파일 필요:

```bash
STUB_DIR=~/anaconda3/envs/protac/lib/python3.10/site-packages/torchdata
mkdir -p $STUB_DIR/datapipes/iter

cat > $STUB_DIR/datapipes/__init__.py << 'EOF'
# Compatibility stub for DGL
EOF

cat > $STUB_DIR/datapipes/iter/__init__.py << 'EOF'
# Compatibility stub for DGL graphbolt
class IterDataPipe:
    pass
EOF
```

> **주의:** DGL 1.1.3은 CPU 전용 빌드 — 그래프 구조 연산은 CPU, 텐서 연산은 torch CUDA 활용  
> AiPROTAC 실행 시 `--device cpu` 옵션 추가 필요 (또는 GPU 연산은 torch tensor로만 처리)

### 환경 재현 명령어

```bash
# 신규 설치 시 (기존 protac env가 없는 경우)
conda create -n protac python=3.10 -y
conda activate protac

# PyTorch + PyG (CUDA 12.8 버전)
pip install torch==2.8.0 torchvision torchaudio --index-url https://download.pytorch.org/whl/cu128
pip install torch-geometric torch-scatter torch-sparse torch-cluster torch-spline-conv \
    -f https://data.pyg.org/whl/torch-2.8.0+cu128.html

# 핵심 패키지
pip install rdkit pandas numpy scipy scikit-learn matplotlib seaborn tqdm
pip install openpyxl habanero python-docx umap-learn networkx

# 모델별 추가 패키지
pip install wandb toml e3nn fair-esm
pip install tabpfn==2.0.3 medchem==2.0.5
pip install dgl==1.1.3
```

---

## 2. `protac_ts` 환경 — PROTAC-TS 전용

chemtsv2는 Python 3.11~3.12 전용 패키지이므로 별도 환경 필요.

### Python 경로
```
~/anaconda3/envs/protac_ts/bin/python
```

### 설치 패키지

| 패키지 | 버전 |
|---|---|
| Python | 3.11 |
| chemtsv2 | 1.1.2 |
| tabpfn | 2.0.3 |
| medchem | 2.0.5 |
| scikit-learn | 1.5.1 |
| rdkit | 2023.9.6 |
| torch | 2.12.0 (CPU) |

### 환경 재현 명령어

```bash
conda create -n protac_ts python=3.11 --no-default-packages -y
~/anaconda3/envs/protac_ts/bin/pip install \
    chemtsv2==1.1.2 tabpfn==2.0.3 medchem==2.0.5 \
    scikit-learn==1.5.1 rdkit
```

---

## 3. `protac_boltz` 환경 — Boltz 공동접힘 예측 (기존)

### Python 경로
```
~/anaconda3/envs/protac_boltz/bin/python
```

| 패키지 | 버전 |
|---|---|
| boltz | 2.2.1 |
| torch | 2.8.0+cu128 |
| rdkit | 2025.9.5 |
| pytorch-lightning | 2.5.0 |

---

## 4. 모델별 실행 환경 및 import 검증

### ✅ DeepPROTACs
```bash
cd ~/DeepPROTACs
~/anaconda3/envs/protac/bin/python -c "
from model import GraphConv, SmilesNet, ProtacModel
from prepare_data import GraphData
print('OK')
"
# 단일 예측 실행
~/anaconda3/envs/protac/bin/python single_prediction.py single_test
```
> **입력:** `ligase_ligand.mol2`, `ligase_pocket.mol2`, `target_ligand.mol2`, `target_pocket.mol2`, `linker.smi`  
> mol2 파일 없을 시 → 웹서버: bailab.siais.shanghaitech.edu.cn/services/deepprotacs/

---

### ✅ DegradeMaster
```bash
cd ~/DegradeMaster
~/anaconda3/envs/protac/bin/python -c "
from model import GraphConv, ProtacModel, EGNNConv
from protacloader import PROTACSet, collater
print('OK')
"
# Case study 실행 (체크포인트 존재: checkpoint/1000/, 2000/)
~/anaconda3/envs/protac/bin/python case_study.py
# 전체 학습/평가 (PROTAC-8K 데이터 필요: zenodo.org/records/14728925)
~/anaconda3/envs/protac/bin/python main.py --config config/config.yml
```
> `process.py`의 openbabel(pybel)은 전처리 전용 — 추론 시 불필요

---

### ✅ AiPROTAC
```bash
cd ~/AiPROTAC
~/anaconda3/envs/protac/bin/python -c "
import dgl
from dgl.dataloading import GraphDataLoader
from model.AiPROTAC import GraphBasedModel
print('OK')
"
# 학습 실행 (CPU 모드)
~/anaconda3/envs/protac/bin/python train_AiPROTAC.py --device cpu
# GPU 모드 (DGL 1.1.3은 CPU 전용 — torch tensor는 GPU 사용 가능)
# DGL graph은 CPU에 두고 model만 GPU로 이동하면 동작
```
> DGL 1.1.3 CPU: 그래프 구조 CPU, torch 텐서 GPU — 학습 가능하나 그래프 조작 속도 제한적

---

### ✅ PROTAC-STAN
```bash
cd ~/PROTAC-STAN
~/anaconda3/envs/protac/bin/python -c "
from data import PROTACData
from data_loader import PROTACLoader
from model import PROTAC_STAN
print('OK')
"
# 데모 실행 (data/demo/ 폴더의 ESM 임베딩 precomputed)
~/anaconda3/envs/protac/bin/python inference.py --config config_demo.toml
# 학습/평가 (data/PROTAC-fine/ 데이터 포함됨)
~/anaconda3/envs/protac/bin/python main.py --config config.toml
```
> ESM 임베딩은 `data/*/esm_s_map.pkl`에 미리 계산되어 있음 — ESM 모델 다운로드 불필요

---

### ✅ DegradoMap
```bash
cd ~/DegradoMap
~/anaconda3/envs/protac/bin/python -c "
from src.models.degradomap import DegradoMap
import e3nn, esm, wandb
print('OK')
"
# 학습 실행 (data acquisition 스크립트로 데이터 구축 후)
~/anaconda3/envs/protac/bin/python scripts/train.py --phase 2 --splits target_unseen
```
> 데이터 수집 스크립트(`src/data/`)는 selenium 필요 (데이터 이미 있으면 불필요)  
> 사전 학습된 결과: `results/*.json` 참조 가능

---

### ✅ PROTAC-TS
```bash
cd ~/PROTAC-TS
~/anaconda3/envs/protac_ts/bin/python -c "
import chemtsv2, tabpfn, medchem
print('OK')
"
# 피처 생성
~/anaconda3/envs/protac_ts/bin/python make_feature.py -c config/setting_feature.yaml
# 모델 학습
~/anaconda3/envs/protac_ts/bin/python make_model.py -c config/setting_model.yaml
# 링커 설계
~/anaconda3/envs/protac_ts/bin/chemtsv2 -c config/setting_protacts.yaml
```
> 데이터: PROTAC-DB 3.0에서 PROTACs + linkers CSV 다운로드 필요 (cadd.zju.edu.cn/protacdb/downloads)

---

## 5. 공간 최적화 전략 요약

| 방법 | 절약 효과 |
|---|---|
| 5개 모델을 `protac` env 1개에 통합 | torch+cuda 중복 설치 4회 방지 (~15GB 절약) |
| DGL 1.1.3 (CPU, pip) 사용 | CUDA DGL 별도 env 불필요 (~3GB 절약) |
| PROTAC-STAN ESM precomputed 활용 | ESM 모델 가중치 다운로드 불필요 (~1.6GB 절약) |
| `protac_ts` 최소 설치 (--no-default-packages) | conda base 패키지 복사 방지 |
| conda 패키지 캐시 공유 | 동일 패키지 재다운로드 방지 |

---

## 6. 트러블슈팅

### DGL `torchdata.datapipes` 오류
```
ModuleNotFoundError: No module named 'torchdata.datapipes'
```
**원인:** torchdata 0.7+ 이후 datapipes 모듈 제거  
**해결:** DGL 1.1.3 사용 (graphbolt 없음, datapipes 불필요)  
만약 DGL 2.x 사용 필요 시 stub 생성:
```bash
STUB=$(python -c "import torchdata,os; print(os.path.dirname(torchdata.__file__))")
mkdir -p $STUB/datapipes/iter
echo "class IterDataPipe: pass" > $STUB/datapipes/iter/__init__.py
echo "" > $STUB/datapipes/__init__.py
```

### AiPROTAC DGL CUDA 오류
```
DGLError: Device API cuda is not enabled
```
**원인:** pip DGL 1.1.3은 CPU 전용 빌드  
**해결:** `--device cpu` 옵션으로 실행, 또는 model 파라미터만 GPU (`model.cuda()`)

### chemtsv2 Python 버전 오류
```
ERROR: No matching distribution found for chemtsv2==1.1.2
```
**원인:** chemtsv2 1.1.x는 Python 3.11~3.12 전용  
**해결:** `protac_ts` env (Python 3.11) 사용
