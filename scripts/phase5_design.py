"""
Phase 5 — 신규 분자 설계 (hypothesis generation, 순환논증 경고).

설계 근거 = 모델 점수가 아니라 Phase 1 SAR/물성 규칙:
  - case study 29/29가 CRBN/글루타리미드 + darolutamide류 AR 워헤드 (C19)
  - leakage-free 활성 앵커 C5(70.85nM)/C1(103.87nM)/C6(199.5nM)의 약리단 보존
  - 말단 위치만 보수적 변형(대사 soft-spot 차단 F, 친유성 미세조정 Me) → SAR 가설
모델(STAN) 출력은 "약한 사전확률(weak prior)"로만, AD 거리 동반, 실험검증 필요 명시.

산출물: outputs/table4.csv, data/processed/design_custom.csv (STAN 입력)
"""
import os
import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors, AllChem, DataStructs

RDLogger.DisableLog("rdApp.*")
ROOT = "/home/kimjisan95/ar_protac_project"
OUT = os.path.join(ROOT, "outputs")
PROC = os.path.join(ROOT, "data", "processed")
STAN = os.path.expanduser("~/PROTAC-STAN")

# leakage-free 활성 앵커 원본
C1 = "O=C(NC1C(NC(CC1)=O)=O)C2=CC=C(N3CCN(C4CCN(C5=CC=C(C(N[C@H]6CC[C@H](OC7=CC(Cl)=C(C#N)C=C7)CC6)=O)C=C5)CC4)CC3)C=C2"
C5 = "O=C(C1=CC=C(C=C1)C2CCN(CC2)CC(C3)CN3C4=CC=C(C(NC5CCC(NC5=O)=O)=O)C=C4)N[C@@H](CC6)CC[C@H]6OC(C=C7)=CC(Cl)=C7C#N"
C6 = "O=C(C1=CC=C(C=C1)C2CCN(CC2)CC3CCCN(C4=CC=C(C(NC5CCC(NC5=O)=O)=O)C=C4)C3)N[C@@H](CC6)CC[C@H]6OC(C=C7)=CC(Cl)=C7C#N"

# 신규 3종 — 약리단(글루타리미드+워헤드) 보존, 말단만 변형
designs = [
    {"ID": "D1", "parent": "C1",
     # 글루타리미드 연결 benzamide 고리에 ortho-F 도입(대사 안정화), 약리단 불변
     "Smiles": "O=C(NC1C(NC(CC1)=O)=O)C2=CC=C(N3CCN(C4CCN(C5=CC=C(C(N[C@H]6CC[C@H](OC7=CC(Cl)=C(C#N)C=C7)CC6)=O)C=C5)CC4)CC3)C(F)=C2",
     "rationale": "C1(103.87nM) 약리단 보존; CRBN benzamide ortho-F → 대사 soft-spot 차단·logD 미세상향. 링커/워헤드 불변."},
    {"ID": "D2", "parent": "C5",
     # 글루타리미드-phenyl에 methyl 도입(입체 미세조정), 약리단 불변
     "Smiles": "O=C(C1=CC=C(C=C1)C2CCN(CC2)CC(C3)CN3C4=CC=C(C(NC5CCC(NC5=O)=O)=O)C(C)=C4)N[C@@H](CC6)CC[C@H]6OC(C=C7)=CC(Cl)=C7C#N",
     "rationale": "C5(70.85nM, 최고활성) 약리단 보존; phenyl-glutarimide ortho-Me → 링커 회전 제약(전엔트로피↓)·삼원복합체 가설. 워헤드 불변."},
    {"ID": "D3", "parent": "C6",
     # AR 워헤드 benzamide 고리에 추가 F(darolutamide류 다중치환 모사)
     "Smiles": "O=C(C1=CC=C(F)C(=C1)C2CCN(CC2)CC3CCCN(C4=CC=C(C(NC5CCC(NC5=O)=O)=O)C=C4)C3)N[C@@H](CC6)CC[C@H]6OC(C=C7)=CC(Cl)=C7C#N",
     "rationale": "C6(199.5nM, leakage-free 리드) 약리단 보존; AR 워헤드 인접 benzoyl 고리 F 도입 → darolutamide류 다중치환 모사, 결합 미세조정 가설."},
]


def descs(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    return {
        "valid": True,
        "MW": Descriptors.MolWt(m), "logP": Crippen.MolLogP(m),
        "TPSA": Descriptors.TPSA(m), "HBD": rdMolDescriptors.CalcNumHBD(m),
        "HBA": rdMolDescriptors.CalcNumHBA(m),
        "RotBond": rdMolDescriptors.CalcNumRotatableBonds(m),
        "glutarimide": m.HasSubstructMatch(Chem.MolFromSmarts("O=C1CCCC(=O)N1")),
        "AR_warhead_CN_Cl": m.HasSubstructMatch(Chem.MolFromSmarts("N#Cc1ccccc1Cl")) or
                            m.HasSubstructMatch(Chem.MolFromSmarts("Clc1ccc(OC)cc1C#N")) or
                            (m.HasSubstructMatch(Chem.MolFromSmarts("C#N")) and
                             m.HasSubstructMatch(Chem.MolFromSmarts("Cl"))),
    }


# STAN 학습셋 FP (AD 거리: Tanimoto 최근접)
train = pd.read_csv(os.path.join(STAN, "data", "PROTAC-fine", "train_compound_smiles.csv"))
train_fps = []
for s in train["SMILES"]:
    m = Chem.MolFromSmiles(s)
    if m:
        train_fps.append(AllChem.GetMorganFingerprintAsBitVect(m, 2, 2048))


def ad_nn(smi):
    m = Chem.MolFromSmiles(smi)
    fp = AllChem.GetMorganFingerprintAsBitVect(m, 2, 2048)
    sims = DataStructs.BulkTanimotoSimilarity(fp, train_fps)
    return max(sims)  # 최근접 유사도(높을수록 AD 내부)


rows = []
for d in designs:
    dd = descs(d["Smiles"])
    assert dd is not None, f"{d['ID']} SMILES 파싱 실패"
    nn = ad_nn(d["Smiles"])
    rows.append({**d, **dd, "AD_nn_Tanimoto": round(nn, 3)})

df = pd.DataFrame(rows)
print("신규 설계 3종 검증:")
print(df[["ID", "parent", "valid", "MW", "logP", "TPSA", "RotBond",
          "glutarimide", "AR_warhead_CN_Cl", "AD_nn_Tanimoto"]].to_string(index=False))

# STAN 입력 csv 작성(weak prior 추론용)
COLS = ["Molecular Weight", "Exact Mass", "XLogP3", "Heavy Atom Count",
        "Ring Count", "Hydrogen Bond Acceptor Count", "Hydrogen Bond Donor Count",
        "Rotatable Bond Count", "Topological Polar Surface Area"]
custom = []
for d in designs:
    m = Chem.MolFromSmiles(d["Smiles"])
    custom.append({
        "ID": d["ID"], "Uniprot": "P10275", "E3 ligase Uniprot": "Q96SW2",
        "Smiles": d["Smiles"],
        "Molecular Weight": Descriptors.MolWt(m), "Exact Mass": Descriptors.ExactMolWt(m),
        "XLogP3": Crippen.MolLogP(m), "Heavy Atom Count": m.GetNumHeavyAtoms(),
        "Ring Count": rdMolDescriptors.CalcNumRings(m),
        "Hydrogen Bond Acceptor Count": rdMolDescriptors.CalcNumHBA(m),
        "Hydrogen Bond Donor Count": rdMolDescriptors.CalcNumHBD(m),
        "Rotatable Bond Count": rdMolDescriptors.CalcNumRotatableBonds(m),
        "Topological Polar Surface Area": Descriptors.TPSA(m),
    })
cdf = pd.DataFrame(custom)[["ID", "Uniprot", "E3 ligase Uniprot", "Smiles"] + COLS]
cdf.to_csv(os.path.join(PROC, "design_custom.csv"), index=False)
df.to_csv(os.path.join(OUT, "table4_predesign.csv"), index=False)
print("\n[저장] data/processed/design_custom.csv (STAN 입력), outputs/table4_predesign.csv")
