"""
PROTAC-STAN 입력 custom.csv (29행) 작성 — C18/C19 반영.

- 입력: data/processed/case_study_clean.csv (29 SMILES)
- 컬럼: data.py의 `columns` 리스트와 1자 단위 일치 + Uniprot/E3 ligase Uniprot
- 전 29행: Uniprot=P10275(AR), E3 ligase Uniprot=Q96SW2(CRBN)  ← C19
  (글루타리미드 29/29 → 전 계열 CRBN. VHL 라벨 금지)
- 9개 물성은 RDKit 재계산: XLogP3←Crippen MolLogP, TPSA←RDKit TPSA 등
- 출력: ~/PROTAC-STAN/data/custom/custom.csv (29행) + 백업본
- 재실행 전 processed/custom 캐시 삭제
"""
import os
import shutil

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors

RDLogger.DisableLog("rdApp.*")

ROOT = "/home/kimjisan95/ar_protac_project"
STAN_CUSTOM = os.path.expanduser("~/PROTAC-STAN/data/custom")
AR_UNIPROT = "P10275"
CRBN_UNIPROT = "Q96SW2"

# data.py columns 와 정확히 일치해야 함
COLS = [
    "Molecular Weight", "Exact Mass", "XLogP3", "Heavy Atom Count",
    "Ring Count", "Hydrogen Bond Acceptor Count", "Hydrogen Bond Donor Count",
    "Rotatable Bond Count", "Topological Polar Surface Area",
]


def props(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return {
        "Molecular Weight": Descriptors.MolWt(mol),
        "Exact Mass": Descriptors.ExactMolWt(mol),
        "XLogP3": Crippen.MolLogP(mol),          # XLogP3 대체: Crippen logP
        "Heavy Atom Count": mol.GetNumHeavyAtoms(),
        "Ring Count": rdMolDescriptors.CalcNumRings(mol),
        "Hydrogen Bond Acceptor Count": rdMolDescriptors.CalcNumHBA(mol),
        "Hydrogen Bond Donor Count": rdMolDescriptors.CalcNumHBD(mol),
        "Rotatable Bond Count": rdMolDescriptors.CalcNumRotatableBonds(mol),
        "Topological Polar Surface Area": Descriptors.TPSA(mol),
    }


cs = pd.read_csv(os.path.join(ROOT, "data", "processed", "case_study_clean.csv"))
rows = []
for _, r in cs.iterrows():
    p = props(r["Smiles"])
    rows.append({
        "ID": r["ID"],
        "Uniprot": AR_UNIPROT,
        "E3 ligase Uniprot": CRBN_UNIPROT,
        "Smiles": r["Smiles"],
        **p,
    })
out = pd.DataFrame(rows)
# data.py가 요구하는 컬럼 순서(ID는 추가 컬럼, 무시됨)
ordered = ["ID", "Uniprot", "E3 ligase Uniprot", "Smiles"] + COLS
out = out[ordered]

# 검증: VHL cue (하이드록시프롤린) / 글루타리미드 매치
glut = Chem.MolFromSmarts("O=C1NCCCC1=O")  # piperidine-2,6-dione 근사
n_glut = sum(Chem.MolFromSmiles(s).HasSubstructMatch(Chem.MolFromSmarts("O=C1CCCC(=O)N1"))
             for s in out["Smiles"])
print(f"글루타리미드(2,6-dione) 매치: {n_glut}/{len(out)} (C19: 전 계열 CRBN 근거)")

dst = os.path.join(STAN_CUSTOM, "custom.csv")
# 기존 데모 백업
if os.path.exists(dst) and not os.path.exists(dst + ".demo_bak"):
    shutil.copy(dst, dst + ".demo_bak")
out.to_csv(dst, index=False)
# 우리 프로젝트에도 사본
out.to_csv(os.path.join(ROOT, "data", "processed", "stan_custom_input.csv"), index=False)

# 캐시 삭제 (재처리 강제)
cache = os.path.join(STAN_CUSTOM, "processed", "custom")
if os.path.exists(cache):
    shutil.rmtree(cache)
    print(f"캐시 삭제: {cache}")

print(f"작성 완료: {dst} ({len(out)}행)")
print(out[["ID", "Uniprot", "E3 ligase Uniprot", "Molecular Weight",
           "XLogP3", "Topological Polar Surface Area"]].to_string(index=False))
