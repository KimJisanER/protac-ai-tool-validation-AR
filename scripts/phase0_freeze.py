"""
Phase 0 — 데이터 동결 스크립트 (§0-END).

모든 카운트의 단일 출처. 보고서/표/그림은 이 스크립트 출력만 인용(하드코딩 금지).

검증 대상 카운트:
  - 구DB(P10275) 행수/unique ik14, 신DB(P10275) 행수/unique ik14
  - Degrader_Class_mask==1 = 507행(양성/음성), AR(P10275) 관측 양성 절대수
  - 구조분할: new-only AR ik14 (= 298 기대)
  - 누수: case study ik14 ∩ STAN train ik14 (= B3, B5 기대)
  - case study 29행(A16/B6/C7), retention 29, DC50/Dmax 13(B/C), 수치형 DC50
  - jm4c AND-양성(DC50<100 AND Dmax>=80) 개수, DC50<100 단독 양성 개수

산출물:
  data/processed/ar_old.csv, ar_new.csv, case_study_clean.csv, leakage_report.csv
  outputs/table1.csv
"""
import os
import re
import sys

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

ROOT = "/home/kimjisan95/ar_protac_project"
STAN = os.path.expanduser("~/PROTAC-STAN")
PROC = os.path.join(ROOT, "data", "processed")
OUT = os.path.join(ROOT, "outputs")
os.makedirs(PROC, exist_ok=True)
os.makedirs(OUT, exist_ok=True)

AR_UNIPROT = "P10275"


def ik14(smiles):
    """RDKit MolToInchiKey 앞 14자(골격키). 실패 시 None."""
    if not isinstance(smiles, str) or not smiles.strip():
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        return Chem.MolToInchiKey(mol)[:14]
    except Exception:
        return None


def parse_numeric(val):
    """DC50/Dmax 문자열 → (point_value or None, censored_flag, raw).
    '÷'(cp949 mojibake)→'±'. '8.45 ± 0.16' → 8.45. '>1000' → censored.
    슬래시 다중값은 첫 값. N.D./공백 → None."""
    if val is None or (isinstance(val, float) and np.isnan(val)):
        return None, None, val
    s = str(val).strip().replace("÷", "±")
    if s == "" or s.upper() in ("N.D.", "ND", "NA", "NAN", "-"):
        return None, None, val
    censored = None
    m = re.match(r"^\s*([<>]=?)", s)
    if m:
        censored = m.group(1)
    # 첫 부동소수 추출
    fm = re.search(r"[-+]?\d*\.?\d+", s)
    pt = float(fm.group()) if fm else None
    return pt, censored, val


# ---------------------------------------------------------------------------
log = []
def P(*a):
    line = " ".join(str(x) for x in a)
    print(line)
    log.append(line)

P("=" * 70)
P("Phase 0 데이터 동결 — 모든 수치의 단일 출처")
P("=" * 70)

# --- 1) 구DB ---------------------------------------------------------------
old = pd.read_csv(os.path.join(ROOT, "ProtacDB_MTL_CLS.csv"))
P(f"\n[구DB] 전체 {old.shape[0]}행 x {old.shape[1]}열")

ar_old = old[old["Uniprot"] == AR_UNIPROT].copy()
ar_old["ik14"] = ar_old["Protac_SMILES"].apply(ik14)
P(f"[구DB] AR(Uniprot=={AR_UNIPROT}) = {len(ar_old)}행, "
  f"unique ik14 = {ar_old['ik14'].nunique()}")

# Degrader_Class mask 적용 (전체)
obs = old[old["Degrader_Class_mask"] == 1]
n_pos = int((obs["Degrader_Class"] == 1).sum())
n_neg = int((obs["Degrader_Class"] == 0).sum())
P(f"[구DB] Degrader_Class_mask==1 관측 라벨 = {len(obs)}행 "
  f"(양성 {n_pos} / 음성 {n_neg})  ← C14: 약 50:50 균형")
miss_rate = 1 - len(obs) / len(old)
P(f"[구DB] 라벨 결손율(미관측) = {(len(old)-len(obs))}/{len(old)} = {miss_rate:.1%}  ← C14 정정")

# AR 관측 양성 절대수 (C4)
ar_obs = ar_old[ar_old["Degrader_Class_mask"] == 1]
ar_pos = int((ar_obs["Degrader_Class"] == 1).sum())
ar_neg = int((ar_obs["Degrader_Class"] == 0).sum())
P(f"[구DB] AR 관측 라벨 = {len(ar_obs)}행 (양성 {ar_pos} / 음성 {ar_neg})  "
  f"← C4: 재학습 포기 사유(양성 절대수 {ar_pos})")

# --- 2) 신DB ---------------------------------------------------------------
new = pd.read_excel(os.path.join(ROOT, "protac.xlsx"))
P(f"\n[신DB] 전체 {new.shape[0]}행 x {new.shape[1]}열")
ar_new = new[new["Uniprot"] == AR_UNIPROT].copy()
ar_new["ik14"] = ar_new["Smiles"].apply(ik14)
P(f"[신DB] AR(Uniprot=={AR_UNIPROT}) = {len(ar_new)}행, "
  f"unique ik14 = {ar_new['ik14'].nunique()}")

# DC50/Dmax 결손율 (신DB AR)
dc50_miss = ar_new["DC50 (nM)"].isna().mean()
dmax_miss = ar_new["Dmax (%)"].isna().mean()
P(f"[신DB] AR DC50 결손율 = {dc50_miss:.1%}, Dmax 결손율 = {dmax_miss:.1%}")

# --- 3) 구조분할 (C13) -----------------------------------------------------
old_iks = set(ar_old["ik14"].dropna())
new_iks = set(ar_new["ik14"].dropna())
new_only = new_iks - old_iks
overlap = new_iks & old_iks
P(f"\n[구조분할] 구DB AR unique={len(old_iks)}, 신DB AR unique={len(new_iks)}, "
  f"overlap={len(overlap)}")
P(f"[구조분할] new-only AR ik14 (test 후보) = {len(new_only)}  ← C13 기대값 298")

# --- 4) case study --------------------------------------------------------
cs = pd.read_csv(os.path.join(ROOT, "jm4c02226_si_002.csv"),
                 encoding="cp949", header=[0, 1])
cols = list(cs.columns)
ids = cs.iloc[:, 0].astype(str).tolist()
smiles = cs.iloc[:, 1].tolist()
retention_raw = cs.iloc[:, 2].tolist()
dc50_raw = cs.iloc[:, 3].tolist()
dmax_raw = cs.iloc[:, 4].tolist()

series = [i[0] for i in ids]  # A/B/C
rows = []
for i in range(len(ids)):
    ret_pt, _, _ = parse_numeric(retention_raw[i])
    dc_pt, dc_cens, _ = parse_numeric(dc50_raw[i])
    dm_pt, _, _ = parse_numeric(dmax_raw[i])
    rows.append({
        "ID": ids[i],
        "series": series[i],
        "Smiles": smiles[i],
        "ik14": ik14(smiles[i]),
        "skin_retention_pct": ret_pt,
        "dc50_nM": dc_pt,
        "dc50_censored": dc_cens,
        "dmax_pct": dm_pt,
        "dc50_raw": dc50_raw[i],
        "dmax_raw": dmax_raw[i],
    })
cs_clean = pd.DataFrame(rows)

n_A = (cs_clean["series"] == "A").sum()
n_B = (cs_clean["series"] == "B").sum()
n_C = (cs_clean["series"] == "C").sum()
P(f"\n[case study] {len(cs_clean)}행 = A{n_A}/B{n_B}/C{n_C}")
P(f"[case study] SMILES 파싱 성공 ik14 = {cs_clean['ik14'].notna().sum()}/{len(cs_clean)}")
P(f"[case study] retention 보유 = {cs_clean['skin_retention_pct'].notna().sum()}")
bc = cs_clean[cs_clean["series"].isin(["B", "C"])]
P(f"[case study] DC50 보유(B/C) = {cs_clean['dc50_nM'].notna().sum()} "
  f"(B/C 행 {len(bc)})")
P(f"[case study] Dmax 보유 = {cs_clean['dmax_pct'].notna().sum()}")
# 수치형 DC50 (비-censored)
num_dc50 = cs_clean[(cs_clean["dc50_nM"].notna()) & (cs_clean["dc50_censored"].isna())]
P(f"[case study] 수치형(비검열) DC50 = {len(num_dc50)}개: "
  f"{dict(zip(num_dc50['ID'], num_dc50['dc50_nM']))}")

# jm4c 양성 정의 (C20)
# STAN AND: DC50<100 AND Dmax>=80
bc_valid = cs_clean[(cs_clean["dc50_nM"].notna()) & (cs_clean["dmax_pct"].notna())]
stan_pos = bc_valid[(bc_valid["dc50_nM"] < 100) & (bc_valid["dmax_pct"] >= 80)]
dc_only_pos = cs_clean[(cs_clean["dc50_nM"].notna()) & (cs_clean["dc50_nM"] < 100)
                       & (cs_clean["dc50_censored"].isna())]
P(f"[case study] STAN AND 양성(DC50<100 AND Dmax>=80) = {len(stan_pos)}개 "
  f"(max Dmax={bc_valid['dmax_pct'].max() if len(bc_valid) else 'NA'})  ← C20 기대 0")
P(f"[case study] DC50<100 단독 양성 = {len(dc_only_pos)}개: "
  f"{dc_only_pos['ID'].tolist()}  ← C20 기대 B3,C5")

# --- 5) 누수 (C5) ----------------------------------------------------------
train = pd.read_csv(os.path.join(STAN, "data", "PROTAC-fine", "train_compound_smiles.csv"))
train_iks = set(train["SMILES"].apply(ik14).dropna())
P(f"\n[누수] STAN train_compound_smiles unique ik14 = {len(train_iks)}")
cs_clean["leaked_stan_train"] = cs_clean["ik14"].isin(train_iks)
cs_clean["leaked_old_db"] = cs_clean["ik14"].isin(old_iks)
leaked = cs_clean[cs_clean["leaked_stan_train"] | cs_clean["leaked_old_db"]]
P(f"[누수] case study ∩ (STAN train ∪ 구DB) = {len(leaked)}개: "
  f"{leaked['ID'].tolist()}  ← C5 기대 B3,B5")
for _, r in leaked.iterrows():
    P(f"        {r['ID']}: ik14={r['ik14']} "
      f"(STAN_train={r['leaked_stan_train']}, old_db={r['leaked_old_db']})")

# --- 산출물 저장 -----------------------------------------------------------
ar_old.to_csv(os.path.join(PROC, "ar_old.csv"), index=False)
ar_new.to_csv(os.path.join(PROC, "ar_new.csv"), index=False)
cs_clean.to_csv(os.path.join(PROC, "case_study_clean.csv"), index=False)
leak_rep = cs_clean[["ID", "series", "ik14", "leaked_stan_train", "leaked_old_db"]]
leak_rep.to_csv(os.path.join(PROC, "leakage_report.csv"), index=False)

# --- 표1 -------------------------------------------------------------------
table1 = pd.DataFrame([
    {"dataset": "구DB AR (P10275)", "n_rows": len(ar_old),
     "unique_ik14": ar_old["ik14"].nunique(),
     "dc50_miss_pct": "-", "dmax_miss_pct": "-",
     "note": f"관측라벨 {len(ar_obs)} (양성{ar_pos}/음성{ar_neg})"},
    {"dataset": "신DB AR (P10275)", "n_rows": len(ar_new),
     "unique_ik14": ar_new["ik14"].nunique(),
     "dc50_miss_pct": f"{dc50_miss:.1%}", "dmax_miss_pct": f"{dmax_miss:.1%}",
     "note": f"new-only ik14={len(new_only)} (구조분할 test 후보)"},
    {"dataset": "case study (jm4c)", "n_rows": len(cs_clean),
     "unique_ik14": cs_clean["ik14"].nunique(),
     "dc50_miss_pct": f"{cs_clean['dc50_nM'].isna().mean():.1%}",
     "dmax_miss_pct": f"{cs_clean['dmax_pct'].isna().mean():.1%}",
     "note": f"A{n_A}/B{n_B}/C{n_C}; 누수 {leaked['ID'].tolist()}"},
])
table1.to_csv(os.path.join(OUT, "table1.csv"), index=False)
P("\n[저장] data/processed/{ar_old,ar_new,case_study_clean,leakage_report}.csv")
P("[저장] outputs/table1.csv")

with open(os.path.join(OUT, "phase0_log.txt"), "w") as f:
    f.write("\n".join(log))
P("\n동결 완료.")
