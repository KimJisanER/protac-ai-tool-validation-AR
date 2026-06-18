"""
Phase 3 — 피부 투과/잔류 (재설계).

3-1. 전신투과 proxy: Potts-Guy logKp = 0.71*logP - 0.0061*MW - 6.3
     (RDKit MolLogP + MolWt 직접 구현). "낮을수록 국소제 유리"로 명명.
3-2. AD 외삽 정량: MW>750 / TPSA>140 플래그, AD 밖 비율 → [표3]
3-3. 2D 분리 시각화 [그림3]: x=logKp(전신투과), y=실측 retention(피부저류).
     B/C 활성 테두리, 이상 사분면(낮은 Kp·높은 retention) 음영.
3-4. 저비용 2D 카멜레온 [그림3b]: TPSA vs RotBond.

⑤: Potts-Guy AD 외삽 → 정량 신뢰 불가; Skin Retention≠Kp; n작음.
산출물: outputs/table3.csv, fig3_separation.png, fig3b_chameleon.png, phase3_log.txt
"""
import os
import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors
from scipy.stats import spearmanr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import plotstyle  # noqa: F401  (한글 폰트)

RDLogger.DisableLog("rdApp.*")
ROOT = "/home/kimjisan95/ar_protac_project"
OUT = os.path.join(ROOT, "outputs")
PROC = os.path.join(ROOT, "data", "processed")
log = []
def P(*a):
    s = " ".join(str(x) for x in a); print(s); log.append(s)


def potts_guy_logkp(logp, mw):
    """Potts-Guy 1992 (Flynn dataset): logKp(cm/h) = 0.71*logP - 0.0061*MW - 6.3"""
    return 0.71 * logp - 0.0061 * mw - 6.3


cs = pd.read_csv(os.path.join(PROC, "case_study_clean.csv"))
probs = pd.read_csv(os.path.join(OUT, "stan_probs.csv"))

rows = []
for _, r in cs.iterrows():
    m = Chem.MolFromSmiles(r["Smiles"])
    logp = Crippen.MolLogP(m)
    mw = Descriptors.MolWt(m)
    tpsa = Descriptors.TPSA(m)
    rotb = rdMolDescriptors.CalcNumRotatableBonds(m)
    rows.append({
        "ID": r["ID"], "series": r["series"],
        "MW": mw, "logP": logp, "TPSA": tpsa, "RotBond": rotb,
        "logKp": potts_guy_logkp(logp, mw),
        "skin_retention_pct": r["skin_retention_pct"],
        "dc50_nM": r["dc50_nM"], "dc50_censored": r["dc50_censored"],
        "dmax_pct": r["dmax_pct"],
        "MW_gt750": mw > 750, "TPSA_gt140": tpsa > 140,
    })
df = pd.DataFrame(rows)

P("=" * 70)
P("Phase 3 — Potts-Guy logKp + AD 외삽 + 2D 분리")
P("=" * 70)
P(f"\nlogKp 범위: {df['logKp'].min():.2f} ~ {df['logKp'].max():.2f} "
  f"(전신투과 proxy, 낮을수록 국소제 유리)")
P(f"B/C logKp 범위: {df[df.series.isin(['B','C'])]['logKp'].min():.2f} ~ "
  f"{df[df.series.isin(['B','C'])]['logKp'].max():.2f}")

# AD 외삽 (Potts-Guy/Flynn: 주로 MW<750 소분자)
n_mw = int(df["MW_gt750"].sum()); n_tpsa = int(df["TPSA_gt140"].sum())
n_out = int((df["MW_gt750"] | df["TPSA_gt140"]).sum())
P(f"\n[AD 외삽] MW>750: {n_mw}/29, TPSA>140: {n_tpsa}/29, "
  f"AD 밖(MW>750 또는 TPSA>140): {n_out}/29 = {n_out/29:.1%}")
bc = df[df.series.isin(["B", "C"])]
P(f"[AD 외삽] B/C 13개 중 MW>750: {int(bc['MW_gt750'].sum())}, "
  f"TPSA>140: {int(bc['TPSA_gt140'].sum())}")
P("⑤: B/C 전부 적합상한(~750) 경계 밖 → MW항 지배로 logKp 균일 압축, 변별력 상실(예상 결과).")

# 변별력: B/C logKp 표준편차
P(f"[변별력] B/C logKp std = {bc['logKp'].std():.3f} (작을수록 압축됨)")

# retention vs logKp 상관 (정성, bootstrap 생략 — n작음 명시)
sub = df.dropna(subset=["skin_retention_pct", "logKp"])
rho = spearmanr(sub["logKp"], sub["skin_retention_pct"]).correlation
P(f"\n[2D 분리] Spearman(logKp, retention) = {rho:.3f}, n={len(sub)} "
  f"(주 결과 아님; 잔류≠통과라 무·음상관이 정상)")

# --- [표3] -----------------------------------------------------------------
t3 = df[["ID", "series", "MW", "logP", "TPSA", "RotBond", "logKp",
         "MW_gt750", "TPSA_gt140", "skin_retention_pct"]].round(3)
t3.to_csv(os.path.join(OUT, "table3.csv"), index=False)
P("\n[저장] outputs/table3.csv")

# --- [그림3] 2D 분리 -------------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 6))
scolors = {"A": "#fbbc04", "B": "#4285f4", "C": "#ea4335"}
# 이상 사분면 음영: 낮은 logKp(좌) · 높은 retention(상)
ymax = df["skin_retention_pct"].max() * 1.1
xmid = df["logKp"].median()
ax.axhline(df["skin_retention_pct"].median(), ls=":", c="gray", lw=0.8)
ax.axvline(xmid, ls=":", c="gray", lw=0.8)
ax.axvspan(df["logKp"].min() - 0.2, xmid, ymin=0.5, ymax=1.0,
           color="green", alpha=0.06)
ax.text(df["logKp"].min(), ymax * 0.80,
        "이상 영역\n(낮은 전신투과·높은 피부잔류)", fontsize=8, color="green")
for sname in ["A", "B", "C"]:
    d = df[df.series == sname]
    active = d["dc50_nM"].notna() & d["dc50_censored"].isna() & (d["dc50_nM"] < 100)
    ax.scatter(d["logKp"], d["skin_retention_pct"], s=70, c=scolors[sname],
               edgecolor=["darkred" if a else "black" for a in active],
               linewidths=[2.0 if a else 0.6 for a in active],
               marker={"A": "o", "B": "s", "C": "^"}[sname], label=f"series {sname}")
for _, r in df.iterrows():
    if r["ID"] in ("C5", "C1", "C6", "B3"):
        ax.annotate(r["ID"], (r["logKp"], r["skin_retention_pct"]),
                    fontsize=8, xytext=(3, 3), textcoords="offset points")
ax.set_xlabel("predicted logKp (Potts-Guy) — 전신투과 proxy (낮을수록 유리)")
ax.set_ylabel("measured Skin Retention Rate (%)")
ax.set_title("2D 분리: 전신투과(logKp) vs 피부잔류(retention)   "
             "[적색 두꺼운 테두리=실측 DC50<100nM]", fontsize=11)
ax.legend(fontsize=8)
plt.tight_layout()
plt.savefig(os.path.join(OUT, "fig3_separation.png"), dpi=150)
P("[저장] outputs/fig3_separation.png")

# --- [그림3b] 2D 카멜레온 (TPSA vs RotBond) --------------------------------
fig, ax = plt.subplots(figsize=(7, 6))
for sname in ["A", "B", "C"]:
    d = df[df.series == sname]
    ax.scatter(d["RotBond"], d["TPSA"], s=70, c=scolors[sname],
               edgecolor="black", marker={"A": "o", "B": "s", "C": "^"}[sname],
               label=f"series {sname}")
ax.set_xlabel("Rotatable Bond Count")
ax.set_ylabel("TPSA (2D, conformer-independent)")
ax.set_title("카멜레온성 2D proxy (TPSA vs RotBond)\n"
             "주의: TPSA는 conformer 무관 고정값 — 3D EPSA는 부록(NICE)")
ax.legend(fontsize=8)
plt.tight_layout()
plt.savefig(os.path.join(OUT, "fig3b_chameleon.png"), dpi=150)
P("[저장] outputs/fig3b_chameleon.png")

with open(os.path.join(OUT, "phase3_log.txt"), "w") as f:
    f.write("\n".join(log))
P("\nPhase 3 완료.")
