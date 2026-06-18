"""
Phase 2 분석 — STAN 활성확률 vs 실측 concordance.

- stan_probs.csv (ID, stan_active_prob) + case_study_clean.csv (DC50/Dmax/leak) 병합
- 주지표: class-1 확률 vs 실측 DC50(역순)/Dmax Spearman 순위상관 + bootstrap 95% CI
- leakage(B3/B5) 사후 마스킹: 전체 vs leakage-free 둘 다 보고
- 산출물: outputs/table2.csv, outputs/fig2_ranking.png, outputs/phase2_log.txt
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import plotstyle  # noqa: F401  (한글 폰트)

ROOT = "/home/kimjisan95/ar_protac_project"
OUT = os.path.join(ROOT, "outputs")
PROC = os.path.join(ROOT, "data", "processed")

probs = pd.read_csv(os.path.join(OUT, "stan_probs.csv"))
cs = pd.read_csv(os.path.join(PROC, "case_study_clean.csv"))
df = cs.merge(probs[["ID", "stan_active_prob", "stan_pred"]], on="ID", how="left")

log = []
def P(*a):
    s = " ".join(str(x) for x in a); print(s); log.append(s)


def boot_spearman(x, y, n=2000, seed=0):
    """bootstrap 95% CI for Spearman rho."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    mask = ~(np.isnan(x) | np.isnan(y))
    x, y = x[mask], y[mask]
    if len(x) < 4:
        return np.nan, (np.nan, np.nan), len(x)
    rho = spearmanr(x, y).correlation
    rng = np.random.RandomState(seed)
    boots = []
    for _ in range(n):
        idx = rng.randint(0, len(x), len(x))
        if len(np.unique(x[idx])) < 2 or len(np.unique(y[idx])) < 2:
            continue
        boots.append(spearmanr(x[idx], y[idx]).correlation)
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return rho, (lo, hi), len(x)


P("=" * 70)
P("Phase 2 — STAN 활성확률 vs 실측 concordance (B/C 13개 활성 데이터)")
P("=" * 70)

bc = df[df["series"].isin(["B", "C"])].copy()
P(f"\nB/C 활성 데이터 행: {len(bc)} (DC50 보유 {bc['dc50_nM'].notna().sum()})")

# 주지표: 확률 vs DC50 (DC50 낮을수록 활성 → 음의 상관 기대; -log10 변환으로 정렬)
# censored '>1000' 등 제외하고 수치형만
num = bc[(bc["dc50_nM"].notna()) & (bc["dc50_censored"].isna())].copy()
P(f"\n[주지표A] 수치형 DC50 {len(num)}개로 Spearman(확률 vs -DC50):")
P(f"          대상: {dict(zip(num['ID'], num['dc50_nM']))}")
rho, ci, nn = boot_spearman(num["stan_active_prob"], -num["dc50_nM"])
P(f"          전체 rho={rho:.3f}, 95%CI[{ci[0]:.3f},{ci[1]:.3f}], n={nn}")

# leakage-free
numlf = num[~num["leaked_stan_train"] & ~num["leaked_old_db"]]
P(f"\n[주지표A-leakage-free] 누수(B3,B5) 제외 {len(numlf)}개: "
  f"{dict(zip(numlf['ID'], numlf['dc50_nM']))}")
rho2, ci2, nn2 = boot_spearman(numlf["stan_active_prob"], -numlf["dc50_nM"])
P(f"          rho={rho2:.3f}, 95%CI[{ci2[0]:.3f},{ci2[1]:.3f}], n={nn2}")

# 확률 vs Dmax
P(f"\n[주지표B] 확률 vs Dmax (B/C {bc['dmax_pct'].notna().sum()}개):")
rhod, cid, nnd = boot_spearman(bc["stan_active_prob"], bc["dmax_pct"])
P(f"          rho={rhod:.3f}, 95%CI[{cid[0]:.3f},{cid[1]:.3f}], n={nnd}")

# leakage-free 주 앵커 정성 일치
P("\n[앵커 concordance] leakage-free 주 앵커:")
for aid in ["C5", "C1", "C6"]:
    r = df[df["ID"] == aid].iloc[0]
    P(f"   {aid}: DC50={r['dc50_nM']}nM, Dmax={r['dmax_pct']}%, "
      f"STAN확률={r['stan_active_prob']:.3f}")
P("[누수 앵커 — 참고용(memorization 위험)]:")
for aid in ["B3", "B5"]:
    r = df[df["ID"] == aid].iloc[0]
    P(f"   {aid}: DC50={r['dc50_nM']}nM, Dmax={r['dmax_pct']}%, "
      f"STAN확률={r['stan_active_prob']:.3f}  [LEAKED]")

# --- [표2] -----------------------------------------------------------------
t2 = df[["ID", "series", "stan_active_prob", "stan_pred",
         "dc50_nM", "dc50_censored", "dmax_pct",
         "skin_retention_pct", "leaked_stan_train", "leaked_old_db"]].copy()
t2["leaked"] = t2["leaked_stan_train"] | t2["leaked_old_db"]
t2 = t2.round(3)
t2.to_csv(os.path.join(OUT, "table2.csv"), index=False)
P("\n[저장] outputs/table2.csv")

# --- [그림2] STAN 확률 랭킹 ------------------------------------------------
fig, ax = plt.subplots(figsize=(11, 5))
d = df.sort_values("stan_active_prob", ascending=False).reset_index(drop=True)
colors = {"A": "#9aa0a6", "B": "#4285f4", "C": "#ea4335"}
bar_colors = [colors[s] for s in d["series"]]
bars = ax.bar(range(len(d)), d["stan_active_prob"], color=bar_colors)
# 누수 표시(해치)
for i, (_, r) in enumerate(d.iterrows()):
    if r["leaked_stan_train"] or r["leaked_old_db"]:
        bars[i].set_hatch("////")
        bars[i].set_edgecolor("black")
# 활성(DC50<100 수치형) 별표
for i, (_, r) in enumerate(d.iterrows()):
    if pd.notna(r["dc50_nM"]) and pd.isna(r["dc50_censored"]) and r["dc50_nM"] < 100:
        ax.text(i, r["stan_active_prob"] + 0.02, "★", ha="center", color="darkred")
ax.set_xticks(range(len(d)))
ax.set_xticklabels(d["ID"], rotation=90, fontsize=8)
ax.set_ylabel("STAN class-1 (active) probability")
ax.set_title("PROTAC-STAN activity-probability ranking (29 case-study compounds)\n"
             "color=series(A/B/C), hatch=leaked(B3/B5), ★=measured DC50<100nM")
from matplotlib.patches import Patch
legend = [Patch(facecolor=colors["A"], label="A (no activity data)"),
          Patch(facecolor=colors["B"], label="B"),
          Patch(facecolor=colors["C"], label="C"),
          Patch(facecolor="white", edgecolor="black", hatch="////", label="leaked")]
ax.legend(handles=legend, fontsize=8)
ax.set_ylim(0, 1.08)
plt.tight_layout()
plt.savefig(os.path.join(OUT, "fig2_ranking.png"), dpi=150)
P("[저장] outputs/fig2_ranking.png")

with open(os.path.join(OUT, "phase2_log.txt"), "w") as f:
    f.write("\n".join(log))
P("\nPhase 2 분석 완료.")
