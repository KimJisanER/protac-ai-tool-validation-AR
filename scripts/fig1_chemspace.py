"""
Phase 1 — Chemical space + bRo5 overlay ([그림1]).

3집합: 구DB AR(P10275), 신DB AR(P10275), case study(A/B/C).
좌: Morgan FP(r=2,2048) PCA, 색=데이터셋, case study는 A/B/C 마커.
우: MW vs TPSA, Ro5(MW=500,TPSA=140) + bRo5 상한(MW=1000,TPSA=250, Doak 2014) 참조선.

⑤ talking point: PCA 분산설명률 낮으면 거리 비해석 → 설명률 출력.
산출물: outputs/fig1_chemspace.png, outputs/phase1_log.txt
"""
import os
import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from sklearn.decomposition import PCA
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


def fp_array(smiles_list):
    fps, mws, tpsas, keep = [], [], [], []
    for s in smiles_list:
        m = Chem.MolFromSmiles(s) if isinstance(s, str) else None
        if m is None:
            continue
        fp = AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048)
        arr = np.zeros((2048,), dtype=np.int8)
        Chem.DataStructs.ConvertToNumpyArray(fp, arr)
        fps.append(arr)
        mws.append(Descriptors.MolWt(m))
        tpsas.append(Descriptors.TPSA(m))
        keep.append(s)
    return np.array(fps), np.array(mws), np.array(tpsas)


ar_old = pd.read_csv(os.path.join(PROC, "ar_old.csv"))
ar_new = pd.read_csv(os.path.join(PROC, "ar_new.csv"))
cs = pd.read_csv(os.path.join(PROC, "case_study_clean.csv"))

fp_old, mw_old, tp_old = fp_array(ar_old["Protac_SMILES"].tolist())
fp_new, mw_new, tp_new = fp_array(ar_new["Smiles"].tolist())
fp_cs, mw_cs, tp_cs = fp_array(cs["Smiles"].tolist())
P(f"FP 행수: 구DB AR {len(fp_old)}, 신DB AR {len(fp_new)}, case study {len(fp_cs)}")

X = np.vstack([fp_old, fp_new, fp_cs])
labels = (["old"] * len(fp_old) + ["new"] * len(fp_new) + ["cs"] * len(fp_cs))
labels = np.array(labels)

pca = PCA(n_components=2, random_state=0)
Z = pca.fit_transform(X)
evr = pca.explained_variance_ratio_
P(f"PCA 분산설명률: PC1={evr[0]:.1%}, PC2={evr[1]:.1%}, 합={evr[:2].sum():.1%}")
P("⑤: 분산설명률이 낮으면 2D 거리는 정성적으로만 해석.")

# --- plot -----------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

ax = axes[0]
ax.scatter(Z[labels == "old", 0], Z[labels == "old", 1], s=12, alpha=0.4,
           c="#9aa0a6", label=f"old DB AR (n={len(fp_old)})")
ax.scatter(Z[labels == "new", 0], Z[labels == "new", 1], s=12, alpha=0.4,
           c="#34a853", label=f"new DB AR (n={len(fp_new)})")
# case study by series
cs_Z = Z[labels == "cs"]
series = cs["series"].values
markers = {"A": "o", "B": "s", "C": "^"}
scolors = {"A": "#fbbc04", "B": "#4285f4", "C": "#ea4335"}
for sname in ["A", "B", "C"]:
    msk = series == sname
    ax.scatter(cs_Z[msk, 0], cs_Z[msk, 1], s=80, marker=markers[sname],
               edgecolor="black", c=scolors[sname], label=f"case {sname}", zorder=5)
ax.set_xlabel(f"PC1 ({evr[0]:.1%})")
ax.set_ylabel(f"PC2 ({evr[1]:.1%})")
ax.set_title("Chemical space (Morgan FP r=2, 2048bit, PCA)")
ax.legend(fontsize=8)

ax = axes[1]
ax.scatter(mw_old, tp_old, s=12, alpha=0.35, c="#9aa0a6", label="old DB AR")
ax.scatter(mw_new, tp_new, s=12, alpha=0.35, c="#34a853", label="new DB AR")
for sname in ["A", "B", "C"]:
    msk = series == sname
    ax.scatter(mw_cs[msk], tp_cs[msk], s=80, marker=markers[sname],
               edgecolor="black", c=scolors[sname], label=f"case {sname}", zorder=5)
# Ro5 / bRo5 참조선
ax.axvline(500, ls="--", c="gray", lw=1); ax.axhline(140, ls="--", c="gray", lw=1)
ax.axvline(1000, ls=":", c="purple", lw=1.2); ax.axhline(250, ls=":", c="purple", lw=1.2)
ax.text(505, ax.get_ylim()[1]*0.95, "Ro5 MW=500", fontsize=7, color="gray")
ax.text(1005, 150, "bRo5 MW≈1000\n(Doak 2014)", fontsize=7, color="purple")
ax.text(50, 143, "TPSA=140", fontsize=7, color="gray")
ax.set_xlabel("Molecular Weight")
ax.set_ylabel("TPSA")
ax.set_title("Property space — Ro5 / bRo5 (Doak 2014) overlay")
ax.legend(fontsize=8)

plt.tight_layout()
plt.savefig(os.path.join(OUT, "fig1_chemspace.png"), dpi=150)
P("[저장] outputs/fig1_chemspace.png")
with open(os.path.join(OUT, "phase1_log.txt"), "w") as f:
    f.write("\n".join(log))
P("Phase 1 완료.")
