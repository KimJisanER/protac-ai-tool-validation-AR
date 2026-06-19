"""
PROTAC-TS Caco-2 CV 재검증:
 (a) 중복(같은 ik14 여러 행) 실체 확인 — SMILES가 완전 동일한지 vs 골격만 같은지
 (b) 더 엄밀/결정적 CV: compound LOGO, scaffold LOGO (leave-one-group-out)
     + 반복 grouped 5-fold(여러 seed) 로 단일분할 변동성 확인
"""
import numpy as np, pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem.Scaffolds import MurckoScaffold
RDLogger.DisableLog('rdApp.*')
from sklearn.model_selection import LeaveOneGroupOut, GroupShuffleSplit
from sklearn.metrics import r2_score
from scipy.stats import spearmanr
from tabpfn import TabPFNRegressor

TS='/home/kimjisan95/PROTAC-TS'
fea=pd.read_csv(f'{TS}/result/feature/test/fea_morganfp.csv')
raw=pd.read_csv(f'{TS}/data/PROTAC_Caco-2_data.csv')
fea=fea[fea['cate']=='PROTAC'].reset_index(drop=True)
fea['Smiles']=fea['Unnamed: 0'].map(lambda i: raw['Smiles'].iloc[int(i)])
def ik14(s):
    m=Chem.MolFromSmiles(str(s)); return Chem.MolToInchiKey(m)[:14] if m else None
def scaf(s):
    m=Chem.MolFromSmiles(str(s))
    try: return MurckoScaffold.MurckoScaffoldSmiles(mol=m) if m else None
    except Exception: return None
fea['ik']=fea['Smiles'].map(ik14); fea['scaf']=fea['Smiles'].map(scaf)
FP=[c for c in fea.columns if c.startswith('FP')]
X=fea[FP].values.astype(np.float32); y=fea['log10Papp'].values.astype(np.float32)

print(f"행 {len(fea)} | 고유 ik14 {fea['ik'].nunique()} | 고유 scaffold {fea['scaf'].nunique()}")
# (a) 중복 실체
vc=fea['ik'].value_counts()
dup=vc[vc>1]
print(f"\n[중복] 2회 이상 측정된 화합물(ik14): {len(dup)}개, 관련 행 {int(dup.sum())}개")
print("  행수 분포:", dup.value_counts().to_dict())
for ik in list(dup.index)[:3]:
    sub=fea[fea['ik']==ik]
    smis=sub['Smiles'].unique()
    print(f"  ik={ik}: {len(sub)}행 | 고유 SMILES {len(smis)}개 | log10Papp={sub['log10Papp'].round(2).tolist()}")
    print(f"     SMILES 동일? {'예(완전 동일)' if len(smis)==1 else '아니오(이성질체/염 변형)'}")
# scaffold가 2개 화합물 묶는 경우 확인
sc=fea.groupby('scaf')['ik'].nunique(); multi=sc[sc>1]
print(f"\n[scaffold] 2개 이상 화합물을 묶는 골격: {len(multi)}개 → {dict(list(multi.items())[:5])}")

import os as _os
_DEV=_os.environ.get('TABPFN_DEVICE','cuda')
def tab(): return TabPFNRegressor(random_state=42, ignore_pretraining_limits=True, device=_DEV)
def run(splits, name):
    yp=np.full(len(y),np.nan)
    for tr,te in splits:
        reg=tab(); reg.fit(X[tr],y[tr]); yp[te]=reg.predict(X[te])
    m=~np.isnan(yp); return r2_score(y[m],yp[m]), spearmanr(y[m],yp[m])[0]

print("\n[결정적 LOGO (leave-one-group-out)]")
for col,nm in [('ik','compound LOGO (화합물 1개씩 통째 제거)'),('scaf','scaffold LOGO (골격 1개씩 통째 제거)')]:
    r2,rho=run(list(LeaveOneGroupOut().split(X,y,fea[col])), nm)
    print(f"  {nm}: R2={r2:+.3f} | Spearman={rho:+.3f}")

print("\n[반복 grouped 5-fold (seed 0~4) — 단일분할 변동성]")
for col,nm in [('ik','compound'),('scaf','scaffold')]:
    r2s,rhos=[],[]
    for sd in range(5):
        gss=GroupShuffleSplit(n_splits=5, test_size=0.2, random_state=sd)
        # pool out-of-sample per seed
        yp=np.full(len(y),np.nan)
        for tr,te in gss.split(X,y,fea[col]):
            reg=tab(); reg.fit(X[tr],y[tr]); yp[te]=reg.predict(X[te])
        m=~np.isnan(yp); r2s.append(r2_score(y[m],yp[m])); rhos.append(spearmanr(y[m],yp[m])[0])
    print(f"  {nm} group: R2={np.mean(r2s):+.3f}±{np.std(r2s):.3f} | Spearman={np.mean(rhos):+.3f}±{np.std(rhos):.3f}")
