"""
PROTAC-TS(Caco-2 투과 TabPFN)를 위한 hold-out 평가 — 직접 train/test 분할.
PROTAC-TS 학습셋이 PROTAC-DB Caco-2를 거의 전부 포함(61/62)해 외부 held-out이 없으므로,
같은 데이터를 분할 기준을 달리해 재학습·재평가한다(전향 일반화 추정).
 1) row-level LOOCV  = PROTAC-TS 자기 보고 재현(중복 화합물 누수 포함)
 2) compound group CV = 같은 화합물 누수 제거(ik14 그룹)
 3) scaffold group CV = 새 골격 일반화(Murcko 그룹)  ← 진짜 hold-out
env: protac_ts (TabPFN 2.0.3). 산출: outputs/protacts_holdout_cv.csv, outputs/protacts_holdout_log.txt
"""
import numpy as np, pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem.Scaffolds import MurckoScaffold
RDLogger.DisableLog('rdApp.*')
from sklearn.model_selection import GroupKFold, LeaveOneOut
from sklearn.metrics import r2_score
from scipy.stats import spearmanr
from tabpfn import TabPFNRegressor

TS='/home/kimjisan95/PROTAC-TS'; OUT='/home/kimjisan95/ar_protac_project/outputs'
log=[]; P=lambda *a:(print(*a), log.append(' '.join(str(x) for x in a)))

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
P(f'PROTAC-TS Caco-2 데이터: {len(fea)}개 측정(행) | 고유 화합물(ik14) {fea["ik"].nunique()} | 고유 Murcko scaffold {fea["scaf"].nunique()}')

def tab(): return TabPFNRegressor(random_state=42, ignore_pretraining_limits=True, device='cpu')
def run(splits, name):
    yp=np.full(len(y), np.nan)
    for tr,te in splits:
        reg=tab(); reg.fit(X[tr], y[tr]); yp[te]=reg.predict(X[te])
    msk=~np.isnan(yp)
    r2=r2_score(y[msk], yp[msk]); rho=spearmanr(y[msk], yp[msk])[0]
    P(f'  {name}: R2={r2:+.3f} | Spearman={rho:+.3f} | n(라벨)={int(msk.sum())}')
    return dict(eval=name, R2=round(r2,3), Spearman=round(rho,3), n=int(msk.sum()))

P('\n[PROTAC-TS Caco-2 — 분할 기준별 hold-out 성능]')
rows=[]
rows.append(run(list(LeaveOneOut().split(X)),                       'row-level LOOCV (PROTAC-TS 재현, 중복화합물 누수 포함)'))
rows.append(run(list(GroupKFold(5).split(X,y,groups=fea['ik'])),    'compound group 5-fold (같은 화합물 누수 제거)'))
rows.append(run(list(GroupKFold(5).split(X,y,groups=fea['scaf'])),  'scaffold group 5-fold (새 골격 일반화 = 진짜 hold-out)'))
pd.DataFrame(rows).to_csv(f'{OUT}/protacts_holdout_cv.csv', index=False)
P(f'\n[저장] {OUT}/protacts_holdout_cv.csv')
open(f'{OUT}/protacts_holdout_log.txt','w').write('\n'.join(log)+'\n')
