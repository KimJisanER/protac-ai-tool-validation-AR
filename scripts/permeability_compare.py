"""
Caco-2 투과 예측 3-way 비교 (동일 실측 Caco-2, 모두 held-out 기준):
 (1) PROTAC-TS  : PROTAC 학습 ML — held-out CV는 protacts_holdout_cv.py / protacts_cv_check.py (Spearman≈0.78, R²≈0.5)
 (2) descriptor 선형식 : 물성(TPSA·cLogP·MW·HBD·RotB·HBA) MLR을 직접 적합 + compound LOGO
 (3) Potts-Guy : 학습 없는 소분자 피부 QSPR 공식 (자동 held-out)
임의 계수를 지어내지 않으며, (2)는 데이터에 직접 적합해 계수를 공개한다.
"""
import numpy as np, pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, Crippen
RDLogger.DisableLog('rdApp.*')
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import LeaveOneGroupOut
from sklearn.metrics import r2_score
from scipy.stats import spearmanr, pearsonr
TS='/home/kimjisan95/PROTAC-TS'
fea=pd.read_csv(f'{TS}/result/feature/test/fea_morganfp.csv'); raw=pd.read_csv(f'{TS}/data/PROTAC_Caco-2_data.csv')
fea=fea[fea['cate']=='PROTAC'].reset_index(drop=True)
fea['Smiles']=fea['Unnamed: 0'].map(lambda i: raw['Smiles'].iloc[int(i)])
def ik14(s):
    m=Chem.MolFromSmiles(str(s)); return Chem.MolToInchiKey(m)[:14] if m else None
def desc(s):
    m=Chem.MolFromSmiles(str(s))
    return [Descriptors.TPSA(m),Crippen.MolLogP(m),Descriptors.MolWt(m),
            Descriptors.NumHDonors(m),Descriptors.NumRotatableBonds(m),Descriptors.NumHAcceptors(m)]
names=['TPSA','cLogP','MW','HBD','RotB','HBA']
D=np.array([desc(s) for s in fea['Smiles']]); y=fea['log10Papp'].values; g=fea['Smiles'].map(ik14).values

# (3) Potts-Guy (피부 공식, 학습 없음)
pg=0.71*D[:,1]-0.0061*D[:,2]-6.3
print(f"[Potts-Guy] vs 실측 Caco-2: Spearman={spearmanr(pg,y)[0]:+.3f} | Pearson={pearsonr(pg,y)[0]:+.3f} "
      f"| AD밖(MW>750|logP>5) {((D[:,2]>750)|(D[:,1]>5)).mean()*100:.0f}%")
# (2) descriptor MLR — compound LOGO
yp=np.full(len(y),np.nan)
for tr,te in LeaveOneGroupOut().split(D,y,g):
    yp[te]=LinearRegression().fit(D[tr],y[tr]).predict(D[te])
print(f"[descriptor MLR] compound LOGO held-out: R2={r2_score(y,yp):+.3f} | Spearman={spearmanr(y,yp)[0]:+.3f}")
lr=LinearRegression().fit(D,y)
print("  전체적합 계수: "+" ".join(f"{n}:{c:+.4f}" for n,c in zip(names,lr.coef_))+f" | 절편 {lr.intercept_:+.3f} | in-sample R2={lr.score(D,y):.3f}")
print("[단일 descriptor Spearman] "+" ".join(f"{n}:{spearmanr(D[:,i],y)[0]:+.2f}" for i,n in enumerate(names)))
print("[PROTAC-TS] held-out CV Spearman≈0.78(LOGO)/0.70(5-fold), R²≈0.5  (protacts_cv_check.py)")
