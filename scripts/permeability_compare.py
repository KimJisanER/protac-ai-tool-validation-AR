"""
투과 도구 머리맞대기: 동일 실측 Caco-2(PROTAC-TS 학습 데이터의 ground truth)에 대해
 - PROTAC-TS : held-out CV(별도 protacts_holdout_cv.py) Spearman≈0.77~0.81, R²≈0.5
 - Potts-Guy : 학습 없는 소분자 피부 QSPR 공식(자동 held-out)
공정 비교를 위해 둘 다 '같은 실측 Caco-2'와의 순위상관으로 평가.
출력(콘솔): Potts-Guy Spearman/Pearson, AD 밖 비율.
"""
import numpy as np, pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, Crippen
RDLogger.DisableLog('rdApp.*')
from scipy.stats import spearmanr, pearsonr
TS='/home/kimjisan95/PROTAC-TS'
fea=pd.read_csv(f'{TS}/result/feature/test/fea_morganfp.csv'); raw=pd.read_csv(f'{TS}/data/PROTAC_Caco-2_data.csv')
fea=fea[fea['cate']=='PROTAC'].reset_index(drop=True)
fea['Smiles']=fea['Unnamed: 0'].map(lambda i: raw['Smiles'].iloc[int(i)])
def feats(s):
    m=Chem.MolFromSmiles(str(s)); return (Descriptors.MolWt(m), Crippen.MolLogP(m)) if m else (None,None)
fea[['MW','logP']]=fea['Smiles'].apply(lambda s: pd.Series(feats(s)))
fea['pottsguy_logKp']=0.71*fea['logP']-0.0061*fea['MW']-6.3   # Potts & Guy 1992 (skin)
y=fea['log10Papp'].values                                     # 실측 Caco-2
ad_out=((fea.MW>750)|(fea.logP>5)).mean()
print(f"n={len(fea)} | MW {fea.MW.min():.0f}~{fea.MW.max():.0f} | Potts-Guy AD 밖(MW>750|logP>5) {ad_out*100:.0f}%")
print(f"Potts-Guy vs 실측 Caco-2:  Spearman={spearmanr(fea.pottsguy_logKp,y)[0]:+.3f} | Pearson={pearsonr(fea.pottsguy_logKp,y)[0]:+.3f}")
print("PROTAC-TS  vs 실측 Caco-2:  held-out CV Spearman≈0.77~0.81, R²≈0.5 (protacts_holdout_cv.py)")
