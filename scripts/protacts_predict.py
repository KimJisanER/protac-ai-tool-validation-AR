"""
PROTAC-TS(세포막 투과도 TabPFN 모델)로 case study 29개 Caco-2 투과도 예측 →
피부 잔류·Potts-Guy logKp·STAN 활성과의 (비)상관 = '3중 도메인 구분' 실증.
env: protac_ts (TabPFN). 모델: ~/PROTAC-TS/model/pre_model/PROTAC/model_all.pkl
산출: outputs/protacts_caco2.csv, outputs/protacts_log.txt, outputs/fig_protacts.png
"""
import os, pickle, numpy as np, pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
RDLogger.DisableLog('rdApp.*')
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.stats import spearmanr

ROOT='/home/kimjisan95/ar_protac_project'; OUT=f'{ROOT}/outputs'
MODEL=os.path.expanduser('~/PROTAC-TS/model/pre_model/PROTAC/model_all.pkl')
log=[]; P=lambda *a:(print(*a), log.append(' '.join(str(x) for x in a)))

def calc_MorganCount(mol, r=2, dimension=500):
    info={}; AllChem.GetMorganFingerprint(mol, r, bitInfo=info)
    cl=[0]*dimension
    for k in info: cl[k%dimension]+=len(info[k])
    return cl

P('='*70); P('PROTAC-TS: case study 29개 Caco-2 세포막 투과도 예측'); P('='*70)
with open(MODEL,'rb') as f: model=pickle.load(f)
P(f'[모델] TabPFN 로드: {MODEL}')

cs=pd.read_csv(f'{ROOT}/data/processed/case_study_clean.csv')
X=np.array([calc_MorganCount(Chem.MolFromSmiles(s),2,500) for s in cs['Smiles']],dtype=np.float32)
cs['caco2_log10Papp_pred']=model.predict(X)            # 모델 native 출력 = log10(μcm/s)
cs['caco2_Papp_pred_nm_s']=(10**cs['caco2_log10Papp_pred'])*10   # 역변환(nm/s)
P(f'[예측] 29개 log10Papp 범위 {cs["caco2_log10Papp_pred"].min():.2f}~{cs["caco2_log10Papp_pred"].max():.2f}')

# 다른 도메인 지표 병합
t3=pd.read_csv(f'{OUT}/table3.csv')[['ID','logKp','MW','TPSA']]
t2=pd.read_csv(f'{OUT}/table2.csv')[['ID','stan_active_prob']]
d=cs.merge(t3,on='ID',how='left').merge(t2,on='ID',how='left')
d.to_csv(f'{OUT}/protacts_caco2.csv',index=False)

def sp(a,b,name):
    m=d[[a,b]].dropna()
    if len(m)>=5:
        r,p=spearmanr(m[a],m[b]); P(f'  {name}: Spearman={r:+.3f} (p={p:.3f}, n={len(m)})'); return r
    return None
P('\n[3중 도메인 구분 — 예측 Caco-2 세포막투과 vs 다른 축]')
sp('caco2_log10Papp_pred','skin_retention_pct','세포막투과 vs 피부잔류(retention)')
sp('caco2_log10Papp_pred','logKp','세포막투과 vs Potts-Guy logKp(피부 전신투과)')
sp('caco2_log10Papp_pred','stan_active_prob','세포막투과 vs STAN 활성확률')
sp('caco2_log10Papp_pred','dmax_pct','세포막투과 vs Dmax')
P('\n[계열별 예측 Caco-2 log10Papp 평균]')
for s,g in d.groupby('series'): P(f'  {s}: {g["caco2_log10Papp_pred"].mean():+.2f} (n={len(g)})')

# 그림: 세포막투과(예측) vs 피부잔류(실측)
fig,ax=plt.subplots(1,2,figsize=(9,4))
cmap={'A':'#1f77b4','B':'#ff7f0e','C':'#2ca02c'}
for s,g in d.groupby('series'):
    ax[0].scatter(g['caco2_log10Papp_pred'],g['skin_retention_pct'],c=cmap.get(s,'#999'),label=f'{s}',s=28)
r1=spearmanr(d['caco2_log10Papp_pred'],d['skin_retention_pct'])[0]
ax[0].set_xlabel('Predicted Caco-2 log10Papp (cell-membrane)');ax[0].set_ylabel('Measured skin retention (%)')
ax[0].set_title(f'Cell-membrane vs skin retention\n(Spearman={r1:+.2f}, n=29)');ax[0].legend(title='series',fontsize=8)
dd=d.dropna(subset=['logKp'])
for s,g in dd.groupby('series'):
    ax[1].scatter(g['caco2_log10Papp_pred'],g['logKp'],c=cmap.get(s,'#999'),label=f'{s}',s=28)
r2=spearmanr(dd['caco2_log10Papp_pred'],dd['logKp'])[0]
ax[1].set_xlabel('Predicted Caco-2 log10Papp');ax[1].set_ylabel('Potts-Guy logKp (skin flux)')
ax[1].set_title(f'Cell-membrane vs skin permeation\n(Spearman={r2:+.2f}, n=29)');ax[1].legend(title='series',fontsize=8)
plt.tight_layout();plt.savefig(f'{OUT}/fig_protacts.png',dpi=150);plt.close()
P('\n[그림] outputs/fig_protacts.png  [저장] outputs/protacts_caco2.csv')
open(f'{OUT}/protacts_log.txt','w').write('\n'.join(log)+'\n')
