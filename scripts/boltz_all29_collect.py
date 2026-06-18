"""
29개 Boltz-2 ternary 예측 신뢰도 집계 + 실측 활성과의 관계.
출력: outputs/boltz_all29_confidence.csv, outputs/boltz_all29_log.txt
"""
import json, glob, os, numpy as np, pandas as pd
from scipy.stats import spearmanr
AP='/home/kimjisan95/ar_protac_project'
PRED='/home/kimjisan95/PROTAC_MTL_v5/boltz_all29/out/boltz_results_inputs/predictions'
log=[]; P=lambda *a:(print(*a), log.append(' '.join(str(x) for x in a)))

rows=[]
for d in sorted(glob.glob(f'{PRED}/*')):
    cid=os.path.basename(d)
    cj=glob.glob(f'{d}/confidence_{cid}_model_0.json')
    if not cj: P(f'  [경고] {cid} confidence 없음'); continue
    c=json.load(open(cj[0]))
    rows.append({'ID':cid,**{k:round(v,3) for k,v in c.items() if isinstance(v,(int,float))}})
b=pd.DataFrame(rows)
P(f'집계된 예측: {len(b)}/29')

# 실측 활성·STAN 병합
cs=pd.read_csv(f'{AP}/data/processed/case_study_clean.csv')[['ID','series','dc50_nM','dmax_pct']]
t2=pd.read_csv(f'{AP}/outputs/table2.csv')[['ID','stan_active_prob','leaked']]
b=b.merge(cs,on='ID',how='left').merge(t2,on='ID',how='left')
cols=['ID','series','confidence_score','ptm','iptm','ligand_iptm','protein_iptm','complex_plddt','dc50_nM','dmax_pct','stan_active_prob','leaked']
b=b[[c for c in cols if c in b.columns]].sort_values(['series','ID'])
b.to_csv(f'{AP}/outputs/boltz_all29_confidence.csv',index=False)
P('\n[표] 29개 Boltz 신뢰도 (요약)')
P(b.to_string(index=False))

P('\n[계열별 평균]')
for s,g in b.groupby('series'):
    P(f'  {s}(n={len(g)}): confidence {g.confidence_score.mean():.3f}, ipTM {g.iptm.mean():.3f}, ligand_ipTM {g.ligand_iptm.mean():.3f}, protein_ipTM {g.protein_iptm.mean():.3f}')

# 예측 신뢰도가 실측 분해 활성과 관계있나? (B/C 수치형 DC50만)
num=b.dropna(subset=['dc50_nM']); num=num[num['dc50_nM']<1000]
if len(num)>=4:
    for met in ['iptm','protein_iptm','ligand_iptm','confidence_score']:
        r,p=spearmanr(num[met], -num['dc50_nM'])
        P(f'  Spearman({met} vs -DC50, n={len(num)}) = {r:+.3f} (p={p:.2f})')
P('\n[저장] outputs/boltz_all29_confidence.csv')
open(f'{AP}/outputs/boltz_all29_log.txt','w').write('\n'.join(log)+'\n')
