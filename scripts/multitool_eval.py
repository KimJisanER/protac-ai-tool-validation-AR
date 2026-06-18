"""
STAN vs Ribes head-to-head — 동일 doubly-held-out 786개, 다중 렌즈(정직 버전).
적대적 검증 반영: (1) 라벨 비중립성(STAN-strict vs Ribes-native), (2) KRAS 단일클래스
블록 교란 → KRAS 제외 + 타겟별 macro, (3) 임계-무관 지표(AUROC/AUPR) 위주, (4) 불일치는
κ·Spearman 주지표. 산출: outputs/multitool_metrics.csv, multitool_pertarget.csv,
multitool_agreement.txt, fig_multitool.png
"""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score, average_precision_score, cohen_kappa_score
from scipy.stats import spearmanr
OUT='/home/kimjisan95/ar_protac_project/outputs'; DP='/home/kimjisan95/ar_protac_project/data/processed'
rng=np.random.default_rng(21332); log=[]; P=lambda *a:(print(*a), log.append(' '.join(str(x) for x in a)))

rib=pd.read_csv(f'{OUT}/ribes_probs.csv')
stan=pd.read_csv(f'{OUT}/stan_ts_new_probs.csv')[['ik14','stan_active_prob']]
tgt=pd.read_csv(f'{DP}/ribes_input.csv')[['ik14','Target','dc50_nM','dmax_pct']]
d=rib.merge(stan,on='ik14').merge(tgt[['ik14','Target']],on='ik14',how='left')
P('='*70); P(f'STAN vs Ribes — doubly-held-out 동일 {len(d)}개 (다중 렌즈)'); P('='*70)

# 두 가지 GT 라벨
y_strict=d['active'].astype(int).values                       # STAN 학습 라벨: DC50<100 & Dmax>=80
y_native=((d['dmax_pct']>=60)&(d['dc50_nM']<=1000)).astype(int).values  # Ribes 학습 라벨: Dmax>=60 & DC50<=1000
P(f'양성률 — STAN-strict GT: {y_strict.mean():.1%}({y_strict.sum()}/{len(d)}) | Ribes-native GT: {y_native.mean():.1%}({y_native.sum()}/{len(d)})')

probs={'PROTAC-STAN (GNN+ESM, 단일 동결)':d['stan_active_prob'].values,
       'Ribes standard (FP+XGB/MLP, random-CV 앙상블)':d['ribes_proba_standard'].values,
       'Ribes target-split (group-CV 앙상블)':d['ribes_proba_target'].values}

def bootci(y,s,n=2000):
    idx=np.arange(len(y));v=[]
    for _ in range(n):
        b=rng.choice(idx,len(idx),replace=True)
        if len(np.unique(y[b]))<2: continue
        try:v.append(roc_auc_score(y[b],s[b]))
        except:pass
    return (round(np.percentile(v,2.5),3),round(np.percentile(v,97.5),3)) if v else (np.nan,np.nan)

kras=d['Target'].astype(str).str.contains('KRAS',case=False,na=False).values
P(f'\n[교란] KRAS 화합물 {int(kras.sum())}/{len(d)} (양성 {int(y_strict[kras].sum())}) — 단일클래스 블록')

rows=[]
for name,s in probs.items():
    lo,hi=bootci(y_strict,s)
    r={'model':name,
       'AUROC_strictGT':round(roc_auc_score(y_strict,s),3),'CI':f'[{lo},{hi}]',
       'AUPR_strictGT':round(average_precision_score(y_strict,s),3),
       'AUROC_nativeGT':round(roc_auc_score(y_native,s),3),
       'AUROC_exKRAS_strictGT':round(roc_auc_score(y_strict[~kras],s[~kras]),3),
       'mean_prob_on_KRAS(all-neg)':round(float(s[kras].mean()),3)}
    rows.append(r)
M=pd.DataFrame(rows); M.to_csv(f'{OUT}/multitool_metrics.csv',index=False)
P('\n[성능 — 다중 렌즈]'); P(M.to_string(index=False))

# 타겟별 AUROC (양·음 둘 다 있는 타겟) + macro 평균
P('\n[타겟별 AUROC — 양·음 둘 다 존재하는 타겟, strict GT]')
pt_rows=[]
for u,g in d.groupby('Target'):
    ys=y_strict[d['Target'].values==u]
    if len(np.unique(ys))<2: continue
    rec={'Target':u,'n':len(g),'pos':int(ys.sum())}
    for name,s in probs.items():
        ss=s[d['Target'].values==u]
        rec[name.split(' (')[0]]=round(roc_auc_score(ys,ss),3)
    pt_rows.append(rec)
PT=pd.DataFrame(pt_rows); PT.to_csv(f'{OUT}/multitool_pertarget.csv',index=False)
P(PT.to_string(index=False))
P(f'\n[타겟별 macro 평균 AUROC] ({len(PT)}개 양·음 타겟)')
for name in probs:
    col=name.split(' (')[0]; P(f'  {col}: {PT[col].mean():.3f}')

# 불일치(주지표 κ·ρ)
P('\n[모델 간 (불)일치 — 라벨-무관, 견고]')
sc=(d['stan_active_prob']>=.5).astype(int); rc=(d['ribes_proba_standard']>=.5).astype(int)
rho,_=spearmanr(d['stan_active_prob'],d['ribes_proba_standard'])
P(f'  Spearman(prob) = {rho:.3f}  |  Cohen κ(0.5컷) = {cohen_kappa_score(sc,rc):.3f}  |  호출일치율 = {(sc==rc).mean():.1%}(0.5컷·불균형 의존, 보조)')
dis=sc!=rc
P(f'  불일치 {int(dis.sum())}건 중 STAN 정답 {((sc[dis]==y_strict[dis]).mean()):.1%} / Ribes 정답 {((rc[dis]==y_strict[dis]).mean()):.1%}')

# 그림
fig,ax=plt.subplots(1,2,figsize=(9.5,4))
for lab,c in [(1,'#2ca02c'),(0,'#999')]:
    m=y_strict==lab
    ax[0].scatter(d['stan_active_prob'][m],d['ribes_proba_standard'][m],s=10,alpha=.5,c=c,label=f'GT(strict)={lab}')
ax[0].axhline(.5,color='k',lw=.6,ls='--');ax[0].axvline(.5,color='k',lw=.6,ls='--')
ax[0].set_xlabel('STAN prob');ax[0].set_ylabel('Ribes prob (standard)')
ax[0].set_title(f'STAN vs Ribes (n={len(d)}, ρ={rho:.2f}, κ=0.13)');ax[0].legend(fontsize=7)
labels=['STAN','Ribes·std','Ribes·tgt']
x=np.arange(3); w=0.38
au_s=[r['AUROC_strictGT'] for r in rows]; au_x=[r['AUROC_exKRAS_strictGT'] for r in rows]
ax[1].bar(x-w/2,au_s,w,label='pooled (strict GT)',color='#1f77b4')
ax[1].bar(x+w/2,au_x,w,label='KRAS excluded',color='#aec7e8')
ax[1].axhline(.5,color='k',ls='--',lw=.8);ax[1].set_xticks(x);ax[1].set_xticklabels(labels,fontsize=8)
ax[1].set_ylim(0,1);ax[1].set_ylabel('AUROC');ax[1].set_title('AUROC: pooled vs KRAS-excluded');ax[1].legend(fontsize=7)
for i,(a,b) in enumerate(zip(au_s,au_x)):
    ax[1].text(i-w/2,a+.02,f'{a:.2f}',ha='center',fontsize=7);ax[1].text(i+w/2,b+.02,f'{b:.2f}',ha='center',fontsize=7)
plt.tight_layout();plt.savefig(f'{OUT}/fig_multitool.png',dpi=150);plt.close()
P('\n[그림] outputs/fig_multitool.png  [저장] multitool_metrics.csv, multitool_pertarget.csv')
open(f'{OUT}/multitool_agreement.txt','w').write('\n'.join(log)+'\n')
