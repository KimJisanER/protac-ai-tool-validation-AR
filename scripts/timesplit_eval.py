"""
Phase TS-C — time-split prospective 성능 평가 + memorization 대조 + 타겟별 + 그림.
입력: outputs/stan_ts_{new,seen}_probs.csv, data/processed/timesplit_compounds.csv
산출: outputs/timesplit_metrics.csv, timesplit_per_target.csv,
      figs: ts_roc.png, ts_pr.png, ts_probdist.png, ts_per_target.png
      outputs/timesplit_eval_log.txt
"""
import os, numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.metrics import (roc_auc_score, average_precision_score, f1_score,
                             matthews_corrcoef, balanced_accuracy_score, accuracy_score,
                             roc_curve, precision_recall_curve)
from scipy.stats import spearmanr

ROOT='/home/kimjisan95/ar_protac_project'; OUT=os.path.join(ROOT,'outputs')
rng=np.random.default_rng(21332)
log=[]; P=lambda *a:(print(*a), log.append(' '.join(str(x) for x in a)))

newp=pd.read_csv(os.path.join(OUT,'stan_ts_new_probs.csv'))
seenp=pd.read_csv(os.path.join(OUT,'stan_ts_seen_probs.csv'))
m1=pd.read_csv(os.path.join(ROOT,'data/processed/timesplit_compounds.csv'))[['ik14','dc50_nM','dmax_pct','Target']]
m2=pd.read_csv(os.path.join(ROOT,'data/processed/timesplit_leaked.csv'))[['ik14','dc50_nM','dmax_pct','Target']]
meta=pd.concat([m1,m2]).drop_duplicates('ik14')
newp=newp.merge(m1,on='ik14',how='left')
seenp=seenp.merge(m2,on='ik14',how='left')

def boot_ci(y,s,fn,n=2000):
    y=np.asarray(y); s=np.asarray(s); vals=[]
    idx=np.arange(len(y))
    for _ in range(n):
        b=rng.choice(idx,len(idx),replace=True)
        if len(np.unique(y[b]))<2: continue
        try: vals.append(fn(y[b],s[b]))
        except Exception: pass
    if not vals: return (np.nan,np.nan)
    return (float(np.percentile(vals,2.5)), float(np.percentile(vals,97.5)))

def metrics(df,name):
    y=df['label'].astype(int).values; s=df['stan_active_prob'].values; pred=(s>=0.5).astype(int)
    row={'cohort':name,'n':len(y),'pos':int(y.sum()),'neg':int((y==0).sum()),'pos_rate':round(y.mean(),3)}
    def add(k,val,ci):
        row[k]=round(val,3); row[k+'_lo']=round(ci[0],3); row[k+'_hi']=round(ci[1],3)
    add('AUROC', roc_auc_score(y,s), boot_ci(y,s,roc_auc_score))
    add('AUPR',  average_precision_score(y,s), boot_ci(y,s,average_precision_score))
    add('F1',    f1_score(y,pred), boot_ci(y,s,lambda yy,ss:f1_score(yy,(ss>=0.5).astype(int))))
    add('MCC',   matthews_corrcoef(y,pred), boot_ci(y,s,lambda yy,ss:matthews_corrcoef(yy,(ss>=0.5).astype(int))))
    add('BalAcc',balanced_accuracy_score(y,pred), boot_ci(y,s,lambda yy,ss:balanced_accuracy_score(yy,(ss>=0.5).astype(int))))
    row['Acc']=round(accuracy_score(y,pred),3)
    return row

P('='*70); P('Phase TS-C — time-split 성능 평가'); P('='*70)
rows=[metrics(newp,'prospective(STAN-unseen)'), metrics(seenp,'leaked(STAN-seen)')]
# AR 서브셋
ar=newp[newp['Uniprot']=='P10275']
if ar['label'].nunique()>1: rows.append(metrics(ar,'prospective-AR(P10275)'))
M=pd.DataFrame(rows); M.to_csv(os.path.join(OUT,'timesplit_metrics.csv'),index=False)
P('\n[지표 표]'); P(M.to_string(index=False))

# Spearman (확률 vs -DC50, vs Dmax) — 수치형만
def sp(df,name):
    if 'dc50_nM' not in df.columns: return
    d=df.dropna(subset=['dc50_nM'])
    if len(d)>=5:
        r1,_=spearmanr(d['stan_active_prob'], -d['dc50_nM'])
        dd=df.dropna(subset=['dmax_pct'])
        r2,_=spearmanr(dd['stan_active_prob'], dd['dmax_pct'])
        P(f'[Spearman/{name}] prob vs -DC50 rho={r1:.3f} (n={len(d)}); prob vs Dmax rho={r2:.3f} (n={len(dd)})')
sp(newp,'prospective'); sp(seenp,'leaked')

# 타겟별 AUROC (양/음 각 >=5)
pt=[]
for u,g in newp.groupby('Uniprot'):
    y=g['label'].astype(int).values
    if (y.sum()>=5) and ((y==0).sum()>=5):
        pt.append({'Uniprot':u,'Target':g['Target'].iloc[0],'n':len(g),'pos':int(y.sum()),
                   'AUROC':round(roc_auc_score(y,g['stan_active_prob']),3)})
PT=pd.DataFrame(pt).sort_values('n',ascending=False)
PT.to_csv(os.path.join(OUT,'timesplit_per_target.csv'),index=False)
P(f'\n[타겟별 AUROC] (양·음 각≥5인 {len(PT)}개 타겟)'); P(PT.to_string(index=False))
P(f'  타겟별 AUROC: 중앙값={PT["AUROC"].median():.3f}, 평균={PT["AUROC"].mean():.3f}, '
  f'>0.5인 타겟={int((PT["AUROC"]>0.5).sum())}/{len(PT)}')

# ===== 그림 (영문 라벨; 폰트 안전) =====
def roc_fig():
    plt.figure(figsize=(4.2,4))
    for df,lab,c in [(newp,'prospective (STAN-unseen)','#1f77b4'),(seenp,'leaked (STAN-seen)','#d62728')]:
        y=df['label'].astype(int); fpr,tpr,_=roc_curve(y,df['stan_active_prob'])
        plt.plot(fpr,tpr,c,label=f'{lab}\nAUROC={roc_auc_score(y,df["stan_active_prob"]):.3f} (n={len(y)})')
    plt.plot([0,1],[0,1],'k--',lw=.8); plt.xlabel('FPR'); plt.ylabel('TPR')
    plt.title('PROTAC-STAN time-split ROC (all targets)'); plt.legend(fontsize=7,loc='lower right')
    plt.tight_layout(); plt.savefig(os.path.join(OUT,'ts_roc.png'),dpi=150); plt.close()
def pr_fig():
    plt.figure(figsize=(4.2,4))
    for df,lab,c in [(newp,'prospective','#1f77b4'),(seenp,'leaked','#d62728')]:
        y=df['label'].astype(int); pr,rc,_=precision_recall_curve(y,df['stan_active_prob'])
        plt.plot(rc,pr,c,label=f'{lab} AUPR={average_precision_score(y,df["stan_active_prob"]):.3f}')
        plt.axhline(y.mean(),color=c,ls=':',lw=.7)
    plt.xlabel('Recall'); plt.ylabel('Precision'); plt.title('Precision-Recall (time-split)')
    plt.legend(fontsize=7); plt.tight_layout(); plt.savefig(os.path.join(OUT,'ts_pr.png'),dpi=150); plt.close()
def dist_fig():
    plt.figure(figsize=(4.6,3.6))
    for lab,c in [(0,'#888'),(1,'#2ca02c')]:
        s=newp[newp['label']==lab]['stan_active_prob']
        plt.hist(s,bins=25,alpha=.6,color=c,label=f'label={lab} (n={len(s)})',density=True)
    plt.axvline(.5,color='k',ls='--',lw=.8); plt.xlabel('STAN active probability'); plt.ylabel('density')
    plt.title('Prospective: prob by true label'); plt.legend(fontsize=8)
    plt.tight_layout(); plt.savefig(os.path.join(OUT,'ts_probdist.png'),dpi=150); plt.close()
def pertarget_fig():
    d=PT.sort_values('AUROC')
    plt.figure(figsize=(5,max(3,0.32*len(d))))
    plt.barh(d['Target'],d['AUROC'],color=np.where(d['AUROC']>=0.5,'#1f77b4','#d62728'))
    plt.axvline(0.5,color='k',ls='--',lw=.8); plt.xlabel('AUROC'); plt.xlim(0,1)
    plt.title('Per-target AUROC (prospective, pos&neg>=5)'); plt.tight_layout()
    plt.savefig(os.path.join(OUT,'ts_per_target.png'),dpi=150); plt.close()
roc_fig(); pr_fig(); dist_fig(); pertarget_fig()
P('\n[그림] ts_roc.png, ts_pr.png, ts_probdist.png, ts_per_target.png')
open(os.path.join(OUT,'timesplit_eval_log.txt'),'w').write('\n'.join(log)+'\n')
P('[저장] outputs/timesplit_eval_log.txt')
