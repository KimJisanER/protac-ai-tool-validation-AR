"""
게재연도 분포 KDE (PNG): AI검증 hold-out(prospective) vs 기존 PROTAC-DB.
DOI→연도(data/processed/doi_years.csv) → 화합물(ik14)별 최초 게재연도 → 그룹별 KDE.
출력: outputs/pubdate_kde.png
"""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from rdkit import Chem, RDLogger; RDLogger.DisableLog('rdApp.*')
ROOT='/home/kimjisan95/ar_protac_project'
def ik(s):
    m=Chem.MolFromSmiles(str(s)); return Chem.MolToInchiKey(m)[:14] if m else None

x=pd.read_excel(f'{ROOT}/protac.xlsx'); x['ik14']=x['Smiles'].map(ik)
yr=pd.read_csv(f'{ROOT}/data/processed/doi_years.csv'); yr['doi']=yr['doi'].astype(str).str.strip()
x['year']=x['Article DOI'].astype(str).str.strip().map(dict(zip(yr['doi'],yr['year'])))
held=set(pd.read_csv(f'{ROOT}/data/processed/timesplit_compounds.csv')['ik14'])
comp=x.dropna(subset=['ik14','year']).groupby('ik14').agg(year=('year','min')).reset_index()
comp['year']=comp['year'].astype(int); comp=comp[(comp['year']>=2000)&(comp['year']<=2026)]
comp['grp']=comp['ik14'].map(lambda i:'held' if i in held else 'rest')
h=comp[comp['grp']=='held']['year'].values; r=comp[comp['grp']=='rest']['year'].values
print(f"held n={len(h)} median={int(np.median(h))} | rest n={len(r)} median={int(np.median(r))}")

XMIN,XMAX=2017,2026
grid=np.linspace(2000,2027,600)
fig,ax=plt.subplots(figsize=(8,4.6))
recent=[2023,2024,2025,2026]
for arr,color,lab in [(r,'#94a3b8','Existing PROTAC-DB'),(h,'#f97316','AI-validation hold-out (prospective)')]:
    kde=gaussian_kde(arr.astype(float)); d=kde(grid)
    ax.fill_between(grid,d,color=color,alpha=0.30)
    ax.plot(grid,d,color=color,lw=2.2,label=f'{lab}  (n={len(arr)}, median={int(np.median(arr))})')
    ax.axvline(np.median(arr),color=color,ls='--',lw=1.2)
    # 최근 연도별 개수 점 + 라벨(KDE 곡선 위)
    for y in recent:
        cnt=int((arr==y).sum())
        if cnt==0: continue
        dy=kde(float(y))[0]
        ax.scatter([y],[dy],s=22,color=color,zorder=5,edgecolor='white',linewidth=0.6)
        ax.annotate(str(cnt),(y,dy),textcoords='offset points',xytext=(0,7),
                    ha='center',fontsize=8,fontweight='bold',color=color)
ax.set_xlabel('Publication year (first appearance per compound)'); ax.set_ylabel('Density (KDE)')
ax.set_title('Publication-year distribution: hold-out vs existing PROTAC-DB\n(zoom 2017–2026; numbers = compound count that year)')
ax.set_xlim(XMIN,XMAX); ax.set_xticks(range(XMIN,XMAX+1)); ax.tick_params(axis='x',labelsize=9)
ax.grid(axis='x',ls=':',alpha=0.4)
# 2023~2026 강조 음영
ax.axvspan(2022.5,2026.5,color='#fde68a',alpha=0.18,zorder=0)
ax.legend(fontsize=9,loc='upper left'); fig.tight_layout()
fig.savefig(f'{ROOT}/outputs/pubdate_kde.png',dpi=160); plt.close()
# 콘솔: 최근 연도별 개수
print('연도별 개수(held / rest):')
for y in recent: print(f'  {y}: {int((h==y).sum())} / {int((r==y).sum())}')
print('[저장] outputs/pubdate_kde.png')
