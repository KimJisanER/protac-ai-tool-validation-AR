"""
Boltz-2 cofold 29개 삼원복합체를 CRBN(chain B) 기준으로 정렬 → AR-LBD(chain A)의 상대 위치로
클러스터링(= 두 단백질의 상대 배치 군집) → 클러스터별 DC50/Dmax 차이 검정.
정렬 PDB와 PyMOL 스크립트도 출력(사용자가 PyMOL로 확인).
env: protac (biopython 1.86, scipy, sklearn)
출력: outputs/cofold_cluster_assignments.csv, outputs/cofold_cluster_summary.csv,
      outputs/fig_cofold_clusters.png, outputs/cofold_aligned/<ID>.pdb, outputs/view_cofold_clusters.pml
"""
import os, glob, numpy as np, pandas as pd
from Bio.PDB import PDBParser, Superimposer, PDBIO
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from scipy.spatial.distance import squareform
from sklearn.metrics import silhouette_score
from scipy.stats import kruskal

PRED='/home/kimjisan95/PROTAC_MTL_v5/boltz_all29/out/boltz_results_inputs/predictions'
AP='/home/kimjisan95/ar_protac_project'; OUT=f'{AP}/outputs'
ALN=f'{OUT}/cofold_aligned'; os.makedirs(ALN, exist_ok=True)
log=[]; P=lambda *a:(print(*a), log.append(' '.join(str(x) for x in a)))
parser=PDBParser(QUIET=True)

ids=sorted([os.path.basename(d) for d in glob.glob(f'{PRED}/*') if os.path.isdir(d)])
def ca(struct, ch):
    return [r['CA'] for r in struct[0][ch] if 'CA' in r]

# 기준 = B3(리드). CRBN(B)로 모든 구조를 동일 프레임에 정렬, AR-LBD(A) CA 좌표 추출
ref=parser.get_structure('B3', f'{PRED}/B3/B3_model_0.pdb'); ref_crbn=ca(ref,'B')
ar_coords={}; structs={}
for cid in ids:
    s=parser.get_structure(cid, f'{PRED}/{cid}/{cid}_model_0.pdb')
    sup=Superimposer(); sup.set_atoms(ref_crbn, ca(s,'B')); sup.apply(s.get_atoms())  # CRBN 정렬→전체 변환
    ar_coords[cid]=np.array([a.coord for a in ca(s,'A')]); structs[cid]=s
    io=PDBIO(); io.set_structure(s); io.save(f'{ALN}/{cid}.pdb')           # 정렬 PDB 저장(PyMOL용)
n=len(ids); P(f'{n}개 정렬 완료(CRBN 기준). AR-LBD CA={ar_coords[ids[0]].shape[0]}개')

# 쌍별 AR-LBD CA RMSD (CRBN 정렬 후, 추가정렬 없음 = 상대배치 차이)
R=np.zeros((n,n))
for i in range(n):
    for j in range(i+1,n):
        d=np.sqrt(((ar_coords[ids[i]]-ar_coords[ids[j]])**2).sum(1).mean())
        R[i,j]=R[j,i]=d
P(f'AR-LBD 상대 RMSD: 평균 {R[np.triu_indices(n,1)].mean():.1f}Å, 최대 {R.max():.1f}Å')

# 계층 클러스터링 + K 선택(silhouette)
Z=linkage(squareform(R), method='average')
best=None
for K in range(2,7):
    lab=fcluster(Z, K, criterion='maxclust')
    if len(set(lab))<2: continue
    sil=silhouette_score(R, lab, metric='precomputed')
    P(f'  K={K}: silhouette={sil:.3f}, 크기={np.bincount(lab)[1:].tolist()}')
    if best is None or sil>best[1]: best=(K,sil,lab)
K,sil,lab=best; P(f'[선택] K={K} (silhouette={sil:.3f})')

df=pd.DataFrame({'ID':ids,'cluster':lab})
cs=pd.read_csv(f'{AP}/data/processed/case_study_clean.csv')[['ID','series','dc50_nM','dc50_censored','dmax_pct']]
df=df.merge(cs,on='ID',how='left').sort_values(['cluster','ID'])
df.to_csv(f'{OUT}/cofold_cluster_assignments.csv',index=False)
P('\n[클러스터별 구성]')
for c,g in df.groupby('cluster'):
    P(f'  C{c} (n={len(g)}): {g["series"].value_counts().to_dict()} | DC50보유 {g["dmax_pct"].notna().sum()}')

# DC50/Dmax 클러스터간 차이 (값 보유 13개만)
val=df.dropna(subset=['dmax_pct']).copy()
P(f'\n[DC50/Dmax 클러스터간 차이 — 값 보유 {len(val)}개]')
summ=[]
for c,g in val.groupby('cluster'):
    dc=g['dc50_nM']; dcmed = dc[dc<1000].median() if (dc<1000).any() else float('nan')
    P(f'  C{c} (n={len(g)}): Dmax 중앙 {g.dmax_pct.median():.0f}% (범위 {g.dmax_pct.min():.0f}-{g.dmax_pct.max():.0f}) | DC50 수치형 {(dc<1000).sum()}개 중앙 {dcmed:.1f}nM')
    summ.append({'cluster':c,'n_total':int((df.cluster==c).sum()),'n_valued':len(g),
                 'dmax_median':round(g.dmax_pct.median(),1),'dc50_numeric_median':round(dc[dc<1000].median(),1) if (dc<1000).any() else None})
pd.DataFrame(summ).to_csv(f'{OUT}/cofold_cluster_summary.csv',index=False)
cls=[g['dmax_pct'].values for _,g in val.groupby('cluster') if len(g)>=2]
if len(cls)>=2:
    H,p=kruskal(*cls); P(f'  Kruskal-Wallis(Dmax ~ cluster): H={H:.2f}, p={p:.3f} (n_clusters={len(cls)})')
    dcl=[g['dc50_nM'].values for _,g in val.groupby('cluster') if len(g)>=2]
    H2,p2=kruskal(*dcl); P(f'  Kruskal-Wallis(DC50[censored=1000] ~ cluster): H={H2:.2f}, p={p2:.3f}')

# 전역 K=2가 chemotype(A vs B/C)·데이터유무와 겹치므로, 값 보유 13개만으로 재클러스터링해 실제 검정
P('\n[값 보유 13개(B/C) 한정 재클러스터링 — 구조 vs DC50/Dmax]')
vids=val['ID'].tolist(); idx=[ids.index(v) for v in vids]
Rv=R[np.ix_(idx,idx)]
P(f'  13개 AR-LBD 상대 RMSD: 평균 {Rv[np.triu_indices(len(idx),1)].mean():.1f}Å (최대 {Rv.max():.1f}Å)')
Zv=linkage(squareform(Rv), method='average')
for Kv in [2,3]:
    lv=fcluster(Zv, Kv, criterion='maxclust'); v2=val.copy(); v2['subcl']=lv
    P(f'  K={Kv}: 크기={np.bincount(lv)[1:].tolist()}')
    for c,g in v2.groupby('subcl'):
        nnum=int((g.dc50_nM<1000).sum())
        P(f'    sub{c} (n={len(g)}): Dmax 중앙 {g.dmax_pct.median():.0f}% | DC50 수치형 {nnum}개 | 멤버 {g.ID.tolist()}')
    grp=[g['dmax_pct'].values for _,g in v2.groupby('subcl') if len(g)>=2]
    if len(grp)>=2:
        H,p=kruskal(*grp); P(f'    KW(Dmax~subcluster): H={H:.2f}, p={p:.3f}')
        grpd=[g['dc50_nM'].values for _,g in v2.groupby('subcl') if len(g)>=2]
        H3,p3=kruskal(*grpd); P(f'    KW(DC50[cens=1000]~subcluster): H={H3:.2f}, p={p3:.3f}')

# 그림: 덴드로그램 + RMSD 히트맵 + 클러스터별 Dmax
fig,ax=plt.subplots(1,3,figsize=(15,4.6))
dendrogram(Z, labels=ids, ax=ax[0], color_threshold=Z[-(K-1),2]); ax[0].set_title(f'AR-LBD relative-pose dendrogram (K={K})',fontsize=10); ax[0].tick_params(axis='x',labelsize=6)
order=dendrogram(Z,no_plot=True)['leaves']
im=ax[1].imshow(R[np.ix_(order,order)],cmap='viridis'); ax[1].set_title('pairwise AR-LBD RMSD (Å)',fontsize=10); plt.colorbar(im,ax=ax[1],fraction=.046)
ax[1].set_xticks([]); ax[1].set_yticks([])
for c,g in val.groupby('cluster'):
    ax[2].scatter([c]*len(g), g['dmax_pct'], s=40, label=f'C{c}')
    for _,r in g.iterrows(): ax[2].annotate(r['ID'],(r['cluster'],r['dmax_pct']),fontsize=6,xytext=(3,0),textcoords='offset points')
ax[2].set_xlabel('cluster'); ax[2].set_ylabel('measured Dmax (%)'); ax[2].set_title('Dmax by structural cluster',fontsize=10)
plt.tight_layout(); plt.savefig(f'{OUT}/fig_cofold_clusters.png',dpi=150); plt.close()

# PyMOL 스크립트(사용자 실행용): 정렬 PDB 로드, 클러스터별 AR-LBD 색, CRBN 회색
colors=['marine','orange','green','magenta','yellow','cyan']
pml=[f'# PyMOL: cofold 클러스터 시각화 (CRBN 정렬). 실행: pymol {OUT}/view_cofold_clusters.pml','bg_color white']
for _,r in df.iterrows():
    pml.append(f'load {ALN}/{r.ID}.pdb, {r.ID}')
pml.append('hide everything')
pml.append('show cartoon')
# CRBN(chain B) 전부 회색, AR-LBD(chain A) 클러스터색, ligand(chain C) 스틱
pml.append('color gray70, chain B')
for _,r in df.iterrows():
    pml.append(f'color {colors[(int(r.cluster)-1)%len(colors)]}, {r.ID} and chain A')
pml.append('show sticks, chain C')
pml.append('set cartoon_transparency, 0.1')
pml.append('orient')
pml.append(f'png {OUT}/cofold_clusters_overview.png, width=1600, height=1200, dpi=150, ray=1')
open(f'{OUT}/view_cofold_clusters.pml','w').write('\n'.join(pml)+'\n')
P(f'\n[저장] 정렬PDB {ALN}/, 클러스터 CSV, fig_cofold_clusters.png, view_cofold_clusters.pml')
open(f'{OUT}/cofold_cluster_log.txt','w').write('\n'.join(log)+'\n')
