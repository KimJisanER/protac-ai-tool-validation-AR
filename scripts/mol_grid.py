"""
case study 29개 분자를 RDKit MCS로 같은 방향 정렬해 2D 구조 그리드(보고서용).
PDB 그리드와 같은 순서(B→C→A). 흰 배경·검은 글씨, ID 좌상단(검정).
각 분자 하단: DC50/Dmax 실측 + good degrader(STAN기준 DC50<100&Dmax≥80) 실측 여부 + STAN 예측 여부.
출력: outputs/mol_grid_BCA.png
"""
import numpy as np, pandas as pd
from rdkit import Chem, RDLogger; RDLogger.DisableLog('rdApp.*')
from rdkit.Chem import rdFMCS, rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import io
from PIL import Image

cs=pd.read_csv('data/processed/case_study_clean.csv')
t2=pd.read_csv('outputs/table2.csv')[['ID','stan_active_prob','stan_pred']]
d=cs.merge(t2,on='ID',how='left')
order=[f'B{i}' for i in range(1,7)]+[f'C{i}' for i in range(1,8)]+[f'A{i}' for i in range(1,17)]
d=d.set_index('ID').loc[order].reset_index()

mols=[Chem.MolFromSmiles(s) for s in d['Smiles']]
# MCS로 공통 코어 → 같은 방향 정렬
res=rdFMCS.FindMCS(mols, ringMatchesRingOnly=True, completeRingsOnly=True, timeout=40)
core=Chem.MolFromSmarts(res.smartsString); rdDepictor.Compute2DCoords(core)
print(f'MCS: {res.numAtoms} atoms / {res.numBonds} bonds')
for m in mols:
    try:
        if core is not None and m.HasSubstructMatch(core):
            rdDepictor.GenerateDepictionMatching2DStructure(m, core)
        else: rdDepictor.Compute2DCoords(m)
    except Exception: rdDepictor.Compute2DCoords(m)

def mol_png(m):
    dr=rdMolDraw2D.MolDraw2DCairo(560,420); dr.drawOptions().bondLineWidth=2
    rdMolDraw2D.PrepareAndDrawMolecule(dr,m); dr.FinishDrawing()
    return np.array(Image.open(io.BytesIO(dr.GetDrawingText())))

def ann(r):
    if pd.isna(r.dc50_nM):
        l1='DC50/Dmax: not measured'; l2='exp good-degrader: n/a'
    else:
        dc=('>1000' if str(r.dc50_censored)=='>' else f'{r.dc50_nM:.1f}')
        l1=f'DC50 {dc} nM | Dmax {r.dmax_pct:.0f}%'
        gd='Yes' if (r.dc50_nM<100 and r.dmax_pct>=80) else 'No'
        l2=f'exp good-degrader: {gd}'
    pred='Yes' if r.stan_pred==1 else 'No'
    l3=f'STAN pred: {pred} (p={r.stan_active_prob:.2f})'
    return l1,l2,l3

ncol,nrow=5,6
fig,axes=plt.subplots(nrow,ncol,figsize=(ncol*3.6,nrow*4.2)); fig.patch.set_facecolor('white')
for k,ax in enumerate(axes.flat):
    ax.set_xlim(0,1); ax.set_ylim(0,1); ax.axis('off')
    if k>=len(order): continue
    r=d.iloc[k]
    ax.imshow(mol_png(mols[k]), extent=(0.0,1.0,0.30,1.0), aspect='auto')
    ax.text(0.02,0.99,r.ID,ha='left',va='top',fontsize=20,fontweight='bold',color='black')
    l1,l2,l3=ann(r)
    ax.text(0.5,0.235,l1,ha='center',va='top',fontsize=14,color='black')
    ax.text(0.5,0.150,l2,ha='center',va='top',fontsize=14,color='black')
    ax.text(0.5,0.065,l3,ha='center',va='top',fontsize=14,color='black')
fig.suptitle('AR-PROTAC case study (B→C→A): 2D structure (MCS-aligned) + measured DC50/Dmax, good-degrader (DC50<100nM & Dmax≥80%), and PROTAC-STAN prediction',
             fontsize=15, y=0.997)
plt.subplots_adjust(left=0.005,right=0.995,top=0.975,bottom=0.005,wspace=0.03,hspace=0.05)
plt.savefig('outputs/mol_grid_BCA.png',dpi=140,facecolor='white'); print('[저장] outputs/mol_grid_BCA.png')
