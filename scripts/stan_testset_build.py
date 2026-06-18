"""
PROTAC-STAN 자체 train/test split을 STAN 입력 포맷(ts_new 동일)으로 빌드.
- test (207행): STAN이 학습 때 보지 않은 자체 held-out 벤치마크 = in-distribution 일반화
- train (1296행): STAN이 학습한 화합물 = 암기(memorization) 상한
출력: ~/PROTAC-STAN/data/{stan_test,stan_train}/<sub>.csv + esm_s_map.pkl
"""
import os, shutil, pandas as pd
from rdkit import Chem, RDLogger; RDLogger.DisableLog('rdApp.*')
def ik14(s):
    m=Chem.MolFromSmiles(str(s)); return Chem.MolToInchiKey(m)[:14] if m else None

STAN='/home/kimjisan95/PROTAC-STAN'; FINE=f'{STAN}/data/PROTAC-fine'
COLS=['Smiles','Molecular Weight','Exact Mass','XLogP3','Heavy Atom Count','Ring Count',
      'Hydrogen Bond Acceptor Count','Hydrogen Bond Donor Count','Rotatable Bond Count',
      'Topological Polar Surface Area','Uniprot','E3 ligase Uniprot','label','ik14']
fine=pd.read_csv(f'{FINE}/protac-fine.csv'); fine['cid']=fine['Compound ID'].astype(str)
te=set(pd.read_csv(f'{FINE}/test_compound_smiles.csv')['Compound ID'].astype(str))
fine['ik14']=fine['Smiles'].map(ik14)
fine['label']=fine['label'].astype(int)
fine['grp']=fine['cid'].map(lambda c:'test' if c in te else 'train')

for grp,sub in [('test','stan_test'),('train','stan_train')]:
    d=f'{STAN}/data/{sub}'; os.makedirs(d,exist_ok=True)
    g=fine[fine.grp==grp][COLS].copy()
    g.to_csv(f'{d}/{sub}.csv',index=False)
    shutil.copy(f'{FINE}/esm_s_map.pkl', f'{d}/esm_s_map.pkl')
    print(f'[{sub}] {len(g)}행 (양성 {int(g.label.sum())}/음성 {int((g.label==0).sum())}) -> {d}/{sub}.csv')
