"""
case study 29개 전부 Boltz-2 ternary 입력(YAML) 생성.
29개 모두 AR-LBD + CRBN + 화합물 → 단백질 MSA는 B3 run 것 재사용(서버 재조회 없음).
출력: ~/PROTAC_MTL_v5/boltz_all29/inputs/<ID>.yaml
"""
import json, os, pandas as pd
from rdkit import Chem, RDLogger; RDLogger.DisableLog('rdApp.*')

AP='/home/kimjisan95/ar_protac_project'; MT='/home/kimjisan95/PROTAC_MTL_v5'
seqs=json.load(open(f'{MT}/boltz_B3_seqs.json'))
ar=seqs['ar_lbd']; crbn=seqs['crbn']
MSA_A=f'{MT}/boltz_B3_ternary/out/boltz_results_B3_ternary/msa/B3_ternary_0.csv'  # AR-LBD
MSA_B=f'{MT}/boltz_B3_ternary/out/boltz_results_B3_ternary/msa/B3_ternary_1.csv'  # CRBN
assert os.path.exists(MSA_A) and os.path.exists(MSA_B), 'B3 MSA 없음'

cs=pd.read_csv(f'{AP}/data/processed/case_study_clean.csv')
outdir=f'{MT}/boltz_all29/inputs'; os.makedirs(outdir, exist_ok=True)
n=0
for _,r in cs.iterrows():
    cid=str(r['ID']); smi=str(r['Smiles'])
    if Chem.MolFromSmiles(smi) is None:
        print('SKIP(파싱실패):', cid); continue
    yaml=f"""version: 1
sequences:
  - protein:
      id: A
      sequence: {ar}
      msa: {MSA_A}
  - protein:
      id: B
      sequence: {crbn}
      msa: {MSA_B}
  - ligand:
      id: C
      smiles: '{smi}'
"""
    open(f'{outdir}/{cid}.yaml','w').write(yaml); n+=1
print(f'[생성] {n}개 YAML → {outdir} (MSA 재사용: AR-LBD/CRBN B3 것)')
