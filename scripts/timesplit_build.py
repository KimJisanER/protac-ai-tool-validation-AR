"""
Phase TS-A — 전체 타겟 time-split prospective 평가 세트 구축.

목적: AR(P10275) 한정이 아니라, 최신 PROTAC-DB 3.0(protac.xlsx)에 수록된
전체 타겟 PROTAC 중 'PROTAC-STAN이 학습 때 보지 못한' 화합물을 prospective
test로 분리하여 도구(STAN)의 일반 분해예측 성능을 검증한다.

산출:
  - data/processed/timesplit_compounds.csv      (prospective, compound-level dedup)
  - data/processed/timesplit_leaked.csv         (대조: STAN 학습셋에 있는 화합물)
  - ~/PROTAC-STAN/data/ts_new/ts_new.csv  + esm_s_map.pkl   (STAN 입력)
  - ~/PROTAC-STAN/data/ts_seen/ts_seen.csv + esm_s_map.pkl  (STAN 입력, 누수 대조)
  - outputs/timesplit_log.txt 에 카운트 동결
"""
import os, re, shutil, pickle
import numpy as np, pandas as pd
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')

ROOT = '/home/kimjisan95/ar_protac_project'
STAN = os.path.expanduser('~/PROTAC-STAN')
FINE = os.path.join(STAN, 'data/PROTAC-fine')

PROP_COLS = ['Molecular Weight','Exact Mass','XLogP3','Heavy Atom Count','Ring Count',
             'Hydrogen Bond Acceptor Count','Hydrogen Bond Donor Count',
             'Rotatable Bond Count','Topological Polar Surface Area']

def num(v):
    if pd.isna(v): return np.nan
    s = str(v).replace(',', '')
    m = re.search(r'[-+]?\d*\.?\d+', s)
    return float(m.group()) if m else np.nan

def is_censored_gt(v):
    """'>1000' 류 = 상한 검열(=비활성 방향)."""
    return bool(re.search(r'>', str(v))) if pd.notna(v) else False

def ik14(s):
    m = Chem.MolFromSmiles(str(s))
    return Chem.MolToInchiKey(m)[:14] if m else None

def rdkit_props(smi):
    from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors as rd
    m = Chem.MolFromSmiles(str(smi))
    if m is None: return None
    return {
        'Molecular Weight': Descriptors.MolWt(m),
        'Exact Mass': Descriptors.ExactMolWt(m),
        'XLogP3': Crippen.MolLogP(m),          # XLogP3 자리에 RDKit MolLogP (data.py columns 일치)
        'Heavy Atom Count': m.GetNumHeavyAtoms(),
        'Ring Count': rd.CalcNumRings(m),
        'Hydrogen Bond Acceptor Count': rd.CalcNumHBA(m),
        'Hydrogen Bond Donor Count': rd.CalcNumHBD(m),
        'Rotatable Bond Count': rd.CalcNumRotatableBonds(m),
        'Topological Polar Surface Area': rd.CalcTPSA(m),
    }

log = []
def P(*a):
    s = ' '.join(str(x) for x in a); print(s); log.append(s)

P('='*70); P('Phase TS-A — 전체 타겟 time-split 세트 구축'); P('='*70)

# --- ESM 커버리지 & E3 name->Uniprot 권위 매핑 (STAN 자체 자산에서) ---
esm = set(pickle.load(open(os.path.join(FINE,'esm_s_map.pkl'),'rb')).keys())
P(f'[ESM] 커버 Uniprot 수 = {len(esm)} (타겟+E3 합산)')

fine = pd.read_csv(os.path.join(FINE,'protac-fine.csv'))
e3map = (fine.dropna(subset=['E3 ligase','E3 ligase Uniprot'])
            .groupby('E3 ligase')['E3 ligase Uniprot']
            .agg(lambda s: s.value_counts().index[0]).to_dict())
P(f'[E3 매핑] {e3map}')

# STAN이 학습/관측한 모든 화합물 골격키(누수 기준)
fine['ik'] = fine['Smiles'].map(ik14)
seen_ik = set(fine['ik'].dropna())
P(f'[누수기준] STAN protac-fine unique ik14 = {len(seen_ik)}')

# --- 신DB 3.0 라벨 후보 ---
x = pd.read_excel(os.path.join(ROOT,'protac.xlsx'))
P(f'[신DB] 전체 {len(x)}행 x {x.shape[1]}열')
x['dc'] = x['DC50 (nM)'].map(num)
x['dm'] = x['Dmax (%)'].map(num)
x['dc_censored'] = x['DC50 (nM)'].map(is_censored_gt)
x['e3u'] = x['E3 ligase'].map(e3map)
x['ik'] = x['Smiles'].map(ik14)

lab = x.dropna(subset=['dc','dm','ik']).copy()
P(f'[신DB] DC50+Dmax 둘 다 보유(파싱) = {len(lab)}행')
cov = lab[lab['Uniprot'].isin(esm) & lab['e3u'].isin(esm)].copy()
P(f'[신DB] +타겟·E3 모두 ESM 커버 = {len(cov)}행 (타겟 {cov["Uniprot"].nunique()}종, E3 {cov["e3u"].nunique()}종)')

# 활성 라벨 정의 = STAN 자기 기준: DC50<100 AND Dmax>=80  (검열 '>'는 비활성)
def active_label(r):
    dc_ok = (not r['dc_censored']) and (r['dc'] < 100)
    return bool(dc_ok and (r['dm'] >= 80))
cov['active'] = cov.apply(active_label, axis=1).astype(int)

# compound-level dedup: ik14 그룹별 median dc/dm, label은 재계산, 대표 SMILES=첫 행
def agg_group(g):
    dc = g['dc'].median(); dm = g['dm'].median()
    cen = bool((g['dc_censored'] & (g['dc']>=1000)).all()) if g['dc_censored'].any() else False
    act = int((not cen) and (dc < 100) and (dm >= 80))
    row = g.iloc[0]
    return pd.Series({'Smiles': row['Smiles'], 'Uniprot': row['Uniprot'], 'Target': row['Target'],
                      'E3 ligase': row['E3 ligase'], 'e3u': row['e3u'],
                      'dc50_nM': dc, 'dmax_pct': dm, 'n_assay': len(g), 'active': act})
comp = cov.groupby('ik', as_index=False).apply(agg_group, include_groups=False)
comp = comp.rename(columns={'ik':'ik14'}) if 'ik' in comp.columns else comp
comp.columns = ['ik14'] + list(comp.columns[1:]) if comp.columns[0]!='ik14' else comp.columns
P(f'[dedup] compound-level unique ik14 = {len(comp)}')

# prospective (STAN 미관측) vs leaked (관측)
comp['leaked'] = comp['ik14'].isin(seen_ik)
newc = comp[~comp['leaked']].copy()
seenc = comp[comp['leaked']].copy()
P(f'[분할] prospective(미관측) = {len(newc)}  /  leaked(관측) = {len(seenc)}')
P(f'[prospective 라벨] 양성 {int(newc["active"].sum())} / 음성 {int((~newc["active"].astype(bool)).sum())} '
  f'(양성률 {newc["active"].mean():.1%})')
P(f'[leaked 라벨]      양성 {int(seenc["active"].sum())} / 음성 {int((~seenc["active"].astype(bool)).sum())} '
  f'(양성률 {seenc["active"].mean():.1%})')
P(f'[prospective] 타겟 수 = {newc["Uniprot"].nunique()}, AR(P10275) 포함 화합물 = {int((newc["Uniprot"]=="P10275").sum())}')

# RDKit 9물성 부착 + STAN 입력 CSV 작성
def write_stan_input(df, sub):
    rows = []
    for _, r in df.iterrows():
        pr = rdkit_props(r['Smiles'])
        if pr is None: continue
        d = {'Smiles': r['Smiles'], 'Uniprot': r['Uniprot'], 'E3 ligase Uniprot': r['e3u'],
             'label': int(r['active']), 'ik14': r['ik14']}
        d.update(pr)
        rows.append(d)
    out = pd.DataFrame(rows)
    cols = ['Smiles'] + PROP_COLS + ['Uniprot','E3 ligase Uniprot','label','ik14']
    out = out[cols]
    dest = os.path.join(STAN, 'data', sub)
    os.makedirs(dest, exist_ok=True)
    # 캐시 제거
    pdir = os.path.join(dest, 'processed', sub)
    if os.path.isdir(pdir): shutil.rmtree(pdir)
    out.to_csv(os.path.join(dest, f'{sub}.csv'), index=False)
    shutil.copy(os.path.join(FINE,'esm_s_map.pkl'), os.path.join(dest,'esm_s_map.pkl'))
    P(f'[STAN입력] {dest}/{sub}.csv  ({len(out)}행)  + esm_s_map.pkl')
    return out

os.makedirs(os.path.join(ROOT,'data/processed'), exist_ok=True)
newc.to_csv(os.path.join(ROOT,'data/processed/timesplit_compounds.csv'), index=False)
seenc.to_csv(os.path.join(ROOT,'data/processed/timesplit_leaked.csv'), index=False)
write_stan_input(newc, 'ts_new')
write_stan_input(seenc, 'ts_seen')

# DOI 연도 분포(가능 범위; ACS 4c=2024 류 휴리스틱은 신뢰 낮음 → 구조적 분할이 주 근거)
P('[참고] 시간성 근거: 분할은 "STAN 학습 미관측" 구조적 기준이 주(主). '
  'protac.xlsx Article DOI 기반 엄밀 연도분할은 저널별 DOI 패턴 불균일로 보조만.')

open(os.path.join(ROOT,'outputs/timesplit_log.txt'),'w').write('\n'.join(log)+'\n')
P('[저장] outputs/timesplit_log.txt')
