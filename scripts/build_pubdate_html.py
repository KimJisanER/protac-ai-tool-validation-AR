"""
게재연도 분포 시각화: AI검증 hold-out(prospective) vs 기존 PROTAC-DB.
DOI→연도(Crossref, data/processed/doi_years.csv) → 화합물(ik14)별 최초 게재연도 →
held vs rest 분포를 인터랙티브 Plotly HTML로.
출력: pubdate_distribution.html (+ 콘솔 요약)
"""
import json, html, numpy as np, pandas as pd
from rdkit import Chem, RDLogger; RDLogger.DisableLog('rdApp.*')
ROOT='/home/kimjisan95/ar_protac_project'
def ik(s):
    m=Chem.MolFromSmiles(str(s)); return Chem.MolToInchiKey(m)[:14] if m else None

x=pd.read_excel(f'{ROOT}/protac.xlsx'); x['ik14']=x['Smiles'].map(ik)
x['doi']=x['Article DOI'].astype(str).str.strip()
yr=pd.read_csv(f'{ROOT}/data/processed/doi_years.csv'); yr['doi']=yr['doi'].astype(str).str.strip()
ymap=dict(zip(yr['doi'],yr['year']))
x['year']=x['doi'].map(ymap)
held=set(pd.read_csv(f'{ROOT}/data/processed/timesplit_compounds.csv')['ik14'])

# 화합물(ik14)별 최초 게재연도 + 그룹
comp=(x.dropna(subset=['ik14','year']).groupby('ik14')
        .agg(year=('year','min')).reset_index())
comp['year']=comp['year'].astype(int)
comp['grp']=comp['ik14'].map(lambda i:'held' if i in held else 'rest')
comp=comp[(comp['year']>=2000)&(comp['year']<=2026)]
held_y=comp[comp['grp']=='held']['year']; rest_y=comp[comp['grp']=='rest']['year']
print(f"연도 확보 화합물: held {len(held_y)} / rest {len(rest_y)}")
print(f"  held 중앙값 {int(held_y.median())}, rest 중앙값 {int(rest_y.median())}")
print(f"  held {held_y.quantile(.25):.0f}~{held_y.quantile(.75):.0f} / rest {rest_y.quantile(.25):.0f}~{rest_y.quantile(.75):.0f}")

years=list(range(int(comp['year'].min()),2027))
def counts(s):
    vc=s.value_counts(); return [int(vc.get(y,0)) for y in years]
hc=counts(held_y); rc=counts(rest_y)
payload={'years':years,'held':hc,'rest':rc,
         'held_n':int(len(held_y)),'rest_n':int(len(rest_y)),
         'held_med':int(held_y.median()),'rest_med':int(rest_y.median())}
js=json.dumps(payload)
page=f"""<!DOCTYPE html><html lang="ko"><head><meta charset="utf-8">
<title>게재연도 분포 — AI검증 hold-out vs 기존 PROTAC-DB</title>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<style>body{{margin:0;font-family:-apple-system,Segoe UI,'Noto Sans KR',sans-serif;background:#f8fafc}}
 header{{padding:14px 20px;background:#fff;border-bottom:1px solid #e2e8f0}} h1{{margin:0;font-size:18px}}
 .sub{{color:#64748b;font-size:13px;margin-top:4px}} #c{{height:46vh}} #c2{{height:40vh}}
 .row{{padding:8px 20px;color:#334155;font-size:13px}} b{{color:#0f172a}}</style></head><body>
<header><h1>게재연도 분포 — AI 검증 hold-out vs 기존 PROTAC-DB</h1>
<div class="sub">DOI→연도(Crossref) · 화합물(ik14)별 최초 게재연도 · 주황=검증 hold-out(prospective) / 회색=나머지(기존)</div></header>
<div class="row">hold-out 중앙값 <b>{payload['held_med']}</b> (n={payload['held_n']}) · 기존 DB 중앙값 <b>{payload['rest_med']}</b> (n={payload['rest_n']}) — 위=정규화(분포 형태 비교), 아래=실제 개수</div>
<div id="c"></div><div id="c2"></div>
<script>
const D={js};
function pct(a){{const s=a.reduce((x,y)=>x+y,0)||1;return a.map(v=>100*v/s);}}
const norm=[
 {{x:D.years,y:pct(D.rest),type:'bar',name:'기존 PROTAC-DB ('+D.rest_n+')',marker:{{color:'rgba(148,163,184,0.75)'}}}},
 {{x:D.years,y:pct(D.held),type:'bar',name:'검증 hold-out ('+D.held_n+')',marker:{{color:'rgba(249,115,22,0.85)'}}}}];
Plotly.newPlot('c',norm,{{barmode:'overlay',bargap:0.05,margin:{{l:50,r:10,t:10,b:30}},
 yaxis:{{title:'그룹 내 %'}},xaxis:{{title:'게재연도',dtick:1}},legend:{{x:0.01,y:0.99}},
 title:''}},{{responsive:true,displaylogo:false}});
const raw=[
 {{x:D.years,y:D.rest,type:'bar',name:'기존 PROTAC-DB',marker:{{color:'rgba(148,163,184,0.75)'}}}},
 {{x:D.years,y:D.held,type:'bar',name:'검증 hold-out',marker:{{color:'rgba(249,115,22,0.85)'}}}}];
Plotly.newPlot('c2',raw,{{barmode:'group',margin:{{l:50,r:10,t:10,b:35}},
 yaxis:{{title:'화합물 수'}},xaxis:{{title:'게재연도',dtick:1}},legend:{{x:0.01,y:0.99}}}},{{responsive:true,displaylogo:false}});
</script></body></html>"""
open(f'{ROOT}/pubdate_distribution.html','w').write(page)
print('[저장] pubdate_distribution.html')
