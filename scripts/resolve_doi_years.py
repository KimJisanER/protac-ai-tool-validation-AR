"""
protac.xlsx의 unique Article DOI(766) → Crossref로 게재연도 조회 (캐시).
출력: data/processed/doi_years.csv  (doi, year)
"""
import os, json, time, urllib.request, urllib.parse, pandas as pd
ROOT='/home/kimjisan95/ar_protac_project'
OUT=f'{ROOT}/data/processed/doi_years.csv'
x=pd.read_excel(f'{ROOT}/protac.xlsx')
dois=sorted(x['Article DOI'].dropna().astype(str).str.strip().unique())
print("unique DOI:", len(dois))

cache={}
if os.path.exists(OUT):
    c=pd.read_csv(OUT); cache=dict(zip(c['doi'].astype(str),c['year']))
    print("캐시 로드:", len(cache))

def cr_year(doi):
    url='https://api.crossref.org/works/'+urllib.parse.quote(doi)+'?mailto=kimjisan915@gmail.com'
    req=urllib.request.Request(url,headers={'User-Agent':'protac-research/1.0 (mailto:kimjisan915@gmail.com)'})
    with urllib.request.urlopen(req,timeout=20) as r:
        m=json.load(r)['message']
        for k in ('issued','published-print','published-online','published','created'):
            try:
                y=m[k]['date-parts'][0][0]
                if y: return int(y)
            except Exception: pass
    return None

n_new=0
for i,d in enumerate(dois):
    if d in cache and pd.notna(cache[d]): continue
    try:
        y=cr_year(d); cache[d]=y; n_new+=1
    except Exception as e:
        cache[d]=cache.get(d);
    if i%50==0:
        pd.DataFrame({'doi':list(cache),'year':list(cache.values())}).to_csv(OUT,index=False)
        print(f"  {i}/{len(dois)} (신규 {n_new})", flush=True)
    time.sleep(0.15)
pd.DataFrame({'doi':list(cache),'year':list(cache.values())}).to_csv(OUT,index=False)
ok=sum(1 for v in cache.values() if pd.notna(v))
print(f"[완료] {OUT} | 연도 확보 {ok}/{len(dois)} ({ok/len(dois):.0%})")
