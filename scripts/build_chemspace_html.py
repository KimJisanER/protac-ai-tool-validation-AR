"""
인터랙티브 chemical space HTML 2종 생성 (Plotly scattergl + RDKit.js 호버 구조).
- View1: 전체 PROTAC-DB(회색) + AR 타겟(P10275, 초록)
- View2: AI 검증 hold-out(prospective 1134, 주황) vs 나머지(회색)   ← View1과 색 구분
각 그룹 체크박스로 on/off(시각화 제외) 가능. 좌표는 캐시(재계산 생략).
출력: chemspace_view1_AR.html, chemspace_view2_split.html (프로젝트 루트)
"""
import os, json, html, numpy as np, pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors as rd
RDLogger.DisableLog('rdApp.*')

ROOT='/home/kimjisan95/ar_protac_project'
CACHE=f'{ROOT}/data/processed/chemspace_coords.csv'

def ik(s):
    m=Chem.MolFromSmiles(str(s)); return Chem.MolToInchiKey(m)[:14] if m else None

if os.path.exists(CACHE):
    print('[1-4] 좌표 캐시 로드:', CACHE)
    df=pd.read_csv(CACHE)
    df['targets']=df['targets'].fillna('')
else:
    print('[1] protac.xlsx 로드 + ik14 dedup')
    x=pd.read_excel(f'{ROOT}/protac.xlsx')
    x['ik14']=x['Smiles'].map(ik); x=x.dropna(subset=['ik14'])
    ar_ik=set(x.loc[x['Uniprot']=='P10275','ik14'])
    g=x.groupby('ik14'); uni=g.first().reset_index()
    uni['targets']=uni['ik14'].map(g['Target'].apply(lambda s:'/'.join(sorted(set(s.dropna().astype(str)))[:3])))
    print(f'    unique 화합물 {len(uni)} (AR {len(ar_ik)})')
    print('[2] held-out(prospective) 골격 로드')
    held=set(pd.read_csv(f'{ROOT}/data/processed/timesplit_compounds.csv')['ik14'])
    seen=set()
    try:
        fine=pd.read_csv('/home/kimjisan95/PROTAC-STAN/data/PROTAC-fine/protac-fine.csv')
        seen=set(fine['Smiles'].map(ik).dropna())
    except Exception: pass
    print('[3] Morgan FP(r2,1024) + 물성')
    fps=[]; keep=[]
    for _,row in uni.iterrows():
        m=Chem.MolFromSmiles(str(row['Smiles']))
        if m is None: continue
        fp=AllChem.GetMorganFingerprintAsBitVect(m,2,1024)
        arr=np.zeros((1024,),dtype=np.int8); AllChem.DataStructs.ConvertToNumpyArray(fp,arr)
        fps.append(arr)
        keep.append({'ik14':row['ik14'],'smiles':str(row['Smiles']),
                     'targets':row['targets'] if isinstance(row['targets'],str) else '',
                     'mw':round(Descriptors.MolWt(m),1),'tpsa':round(rd.CalcTPSA(m),1),
                     'is_ar':int(row['ik14'] in ar_ik),'is_held':int(row['ik14'] in held),
                     'is_seen':int(row['ik14'] in seen)})
    X=np.vstack(fps); df=pd.DataFrame(keep)
    print(f'    FP {X.shape}')
    print('[4] UMAP(Jaccard) 2D 투영')
    import umap
    emb=umap.UMAP(n_neighbors=15,min_dist=0.1,metric='jaccard',random_state=42).fit_transform(X)
    df['x']=np.round(emb[:,0],3); df['y']=np.round(emb[:,1],3)
    df.to_csv(CACHE,index=False); print('    [캐시 저장]', CACHE)

# ===== HTML 빌더 =====
def build_html(title, subtitle, group_col, hi_label, gray_label, out, hi_rgba, hi_line):
    hi=df[df[group_col]==1]; gr=df[df[group_col]==0]
    def pack(d):
        return {'x':d['x'].tolist(),'y':d['y'].tolist(),
                'cd':[[r.smiles,r.targets,r.mw,r.tpsa,r.ik14] for r in d.itertuples()]}
    payload={'hi':pack(hi),'gray':pack(gr),
             'hi_label':f'{hi_label} (n={len(hi)})','gray_label':f'{gray_label} (n={len(gr)})'}
    js=json.dumps(payload)
    page=f"""<!DOCTYPE html><html lang="ko"><head><meta charset="utf-8">
<title>{html.escape(title)}</title>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<script src="https://unpkg.com/@rdkit/rdkit/dist/RDKit_minimal.js"></script>
<style>
 body{{margin:0;font-family:-apple-system,Segoe UI,Roboto,'Noto Sans KR',sans-serif;background:#f8fafc;color:#0f172a}}
 header{{padding:14px 20px;background:#fff;border-bottom:1px solid #e2e8f0}}
 h1{{margin:0;font-size:18px}} .sub{{color:#64748b;font-size:13px;margin-top:4px}}
 #wrap{{display:flex;height:calc(100vh - 64px)}} #plot{{flex:1;min-width:0}}
 #side{{width:380px;border-left:1px solid #e2e8f0;background:#fff;padding:16px;overflow:auto}}
 .panel{{border:1px solid #e2e8f0;border-radius:10px;padding:12px;margin-bottom:14px;background:#fff}}
 .ptitle{{font-weight:700;font-size:13px;margin-bottom:8px}}
 .tog{{display:flex;align-items:center;gap:8px;font-size:13px;margin:6px 0;cursor:pointer}}
 .tog input{{width:16px;height:16px;cursor:pointer}} .dot{{font-size:15px;line-height:1}}
 #mol{{min-height:270px;display:flex;align-items:center;justify-content:center}} #mol svg{{max-width:100%;height:auto}}
 .info{{margin-top:10px;font-size:13px;line-height:1.6}} .info b{{color:#334155}}
 .smi{{word-break:break-all;color:#475569;font-size:11px;background:#f1f5f9;padding:8px;border-radius:8px;margin-top:8px}}
 .hint{{color:#94a3b8;font-size:12px;margin-top:10px}}
</style></head><body>
<header><h1>{html.escape(title)}</h1><div class="sub">{html.escape(subtitle)}</div></header>
<div id="wrap"><div id="plot"></div>
 <div id="side">
   <div class="panel"><div class="ptitle">그룹 표시 on/off</div>
     <label class="tog"><input type="checkbox" checked onchange="toggle(1,this.checked)"><span class="dot" style="color:{hi_line}">●</span><span>{html.escape(payload['hi_label'])}</span></label>
     <label class="tog"><input type="checkbox" checked onchange="toggle(0,this.checked)"><span class="dot" style="color:#94a3b8">●</span><span>{html.escape(payload['gray_label'])}</span></label>
   </div>
   <div class="panel"><div class="ptitle">마우스 호버 화합물</div>
     <div id="mol"><span style="color:#94a3b8;font-size:13px">점 위에 마우스를 올리세요</span></div>
     <div class="info" id="info"></div>
   </div>
   <div class="hint">스크롤=확대 · 드래그=이동 · 더블클릭=리셋 · 체크박스/범례=그룹 표시 토글</div>
 </div></div>
<script>
const DATA={js};
let RD=null;
window.initRDKitModule({{locateFile:p=>p.endsWith('.wasm')?'https://unpkg.com/@rdkit/rdkit/dist/RDKit_minimal.wasm':p}}).then(m=>{{RD=m;}});
const traces=[
 {{x:DATA.gray.x,y:DATA.gray.y,customdata:DATA.gray.cd,mode:'markers',type:'scattergl',
   name:DATA.gray_label,marker:{{size:3,color:'rgba(148,163,184,0.45)'}},hoverinfo:'none'}},
 {{x:DATA.hi.x,y:DATA.hi.y,customdata:DATA.hi.cd,mode:'markers',type:'scattergl',
   name:DATA.hi_label,marker:{{size:6,color:'{hi_rgba}',line:{{width:0.5,color:'{hi_line}'}}}},hoverinfo:'none'}}
];
const layout={{margin:{{l:30,r:10,t:10,b:30}},xaxis:{{title:'UMAP-1',zeroline:false}},
 yaxis:{{title:'UMAP-2',zeroline:false}},legend:{{x:0.01,y:0.99,bgcolor:'rgba(255,255,255,0.7)'}},
 hovermode:'closest',dragmode:'pan'}};
Plotly.newPlot('plot',traces,layout,{{responsive:true,scrollZoom:true,displaylogo:false}});
const plot=document.getElementById('plot');
function toggle(i,on){{Plotly.restyle('plot',{{visible:on?true:false}},[i]);}}
function show(cd){{
 const [smi,tg,mw,tpsa,ikk]=cd;
 if(RD){{try{{const mol=RD.get_mol(smi); if(mol){{document.getElementById('mol').innerHTML=mol.get_svg(360,260); mol.delete();}}}}catch(e){{}}}}
 document.getElementById('info').innerHTML='<b>Target(s):</b> '+(tg||'-')+'<br><b>MW:</b> '+mw+' &nbsp; <b>TPSA:</b> '+tpsa+'<br><b>InChIKey14:</b> '+ikk+'<div class="smi">'+smi+'</div>';
}}
plot.on('plotly_hover',e=>{{if(e.points&&e.points.length)show(e.points[0].customdata);}});
</script></body></html>"""
    with open(f'{ROOT}/{out}','w') as f: f.write(page)
    print(f'    [저장] {out}  (hi {len(hi)} / gray {len(gr)})')

print('[5] HTML 생성')
build_html('PROTAC-DB 화학공간 — AR 타겟 강조',
  '회색=전체 PROTAC-DB 3.0 · 초록=AR(UniProt P10275) 타겟 PROTAC · Morgan FP(r2,1024)→UMAP(Jaccard)',
  'is_ar','AR 타겟 PROTAC (P10275)','전체 PROTAC-DB','chemspace_view1_AR.html',
  'rgba(34,197,94,0.85)','#15803d')                      # 초록
build_html('PROTAC-DB 화학공간 — AI 검증 hold-out 분포',
  '회색=나머지(학습/기타) · 주황=AI 모델 검증용 hold-out(prospective; STAN·Ribes 미관측) · View1과 동일 UMAP 좌표',
  'is_held','검증 hold-out (prospective)','나머지 PROTAC-DB','chemspace_view2_split.html',
  'rgba(249,115,22,0.85)','#c2410c')                     # 주황 (View1 초록과 구분)
print('[완료]')
