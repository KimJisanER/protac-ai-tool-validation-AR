"""
Boltz cofold 29개(CRBN 정렬, 동일 시점) 타일을 B→C→A 순서로 3×10 바둑판 합성.
검은 배경 + 각 타일 좌상단 흰색 ID + 색 범례(AR 하늘색/CRBN 주황/ligand 초록).
타일은 render_tiles.py(PyMOL)가 outputs/tiles/<ID>.png 로 생성.
출력: outputs/cofold_grid_BCA.png
"""
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, matplotlib.image as mpimg
TILE='/home/kimjisan95/ar_protac_project/outputs/tiles'
OUT='/home/kimjisan95/ar_protac_project/outputs/cofold_grid_BCA.png'
order=[f'B{i}' for i in range(1,7)]+[f'C{i}' for i in range(1,8)]+[f'A{i}' for i in range(1,17)]  # 29
ncol,nrow=10,3
fig,axes=plt.subplots(nrow,ncol,figsize=(ncol*1.7,nrow*1.7+0.4))
fig.patch.set_facecolor('black')
for k,ax in enumerate(axes.flat):
    ax.set_facecolor('black'); ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_visible(False)
    if k<len(order):
        cid=order[k]; ax.imshow(mpimg.imread(f'{TILE}/{cid}.png'))
        ax.text(0.04,0.96,cid,color='white',fontsize=12,fontweight='bold',va='top',ha='left',transform=ax.transAxes)
    else:
        ax.axis('off')
# 색 범례(흰 글씨)
fig.text(0.5,0.985,'AR-LBD = skyblue    CRBN = orange    PROTAC = green   (CRBN-aligned, common view)',
         color='white',ha='center',va='top',fontsize=11)
plt.subplots_adjust(left=0.002,right=0.998,top=0.965,bottom=0.005,wspace=0.02,hspace=0.02)
plt.savefig(OUT,dpi=150,facecolor='black'); print('[저장]',OUT)
