# 실험 PDB는 RCSB에서 다운로드: AR(2AM9,1E3G,1T7R,2PNU,5T8E) CRBN(4CI1,4CI3,4TZ4,5FQD,8D7U)
#   curl -s -o <ID>.pdb https://files.rcsb.org/download/<ID>.pdb  → outputs/ref_pdb/
# 실행: pymol -cq scripts/align_b3_rank.py (RMSD 순위) / align_b3_render.py (스크린샷)
from pymol import cmd
cmd.load("/home/kimjisan95/PROTAC_MTL_v5/boltz_all29/out/boltz_results_inputs/predictions/B3/B3_model_0.pdb","B3")
ar=['2AM9','1E3G','1T7R','2PNU','5T8E']; crbn=['4CI1','4CI3','4TZ4','5FQD','8D7U']
for x in ar+crbn: cmd.load(x+'.pdb', x)
print("RANKAR")
for x in ar:
    r=cmd.super(x, "B3 and chain A and polymer")
    print(f"AR {x} RMSD {r[0]:.2f} atoms {r[1]}")
print("RANKCRBN")
for x in crbn:
    r=cmd.super(x, "B3 and chain B and polymer")
    print(f"CRBN {x} RMSD {r[0]:.2f} atoms {r[1]}")