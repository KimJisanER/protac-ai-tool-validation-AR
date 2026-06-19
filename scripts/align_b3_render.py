from pymol import cmd
cmd.load("/home/kimjisan95/PROTAC_MTL_v5/boltz_all29/out/boltz_results_inputs/predictions/B3/B3_model_0.pdb","B3")
cmd.load("2PNU.pdb","AR_exp")     # best AR-LBD (RMSD 0.38)
cmd.load("4CI3.pdb","CRBN_exp")   # best CRBN (RMSD 0.53, pomalidomide)
cmd.super("AR_exp","B3 and chain A and polymer")
cmd.super("CRBN_exp","B3 and chain B and polymer")
# 4CI3는 DDB1+CRBN → B3 CRBN 근처 폴리머만 = 실제 CRBN
cmd.select("crbn_only","CRBN_exp and polymer within 12 of (B3 and chain B)")
cmd.select("crbn_chain","byres (CRBN_exp and polymer within 6 of crbn_only)")
cmd.bg_color("white")
cmd.hide("everything")
cmd.set("ray_opaque_background",0); cmd.set("cartoon_transparency",0.0); cmd.set("ray_shadows",0); cmd.set("antialias",2)
def ligsel(obj): return f"{obj} and not polymer and not solvent and not inorganic"
# ---- overview ----
cmd.show("cartoon","B3 and chain A"); cmd.color("marine","B3 and chain A")
cmd.show("cartoon","B3 and chain B"); cmd.color("orange","B3 and chain B")
cmd.show("sticks","B3 and chain C"); cmd.color("yellow","B3 and chain C")
cmd.show("cartoon","AR_exp and polymer"); cmd.color("palegreen","AR_exp and polymer"); cmd.set("cartoon_transparency",0.55,"AR_exp")
cmd.show("cartoon","crbn_chain"); cmd.color("lightpink","crbn_chain"); cmd.set("cartoon_transparency",0.55,"crbn_chain")
cmd.orient("B3 and (chain A or chain B)")
cmd.ray(1500,1100); cmd.png(f"/home/kimjisan95/ar_protac_project/outputs/align_B3_overview.png", dpi=150)
# ---- AR pocket zoom: B3 워헤드 vs 2PNU 리간드 ----
cmd.hide("everything")
cmd.show("cartoon","B3 and chain A"); cmd.color("marine","B3 and chain A"); cmd.set("cartoon_transparency",0.6,"B3")
cmd.show("cartoon","AR_exp and polymer"); cmd.color("palegreen","AR_exp"); cmd.set("cartoon_transparency",0.6,"AR_exp")
cmd.show("sticks", ligsel("AR_exp")); cmd.color("cyan", ligsel("AR_exp"))     # 실험 AR 리간드
cmd.show("sticks","B3 and chain C"); cmd.color("yellow","B3 and chain C")     # B3 PROTAC
cmd.orient(ligsel("AR_exp"))
cmd.zoom(ligsel("AR_exp"), 8)
cmd.ray(1300,1000); cmd.png(f"/home/kimjisan95/ar_protac_project/outputs/align_AR_pocket.png", dpi=150)
# ---- CRBN pocket zoom: B3 글루타리미드 vs 4CI3 pomalidomide ----
cmd.hide("everything")
cmd.show("cartoon","crbn_chain"); cmd.color("lightpink","crbn_chain"); cmd.set("cartoon_transparency",0.6,"CRBN_exp")
cmd.show("cartoon","B3 and chain B"); cmd.color("orange","B3 and chain B"); cmd.set("cartoon_transparency",0.6,"B3")
cmd.show("sticks", ligsel("CRBN_exp")); cmd.color("magenta", ligsel("CRBN_exp"))
cmd.show("sticks","B3 and chain C"); cmd.color("yellow","B3 and chain C")
cmd.orient(ligsel("CRBN_exp"))
cmd.zoom(ligsel("CRBN_exp"), 8)
cmd.ray(1300,1000); cmd.png(f"/home/kimjisan95/ar_protac_project/outputs/align_CRBN_pocket.png", dpi=150)
print("RENDER_DONE")
