from pymol import cmd
ALN='/home/kimjisan95/ar_protac_project/outputs/cofold_aligned'
OUT='/home/kimjisan95/ar_protac_project/outputs/tiles'
order=[f'B{i}' for i in range(1,7)]+[f'C{i}' for i in range(1,8)]+[f'A{i}' for i in range(1,17)]
for cid in order: cmd.load(f'{ALN}/{cid}.pdb', cid)
cmd.hide('everything'); cmd.show('cartoon','chain A or chain B'); cmd.show('sticks','chain C')
cmd.color('skyblue','chain A'); cmd.color('orange','chain B'); cmd.color('green','chain C')
cmd.set('stick_radius',0.25)
cmd.bg_color('black'); cmd.set('ray_opaque_background',1); cmd.set('ray_shadows',0); cmd.set('antialias',2)
cmd.enable('all'); cmd.orient('all'); cmd.zoom('all', 3)
view=cmd.get_view()
for cid in order:
    cmd.disable('all'); cmd.enable(cid); cmd.set_view(view)
    cmd.ray(560,560); cmd.png(f'{OUT}/{cid}.png', dpi=100)
print('TILES_DONE', len(order))
