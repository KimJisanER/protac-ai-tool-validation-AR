bg_color white
load /home/kimjisan95/ar_protac_project/outputs/cofold_aligned/A15.pdb, A15
load /home/kimjisan95/ar_protac_project/outputs/cofold_aligned/A8.pdb, A8
load /home/kimjisan95/ar_protac_project/outputs/cofold_aligned/B1.pdb, B1
load /home/kimjisan95/ar_protac_project/outputs/cofold_aligned/B2.pdb, B2
load /home/kimjisan95/ar_protac_project/outputs/cofold_aligned/B3.pdb, B3
load /home/kimjisan95/ar_protac_project/outputs/cofold_aligned/B4.pdb, B4
load /home/kimjisan95/ar_protac_project/outputs/cofold_aligned/B5.pdb, B5
load /home/kimjisan95/ar_protac_project/outputs/cofold_aligned/B6.pdb, B6
load /home/kimjisan95/ar_protac_project/outputs/cofold_aligned/C1.pdb, C1
load /home/kimjisan95/ar_protac_project/outputs/cofold_aligned/C2.pdb, C2
load /home/kimjisan95/ar_protac_project/outputs/cofold_aligned/C3.pdb, C3
load /home/kimjisan95/ar_protac_project/outputs/cofold_aligned/C5.pdb, C5
load /home/kimjisan95/ar_protac_project/outputs/cofold_aligned/C6.pdb, C6
load /home/kimjisan95/ar_protac_project/outputs/cofold_aligned/C7.pdb, C7
load /home/kimjisan95/ar_protac_project/outputs/cofold_aligned/A1.pdb, A1
load /home/kimjisan95/ar_protac_project/outputs/cofold_aligned/A10.pdb, A10
load /home/kimjisan95/ar_protac_project/outputs/cofold_aligned/A11.pdb, A11
load /home/kimjisan95/ar_protac_project/outputs/cofold_aligned/A12.pdb, A12
load /home/kimjisan95/ar_protac_project/outputs/cofold_aligned/A13.pdb, A13
load /home/kimjisan95/ar_protac_project/outputs/cofold_aligned/A14.pdb, A14
load /home/kimjisan95/ar_protac_project/outputs/cofold_aligned/A16.pdb, A16
load /home/kimjisan95/ar_protac_project/outputs/cofold_aligned/A2.pdb, A2
load /home/kimjisan95/ar_protac_project/outputs/cofold_aligned/A3.pdb, A3
load /home/kimjisan95/ar_protac_project/outputs/cofold_aligned/A4.pdb, A4
load /home/kimjisan95/ar_protac_project/outputs/cofold_aligned/A5.pdb, A5
load /home/kimjisan95/ar_protac_project/outputs/cofold_aligned/A6.pdb, A6
load /home/kimjisan95/ar_protac_project/outputs/cofold_aligned/A7.pdb, A7
load /home/kimjisan95/ar_protac_project/outputs/cofold_aligned/A9.pdb, A9
load /home/kimjisan95/ar_protac_project/outputs/cofold_aligned/C4.pdb, C4
hide everything
show cartoon
set cartoon_transparency, 0.2
color gray80, chain B
color marine, A15 and chain A
color marine, A8 and chain A
color marine, B1 and chain A
color marine, B2 and chain A
color marine, B3 and chain A
color marine, B4 and chain A
color marine, B5 and chain A
color marine, B6 and chain A
color marine, C1 and chain A
color marine, C2 and chain A
color marine, C3 and chain A
color marine, C5 and chain A
color marine, C6 and chain A
color marine, C7 and chain A
color orange, A1 and chain A
color orange, A10 and chain A
color orange, A11 and chain A
color orange, A12 and chain A
color orange, A13 and chain A
color orange, A14 and chain A
color orange, A16 and chain A
color orange, A2 and chain A
color orange, A3 and chain A
color orange, A4 and chain A
color orange, A5 and chain A
color orange, A6 and chain A
color orange, A7 and chain A
color orange, A9 and chain A
color orange, C4 and chain A
hide cartoon, chain C
show sticks, chain C
color yellow, chain C
set stick_radius, 0.12
set ray_opaque_background, 0
set antialias, 2
set ray_shadows, 0
orient
ray 1500, 1100
png /home/kimjisan95/ar_protac_project/outputs/cofold_clusters_overview.png, dpi=150
disable all
enable A15
enable A8
enable B1
enable B2
enable B3
enable B4
enable B5
enable B6
enable C1
enable C2
enable C3
enable C5
enable C6
enable C7
orient (enabled)
ray 1300,1000
png /home/kimjisan95/ar_protac_project/outputs/cofold_cluster1_BC.png, dpi=150
disable all
enable A1
enable A10
enable A11
enable A12
enable A13
enable A14
enable A16
enable A2
enable A3
enable A4
enable A5
enable A6
enable A7
enable A9
enable C4
orient (enabled)
ray 1300,1000
png /home/kimjisan95/ar_protac_project/outputs/cofold_cluster2_A.png, dpi=150
