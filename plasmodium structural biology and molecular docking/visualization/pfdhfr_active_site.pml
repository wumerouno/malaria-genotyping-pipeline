reinitialize
set retain_order, 1
set antialias, 2
set ray_trace_mode, 1
set ray_opaque_background, off
set cartoon_fancy_helices, 1
set stick_radius, 0.18
set ambient, 0.45
set spec_reflect, 0.2
set label_size, 12
set label_position, (1.4, 1.2, 1.0)
set label_color, black
bg_color white
load data/prepared/receptor/1J3I_chainA_protein.pdb, receptor
hide everything, all
show cartoon, receptor
color gray80, receptor
select active_site, receptor and chain A and resi 14+15+16+46+48+49+54+55+57+58+108+111+112+113+116+119+164+165+170+185
show sticks, active_site
color teal, active_site
label active_site and name CA and resi 14+15+54+58+108+111+164+170, resn + ' ' + resi
load data/prepared/receptor/1J3I_native_wra_chainA.pdb, native_wr99210
show sticks, native_wr99210
color yelloworange, native_wr99210
set stick_radius, 0.28, native_wr99210
load results/poses/wr99210_docked.pdb, docked_wr99210
show sticks, docked_wr99210
color tv_orange, docked_wr99210
set stick_radius, 0.24, docked_wr99210
load results/poses/pyrimethamine_docked.pdb, docked_pyrimethamine
show sticks, docked_pyrimethamine
color magenta, docked_pyrimethamine
set stick_radius, 0.24, docked_pyrimethamine
load results/poses/cycloguanil_docked.pdb, docked_cycloguanil
show sticks, docked_cycloguanil
color marine, docked_cycloguanil
set stick_radius, 0.24, docked_cycloguanil
load results/poses/methotrexate_docked.pdb, docked_methotrexate
show sticks, docked_methotrexate
color forest, docked_methotrexate
set stick_radius, 0.24, docked_methotrexate
remove hydro
select ligand_cluster, native_wr99210 or docked_wr99210 or docked_pyrimethamine or docked_cycloguanil or docked_methotrexate
set transparency, 0.86
show surface, active_site
color palecyan, active_site
show sticks, active_site
color teal, active_site
orient ligand_cluster or active_site
center ligand_cluster
zoom (ligand_cluster or active_site), 4
clip slab, 28
turn x, 15
turn y, -18
ray 1400, 1000
png results/figures/pfdhfr_active_site_docking.png, dpi=300
save results/figures/pfdhfr_active_site_docking.pse
quit
