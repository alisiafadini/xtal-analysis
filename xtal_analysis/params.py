import gemmi as gm

#### For entire dataset ####

h_res         = 1.63
l_res         = 20
alpha         = 0.05
sampling      = 0.41  ##specifies map grid
map_res       = 4 ## specifies map grid for GEMMI map_from_Fs function
name          = '100ps-400nm' #dataset name

####  For map masks #####

chrom_center  = [-22.46, 2.33, -3.41]
chrom_radius  = 5.4


####  For solvent masks #####

masker         = gm.SolventMasker(gm.AtomicRadiiSet.Constant, 1.5)
masker.rprobe  = 0.9
masker.rshrink = 1.1


#### Writing out maps ###

cell         = [51.99, 62.91, 72.03, 90, 90, 90]
space_group  = 'P212121'
#grid_size    = [150, 180, 216]
grid_size    = [128, 160, 180]

#### Paths for correlation pipeline ####

loop_path    = '/Users/alisia/Desktop/SACLA_JUN2021/final_scripts/background_subtracted_maps'
pdb_path     = '/Users/alisia/Desktop/SACLA_JUN2021/final_scripts/background_subtracted_maps/dark.pdb'



