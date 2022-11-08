# Generate selected figures for SC ISV and FC gradients

# Settings
PROJ_DIR = '~/Documents/CUBRIC/ALSPAC/proj_cn/'
SC_ISV_FIL = PROJ_DIR + 'data/isv_sc_z_noarea.csv'

session = 'rest'
matrix = 'FC'
space = 'fsaverage5'
parcellation = 'glasser-360'
approach = 'dm'
kernel = 'normalized_angle'
alignment = 'procrustes'
threshold = 0.9


FC_GRAD_FIL = PROJ_DIR + 'data/gradients/' + space + '/' + space+'_'+parcellation+'_'+session+'_'+matrix+'_'+approach+'_'+kernel+'_'+alignment+'_'+str(int(threshold*100))+'_aligned-grads.npy'

