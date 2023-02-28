import numpy as np
from brainspace.gradient import GradientMaps
from corr_functions import corr_single, corr_surrogates, remove_unconnected_rois, sig_test_surrogates, get_mean_grads
from visualisation import plot_methods_grid

###### SETTINGS ######
n_surr = 5000
n_grads = 6
space = 'fsaverage5'
parcellation = 'glasser-360'
session = 'rest'
matrix = 'FC'
alignment = 'procrustes'
approaches = ['pca','dm','le']
kernels = ['cosine','normalized_angle','gaussian']
thresholds = [0,0.8,0.9]
random_seed = 0

PROJ_DIR = '/home/cnorman/Documents/CUBRIC/ALSPAC/CardiffFC/'
sc_isv_fil = PROJ_DIR+'data/structural/alspac_isv_sc_z_noarea.csv'
sc_isv_surr_fil = PROJ_DIR+'data/structural/alspac_isv_sc_z_noarea_'+str(n_surr)+'_surrogates.csv'
fc_mean_fil = PROJ_DIR+'data/mica_processed/micapipe/mean_surfs/func/rest/fsaverage5_glasser-360_FC.txt'
grad_dir = PROJ_DIR+'data/gradients/'+space+'/'

# Load SC ISV data
sc_isv = np.loadtxt(open(sc_isv_fil,"r"),delimiter=' ')
n_reg = np.shape(sc_isv)[0]
sc_isv_surr = np.loadtxt(open(sc_isv_surr_fil,"r"),delimiter=',')

###### AVERAGE FC METHOD ######
print('------Average FC Method------')
fc_mean_raw = np.loadtxt(open(fc_mean_fil,"r"),delimiter=',')
# Remove non-cortical regions
fc_mean = fc_mean_raw[49:,49:]
# Zero diagonal
np.fill_diagonal(fc_mean,0)

# Loop over methods
r_all_methods = np.zeros((len(approaches),len(kernels),len(thresholds),n_grads))
r_surr_all_methods = np.zeros((len(approaches),len(kernels),len(thresholds),n_surr,n_grads))
p_all_methods = np.zeros((len(approaches),len(kernels),len(thresholds),n_grads))
for idx_threshold, threshold in enumerate(thresholds):

    # Apply whole-matrix thresholding
    thresh_val = np.sort(fc_mean.flatten())[int(threshold*n_reg**2)]
    fc_mean[fc_mean < thresh_val] = 0
    # Remove unocnnected rois (if any were removed by thresholding)
    fc_mean, ind_zero = remove_unconnected_rois(fc_mean)

    for idx_approach, approach in enumerate(approaches):
        for idx_kernel, kernel in enumerate(kernels):

            # Get gradients
            gm = GradientMaps(n_components = n_grads,
                    approach = approach,
                    kernel = kernel,
                    random_state = random_seed)
            gm .fit(fc_mean)
            # replace unconnected indices
            grads = np.zeros((n_reg,n_grads)) * np.nan
            ind_valid = np.setdiff1d(np.arange(0, n_reg), ind_zero)   
            grads[np.ix_(ind_valid, np.arange(0, n_grads))] = gm.gradients_
            # standardise
            for idx_grad in range(n_grads):
                col = grads[:,idx_grad]
                grads[:,idx_grad] = (col - np.mean(col)) / np.std(col)

            # Get raw correlations
            r, p = corr_single(sc_isv,grads,n_grads=n_grads)
            r = np.asarray(r)
            p = np.asarray(p)

            # Get correlations using surrogate maps
            r_surr, p_surr = corr_surrogates(sc_isv_surr,grads,n_grads=n_grads)

            # Get significance
            p_SA_removed = sig_test_surrogates(r,r_surr,test='normal')

            r_all_methods[idx_approach,idx_kernel,idx_threshold] = r
            r_surr_all_methods[idx_approach,idx_kernel,idx_threshold] = r_surr
            p_all_methods[idx_approach,idx_kernel,idx_threshold] = p_SA_removed

# Illustrate results
fig, ax = plot_methods_grid(r_all_methods,r_surr_all_methods,p_all_methods,approaches=approaches,kernels=kernels,thresholds=[str(thresh) for thresh in thresholds])

fig.savefig(PROJ_DIR+'gradients/figures/alspac_corr_grid_mean_FC_'+str(n_surr)+'_surrogates.png',bbox_inches="tight")


###### AVERAGE GRAD METHOD ######
print('\n------Average Grad Method------')

# Loop over methods
r_all_methods = np.zeros((len(approaches),len(kernels),len(thresholds),n_grads))
r_surr_all_methods = np.zeros((len(approaches),len(kernels),len(thresholds),n_surr,n_grads))
p_all_methods = np.zeros((len(approaches),len(kernels),len(thresholds),n_grads))
for idx_threshold in range(len(thresholds)):
    threshold = thresholds[idx_threshold]
    for idx_approach in range(len(approaches)):
        approach = approaches[idx_approach]
        for idx_kernel in range(len(kernels)):
            kernel = kernels[idx_kernel]

            # Get gradients
            grad_fil = grad_dir + space+'_'+parcellation+'_'+session+'_'+matrix+'_'+approach+'_'+kernel+'_'+alignment+'_'+str(int(threshold*100))+'_aligned-grads.npy'
            grads = get_mean_grads(grad_fil,standardize=True)

            # Get raw correlations
            r, p = corr_single(sc_isv,grads,n_grads=n_grads)
            r = np.asarray(r)
            p = np.asarray(p)

            # Get correlations using surrogate maps
            r_surr, p_surr = corr_surrogates(sc_isv_surr,grads,n_grads=n_grads)

            # Get significance
            p_SA_removed = sig_test_surrogates(r,r_surr,test='normal')

            r_all_methods[idx_approach,idx_kernel,idx_threshold] = r
            r_surr_all_methods[idx_approach,idx_kernel,idx_threshold] = r_surr
            p_all_methods[idx_approach,idx_kernel,idx_threshold] = p_SA_removed

# Illustrate results
fig, ax = plot_methods_grid(r_all_methods,r_surr_all_methods,p_all_methods,approaches=approaches,kernels=kernels,thresholds=[str(thresh) for thresh in thresholds])

fig.savefig(PROJ_DIR+'gradients/figures/alspac_corr_grid_mean_grads_'+str(n_surr)+'_surrogates.png',bbox_inches="tight")

