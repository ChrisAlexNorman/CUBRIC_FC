import numpy as np
from copy import deepcopy
from brainspace.gradient import GradientMaps
from corr_functions import corr_single, corr_surrogates, remove_unconnected_rois, sig_test_surrogates, get_mean_grads
from visualisation import plot_corr_coeff_variance

###### SETTINGS ######
n_surr = 1000
n_grads = 5
space = 'fsaverage5'
parcellation = 'glasser-360'
session = 'rest'
matrix = 'FC'
alignment = 'procrustes'
approaches = ['pca','dm','le']
kernels = ['cosine','normalized_angle','gaussian']
thresholds = [0,0.8,0.9]
random_seed = 0

PROJ_DIR = 'C:/Users/Chris/OneDrive/Documents/CUBRIC/ALSPAC/CardiffFC/'
sc_sample_names = ['alspac','camcan_gr_0','camcan_gr_1','camcan_gr_2']
fc_mean_fil = PROJ_DIR+'data/mica_processed/micapipe/mean_surfs/func/rest/fsaverage5_glasser-360_FC.txt'
grad_dir = PROJ_DIR+'data/gradients/'+space+'/'

###### AVERAGE FC METHOD ######
print('------Average FC Method------')
fc_mean_raw = np.loadtxt(open(fc_mean_fil,"r"),delimiter=',')
# Remove non-cortical regions
fc_mean_raw = fc_mean_raw[49:,49:]
# Zero diagonal
np.fill_diagonal(fc_mean_raw,0)
# Record number of ROIs
n_reg = np.shape(fc_mean_raw)[0]

# New fig for each threshold
for idx_threshold, threshold in enumerate(thresholds):

    # Initialise results objects
    r_all_methods = np.zeros((len(approaches),len(kernels),len(sc_sample_names),n_grads))
    r_surr_all_methods = np.zeros((len(approaches),len(kernels),len(sc_sample_names),n_surr,n_grads))
    p_all_methods = np.zeros((len(approaches),len(kernels),len(sc_sample_names),n_grads))

    # Apply whole-matrix thresholding of FC
    fc_mean = deepcopy(fc_mean_raw)
    thresh_val = np.sort(fc_mean.flatten())[int(threshold*n_reg**2)]
    fc_mean[fc_mean < thresh_val] = 0
    # Remove unocnnected rois (if any were removed by thresholding)
    fc_mean, ind_zero = remove_unconnected_rois(fc_mean)

    # Loop over gradient calculation methods
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
            
            # Loop over SC data samples
            for idx_sample, sample in enumerate(sc_sample_names):

                # Load SC ISV data
                sc_isv_fil = PROJ_DIR+'data/structural/'+sample+'_isv_sc_z_noarea.csv'
                sc_isv = np.loadtxt(open(sc_isv_fil,"r"),delimiter=' ')
                sc_isv_surr_fil = PROJ_DIR+'data/structural/'+sample+'_isv_sc_z_noarea_'+str(n_surr)+'_surrogates.csv'
                sc_isv_surr = np.loadtxt(open(sc_isv_surr_fil,"r"),delimiter=',')

                # Get raw correlations
                r, p = corr_single(sc_isv,grads,n_grads=n_grads)
                r = np.asarray(r)
                p = np.asarray(p)

                # Get correlations using surrogate maps
                r_surr, p_surr = corr_surrogates(sc_isv_surr,grads,n_grads=n_grads)

                # Get significance
                p_SA_removed = sig_test_surrogates(r,r_surr,test='normal')

                r_all_methods[idx_approach,idx_kernel,idx_sample] = r
                r_surr_all_methods[idx_approach,idx_kernel,idx_sample] = r_surr
                p_all_methods[idx_approach,idx_kernel,idx_sample] = p_SA_removed

    # Illustrate results
    fig = plot_corr_coeff_variance(r_all_methods,labels=sc_sample_names,title=str(int(threshold*100))+'% Thresholding')
    fig.savefig(PROJ_DIR+'gradients/figures/aggregate_corr_grid_samples_threshold_'+str(threshold)+'_mean_FC_'+str(n_surr)+'_surrogates.png',bbox_inches="tight")


###### AVERAGE GRAD METHOD ######
print('\n------Average Grad Method------')

# New fig for each threshold
for idx_threshold, threshold in enumerate(thresholds):

    # Initialise results objects
    r_all_methods = np.zeros((len(approaches),len(kernels),len(sc_sample_names),n_grads))
    r_surr_all_methods = np.zeros((len(approaches),len(kernels),len(sc_sample_names),n_surr,n_grads))
    p_all_methods = np.zeros((len(approaches),len(kernels),len(sc_sample_names),n_grads))

    # Loop over gradient calculation methods
    for idx_approach, approach in enumerate(approaches):
        for idx_kernel, kernel in enumerate(kernels):

            # Get gradients
            grad_fil = grad_dir + space+'_'+parcellation+'_'+session+'_'+matrix+'_'+approach+'_'+kernel+'_'+alignment+'_'+str(int(threshold*100))+'_aligned-grads.npy'
            grads = get_mean_grads(grad_fil,standardize=True)

            # Loop over SC data samples
            for idx_sample, sample in enumerate(sc_sample_names):

                # Load SC ISV data
                sc_isv_fil = PROJ_DIR+'data/structural/'+sample+'_isv_sc_z_noarea.csv'
                sc_isv = np.loadtxt(open(sc_isv_fil,"r"),delimiter=' ')
                sc_isv_surr_fil = PROJ_DIR+'data/structural/'+sample+'_isv_sc_z_noarea_'+str(n_surr)+'_surrogates.csv'
                sc_isv_surr = np.loadtxt(open(sc_isv_surr_fil,"r"),delimiter=',')

                # Get raw correlations
                r, p = corr_single(sc_isv,grads,n_grads=n_grads)
                r = np.asarray(r)
                p = np.asarray(p)

                # Get correlations using surrogate maps
                r_surr, p_surr = corr_surrogates(sc_isv_surr,grads,n_grads=n_grads)

                # Get significance
                p_SA_removed = sig_test_surrogates(r,r_surr,test='normal')

                r_all_methods[idx_approach,idx_kernel,idx_sample] = r
                r_surr_all_methods[idx_approach,idx_kernel,idx_sample] = r_surr
                p_all_methods[idx_approach,idx_kernel,idx_sample] = p_SA_removed

    # Illustrate results
    fig = plot_corr_coeff_variance(r_all_methods,labels=sc_sample_names,title=str(int(threshold*100))+'% Thresholding')
    fig.savefig(PROJ_DIR+'gradients/figures/aggregate_corr_grid_samples_threshold_'+str(threshold)+'_mean_grads_'+str(n_surr)+'_surrogates.png',bbox_inches="tight")

