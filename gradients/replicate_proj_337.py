import numpy as np
from brainspace.gradient import GradientMaps
from corr_functions import corr_single, corr_surrogates, remove_unconnected_rois
from visualisation import basic_one_plot

# Settings
n_surr = int(10**4)
ci = 95
sig_level = 0.01
PROJ_DIR = '/home/cnorman/Documents/CUBRIC/ALSPAC/CardiffFC/'
sc_isv_fil = PROJ_DIR + 'data/structural/isv_sc_z_noarea.csv'
sc_isv_surr_fil = PROJ_DIR + 'data/structural/isv_sc_z_noarea_'+str(n_surr)+'_surrogates.csv'
fc_337_fil = PROJ_DIR + 'data/proj_337/FC_fmri_29subj_thr.csv'
fc_alspac_fil = PROJ_DIR + 'data/mica_processed/micapipe/mean_surfs/func/rest/fsaverage5_glasser-360_FC.txt'

# Proj 337 settings
fc_threshold = 0.0 # Trying to match proj 337
n_components = 10
approach = 'dm'
kernel = 'normalized_angle'
random_state = 0

# Load SC ISV data
sc_isv = np.loadtxt(open(sc_isv_fil,"r"),delimiter=' ')
n_reg = np.shape(sc_isv)[0]
sc_isv_surr = np.loadtxt(open(sc_isv_surr_fil,"r"),delimiter=',')


############ Proj 337 FC ############
print('------Proj. 337 Correlations------')
fc_mean_337_raw = np.loadtxt(open(fc_337_fil,"r"),delimiter=',')
# Remove unconnected regions in proj 337 matrix to protect gradient estimation
fc_mean_337, ind_zero_337 = remove_unconnected_rois(fc_mean_337_raw)

# Get gradients
gm_337 = GradientMaps(n_components = n_components,
        approach = approach,
        kernel = kernel,
        random_state = random_state)
gm_337.fit(fc_mean_337)
# replace unconnected indices
grads_337 = np.zeros((n_reg,n_components)) * np.nan
ind_valid_alspac = np.setdiff1d(np.arange(0, n_reg), ind_zero_337)   
grads_337[np.ix_(ind_valid_alspac, np.arange(0, n_components))] = gm_337.gradients_  

# Get raw correlations
corrs_r_337, corrs_p_337 = corr_single(sc_isv,grads_337,n_grads=n_components)
corrs_r_337 = np.asarray(corrs_r_337)
corrs_p_337 = np.asarray(corrs_p_337)
sig_grads_337 = [ind+1 for ind in range(len(corrs_p_337)) if corrs_p_337[ind] < sig_level]
print('rho:',corrs_r_337)
print('p:',corrs_p_337)

# Check correlations using surrogate maps
corrs_r_337_surr, corrs_p_337_surr = corr_surrogates(sc_isv_surr,grads_337,n_grads=n_components)

# Illustrate
fig_337, ax_337 = basic_one_plot(corr=corrs_r_337,p=corrs_p_337,corr_surr=corrs_r_337_surr,ci=ci,sig=sig_level)
ax_337.set_title('Proj. 337')
fig_337.savefig(PROJ_DIR+'gradients/figures/Corr_from_mean_FCs_proj_337_'+str(n_surr)+'_surrogates.png',bbox_inches="tight")


############ Alspac FC ############
print('\n------Alspac Correlations------')
fc_mean_alspac_raw = np.loadtxt(open(fc_alspac_fil,"r"),delimiter=',')
# Remove non-cortical regions
fc_mean_alspac = fc_mean_alspac_raw[49:,49:]
# Zero diagonal
np.fill_diagonal(fc_mean_alspac,0)
# Apply whole-matrix thresholding
threshold = np.sort(fc_mean_alspac.flatten())[int(fc_threshold*n_reg**2)]
fc_mean_alspac[fc_mean_alspac < threshold] = 0
# Remove unocnnected rois (if any were removed by thresholding)
fc_mean_alspac, ind_zero_alspac = remove_unconnected_rois(fc_mean_alspac)

# Get gradients
gm_alspac = GradientMaps(n_components = n_components,
        approach = approach,
        kernel = kernel,
        random_state = random_state)
gm_alspac .fit(fc_mean_alspac)
# replace unconnected indices
grads_alspac = np.zeros((n_reg,n_components)) * np.nan
ind_valid_alspac = np.setdiff1d(np.arange(0, n_reg), ind_zero_alspac)   
grads_alspac[np.ix_(ind_valid_alspac, np.arange(0, n_components))] = gm_alspac.gradients_

# Get raw correlations
corrs_r_alspac, corrs_p_alspac = corr_single(sc_isv,grads_alspac,n_grads=n_components)
corrs_r_alspac = np.asarray(corrs_r_alspac)
corrs_p_alspac = np.asarray(corrs_p_alspac)
sig_grads_alspac = [ind+1 for ind in range(len(corrs_p_alspac)) if corrs_p_alspac[ind] < sig_level]
print('rho:',corrs_r_alspac)
print('p:',corrs_p_alspac)

# Check correlations using surrogate maps
corrs_r_alspac_surr, corrs_p_alspac_surr = corr_surrogates(sc_isv_surr,grads_alspac,n_grads=n_components)

# Illustrate
fig_alspac, ax_alspac = basic_one_plot(corr=corrs_r_alspac,p=corrs_p_alspac,corr_surr=corrs_r_alspac_surr,ci=ci,sig=sig_level)
ax_alspac.set_title('ALSPAC')
fig_alspac.savefig(PROJ_DIR+'gradients/figures/Corr_from_mean_FCs_ALSPAC_thresh_'+str(fc_threshold)+'_'+str(n_surr)+'_surrogates.png',bbox_inches="tight")

