import numpy as np
import matplotlib.pyplot as plt
from brainspace.gradient import GradientMaps
from scipy.stats import spearmanr

# Settings
n_surr = int(10**4)
ci = 95
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

def remove_unconnected_rois(FC):
    """Remove, from both hemispheres, regions that are not connected to any other region, return reduced connectivity matrix and array of unconnected roi indices"""

    nreg = np.shape(FC)[0]
    ind_zero = np.nonzero(FC.sum(axis=0) == 0)[0]

    # Remove corresponding region from other hemisphere
    ind_zero_pair = []
    for a in ind_zero:
        ind_zero_pair.append(a)
        if a < nreg / 2:
            ind_zero_pair.append(a + int(nreg / 2))
        else:
            ind_zero_pair.append(a - int(nreg / 2))

    ind_zero_sym = np.unique(np.array(ind_zero_pair))
    print('Discarding indices: ', ind_zero_sym)
    valid_ind = np.setdiff1d(np.arange(0, nreg), ind_zero_sym)
    FC_reduced = FC[np.ix_(valid_ind, valid_ind)]

    return FC_reduced, ind_zero_sym

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
corrs_r_337 = np.zeros(n_components)
corrs_p_337 = np.zeros(n_components)
for grad in range(0,n_components):
    corrs_r_337[grad], corrs_p_337[grad] = spearmanr(sc_isv,grads_337[:,grad],nan_policy='omit')
print('rho:',corrs_r_337)
print('p:',corrs_p_337)
# Illustrate
fig_337, ax_337 = plt.subplots(3,1,figsize=(8,6))
ax = ax_337[0]
ax.plot(range(1,n_components+1),corrs_r_337,'k.-',ms=10)
ax.plot([1,n_components],[0,0],'k.--')
ax.set_xticks(range(1,n_components+1))
yl1 = ax.get_ylim()
ax.set_ylabel('Corr. Coeff.')
ax = ax_337[1]
ax.plot(range(1,n_components+1),corrs_p_337,'k.',ls='-',ms=10)
ax.set_xticks(range(1,n_components+1))
ax.set_ylim([0,1])
ax.set_ylabel('p value')
# ax.set_title('Proj. 337 - Raw SC-ISV')

# Check correlations using surrogate maps
corrs_r_337_surr = np.zeros((n_surr,n_components))
corrs_p_337_surr = np.zeros((n_surr,n_components))
for surr in range(0,n_surr):
    for grad in range(0,n_components):
        corrs_r_337_surr[surr,grad], corrs_p_337_surr[surr,grad] = spearmanr(sc_isv_surr[:,surr],grads_337[:,grad],nan_policy='omit')
corrs_r_337_mean = np.mean(corrs_r_337_surr,axis=0)
corrs_r_337_lower = np.percentile(corrs_r_337_surr,(100-ci)/2,axis=0)
corrs_r_337_upper = np.percentile(corrs_r_337_surr,(100+ci)/2,axis=0)
corrs_r_337_bounds = np.row_stack((corrs_r_337_mean-corrs_r_337_lower,corrs_r_337_upper-corrs_r_337_mean))
# Illustrate
ax = ax_337[2]
ax.errorbar(np.asarray(range(1,n_components+1)),corrs_r_337_mean,corrs_r_337_bounds,ls='none',color='k')
ax.plot(np.asarray(range(1,n_components+1)),corrs_r_337_mean,marker='.',markersize=10,ls='none',color='k')
ax.plot([1,n_components],[0,0],'k--')
ax.set_xticks(range(1,n_components+1))
yl2 = ax.get_ylim()
ax.set_xlabel('Component')
ax.set_ylabel('Corr. Coeff.')
# ax.set_title('Proj. 337 - Surrogate SC-ISV')
yl = [min([yl1[0],yl2[0]]),max([yl1[1],yl2[1]])]
for ax in [ax_337[0],ax_337[2]]:
    ax.set_ylim(yl)
plt.savefig(PROJ_DIR+'/gradients/Figures/Corr_from_mean_FCs_proj_337_'+str(n_surr)+'_surrogates.png')

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

# Get correlations
corrs_r_alspac  = np.zeros(n_components)
corrs_p_alspac = np.zeros(n_components)
for grad in range(0,n_components):
    corrs_r_alspac[grad], corrs_p_alspac[grad] = spearmanr(sc_isv,grads_alspac[:,grad],nan_policy='omit')
print('rho:',corrs_r_alspac)
print('p:',corrs_p_alspac)
# Illustrate
fig_alspac, ax_alspac = plt.subplots(3,1,figsize=(8,6))
ax = ax_alspac[0]
ax.plot(range(1,n_components+1),corrs_r_alspac,'k.-',ms=10)
ax.plot([1,n_components],[0,0],'k.--')
ax.set_xticks(range(1,n_components+1))
yl1 = ax.get_ylim()
ax.set_ylabel('Corr. Coeff.')
ax = ax_alspac[1]
ax.plot(range(1,n_components+1),corrs_p_alspac,'k.',ls='-',ms=10)
ax.set_xticks(range(1,n_components+1))
ax.set_ylim([0,1])
ax.set_ylabel('p value')
# ax.set_title('Proj. 337 - Raw SC-ISV')

# Check correlations using surrogate maps
corrs_r_alspac_surr = np.zeros((n_surr,n_components))
corrs_p_alspac_surr = np.zeros((n_surr,n_components))
for surr in range(0,n_surr):
    for grad in range(0,n_components):
        corrs_r_alspac_surr[surr,grad], corrs_p_alspac_surr[surr,grad] = spearmanr(sc_isv_surr[:,surr],grads_alspac[:,grad],nan_policy='omit')
corrs_r_alspac_mean = np.mean(corrs_r_alspac_surr,axis=0)
corrs_r_alspac_lower = np.percentile(corrs_r_alspac_surr,(100-ci)/2,axis=0)
corrs_r_alspac_upper = np.percentile(corrs_r_alspac_surr,(100+ci)/2,axis=0)
corrs_r_alspac_bounds = np.row_stack((corrs_r_alspac_mean-corrs_r_alspac_lower,corrs_r_alspac_upper-corrs_r_alspac_mean))
# Illustrate
ax = ax_alspac[2]
ax.errorbar(np.asarray(range(1,n_components+1)),corrs_r_alspac_mean,corrs_r_alspac_bounds,ls='none',color='k')
ax.plot(np.asarray(range(1,n_components+1)),corrs_r_alspac_mean,marker='.',markersize=10,ls='none',color='k')
ax.plot([1,n_components],[0,0],'k--')
ax.set_xticks(range(1,n_components+1))
yl2 = ax.get_ylim()
ax.set_xlabel('Component')
ax.set_ylabel('Corr. Coeff.')
# ax.set_title('Proj. 337 - Surrogate SC-ISV')
yl = [min([yl1[0],yl2[0]]),max([yl1[1],yl2[1]])]
for ax in [ax_alspac[0],ax_alspac[2]]:
    ax.set_ylim(yl)
plt.savefig(PROJ_DIR+'/gradients/Figures/Corr_from_mean_FCs_Alspac_'+str(n_surr)+'_surrogates.png')
