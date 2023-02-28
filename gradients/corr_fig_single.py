import numpy as np
from brainspace.gradient import GradientMaps
import matplotlib.pyplot as plt
from corr_functions import corr_single, corr_surrogates, remove_unconnected_rois, sig_test_surrogates, get_mean_grads
from visualisation import basic_one_plot

###### SETTINGS ######
n_surr = 1000
n_grads = 8
space = 'fsaverage5'
parcellation = 'glasser-360'
session = 'rest'
matrix = 'FC'
alignment = 'procrustes'
approach = 'dm'
kernel = 'normalized_angle'
threshold = 0
random_seed = 0

PROJ_DIR = 'C:/Users/Chris/OneDrive/Documents/CUBRIC/ALSPAC/CardiffFC/'
sc_sample_names = ['alspac','camcan_gr_0','camcan_gr_1','camcan_gr_2']
fc_mean_fil = PROJ_DIR+'data/mica_processed/micapipe/mean_surfs/func/rest/fsaverage5_glasser-360_FC.txt'
grad_dir = PROJ_DIR+'data/gradients/'+space+'/'

fc_mean = np.loadtxt(open(fc_mean_fil,"r"),delimiter=',')
# Remove non-cortical regions
fc_mean = fc_mean[49:,49:]
# Zero diagonal
np.fill_diagonal(fc_mean,0)
# Record number of ROIs
n_reg = np.shape(fc_mean)[0]
# Apply whole-matrix thresholding of FC
thresh_val = np.sort(fc_mean.flatten())[int(threshold*n_reg**2)]
fc_mean[fc_mean < thresh_val] = 0
# Remove unocnnected rois (if any were removed by thresholding)
fc_mean, ind_zero = remove_unconnected_rois(fc_mean)

# Get Gradients
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

# Plot Lambdas:
fig = plt.figure(figsize=(4,2))
ax = plt.axes()
ax.plot(range(1,n_grads+1),gm.lambdas_,'.',color='k',ms=10)
ylim = ax.get_ylim()
ax.set_ylim([0,ylim[1]])
ax.set_xticks(range(1,n_grads+1))
ax.set_xlabel('Component')
ax.set_ylabel('Eigenvalue')
fig.savefig(PROJ_DIR+'gradients/figures/lambdas_'+
'_'.join((space,parcellation,session,matrix,approach,kernel))+
'_threshold'+str(int(threshold*100))+
'_grads'+str(n_grads)+
'.png',bbox_inches="tight")

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

    # Illustrate results
    fig, ax = basic_one_plot(r,p_SA_removed,[])
    ax.set_title(sample.capitalize())
    fig.savefig(PROJ_DIR+'gradients/figures/corr_'+
    '_'.join((sample,space,parcellation,session,matrix,approach,kernel))+
    '_threshold'+str(int(threshold*100))+
    '_surrogates'+str(n_surr)+
    '_grads'+str(n_grads)+
    '.png',bbox_inches="tight")
