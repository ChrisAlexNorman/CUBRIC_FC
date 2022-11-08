# This is a docstring

import os
import sys
import csv
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.stats import spearmanr
from scipy.stats import pearsonr
from brainspace.utils.parcellation import map_to_labels
from brainspace.gradient import GradientMaps
from brainspace.plotting import plot_hemispheres
from brainspace.datasets import load_conte69
from brainspace.datasets import load_parcellation
from brainspace.datasets import load_group_fc
from bootstrap_corr import bootstrap_corr
from cn_load_hcpmmp import load_hcpmmp

from brainspace.datasets import load_conte69
from wbplot import pscaler

def load_surf_parc(parc_dir,space='fsaverage5',atlas='glasser-360'):
    """load requirements for illustration across specified surface

    Parameters
    ----------
    parc_dir : str
        Path to directory containing surface parcellation files
    space : {'conte69-32k', 'fsaverage5', 'fsnative'}, optional
        Surface mesh name, default is 'fsaverage5'.
    atlas : str, optional
        Parcellation name. Check subject data for available options.
        Default is 'glasser-360'
    matrix : {'FC' or 'timeseries'}, optional
        Matrix type, default is 'FC'.
    n_components : int, optional
        Number of gradients. Default is 10.
    approach : {'dm', 'le', 'pca'} or object, optional
        Embedding approach. Default is 'dm'. It can be a string or instance:
        'dm' or DiffusionMaps: embedding using diffusion maps.
        'le' or LaplacianEigenmaps: embedding using Laplacian eigenmaps.
        'pca' or PCAMaps: embedding using PCA.
    kernel : str, callable or None, optional
        Kernel function to build the affinity matrix. Possible options: {'pearson', 'spearman', 'cosine', 'normalized_angle', 'gaussian'}.
        If callable, must receive a 2D array and return a 2D square array.
        If None, use input matrix. Default is None.
    alignment : {'procrustes', 'joint'}, object or None
        Alignment approach. Only used when two or more datasets are provided. If None, no alignment is performed. If object, it accepts an instance of ProcrustesAlignment. Default is None.
        If 'procrustes', datasets are aligned using generalized procrustes analysis.
        If 'joint', datasets are embedded simultaneously based on a joint affinity matrix built from the individual datasets. This option is only available for 'dm' and 'le' approaches.
    random_state : int or None, optional
        Random state. Default is None.
    threshold : float, optional
        Proportion of the smallest elements to zero-out for each row.
        Default is 0.9.

    Returns
    -------
    surf_lh : object
        Left hemisphere surface.
    surf_rh : object
        Right hemisphere surface.
    labelling : list
    settings : dict
        Dictionary of input variables used to fit gradients.
    """
    


# Paths
PROJ_DIR = '~/Documents/CUBRIC/ALSPAC/proj_cn/'


SC_ISV_Fil = Parent_Dir + 'isv_sc_z_noarea.csv'
atlas = 'schaefer-400'
FC_Fil = Parent_Dir + 'FC_rest_conte69-32k_' + atlas + '.txt'
HCP_Labs_Dir = Parent_Dir + 'HCP_MMP'
Fig_Dir = Parent_Dir + 'Results/Builtin_schaefer-400'

threshold = 0.9
grad_max = 4

# Turn off interactive plots
plt.ioff()

# Load and visualise pre-calculated SC_ISV
SC_ISV = np.loadtxt(SC_ISV_Fil)

surf_lh, surf_rh, labeling = load_hcpmmp(HCP_Labs_Dir)
surf_lh, surf_rh = load_conte69()
# labeling = load_parcellation('schaefer', scale=400, join=True)
mask = labeling > 0

# Load mean FC map
FC_Mat_Raw = list(csv.reader(open(FC_Fil, "r"), delimiter=','))
FC_Mat = np.array(FC_Mat_Raw).astype('float')[49:,49:]
nreg = FC_Mat.shape[0]

plt.imshow(FC_Mat,cmap='hot')
plt.colorbar()
plt.savefig(Fig_Dir+'/FC_mean_'+atlas+'.png')

# Find Gradients
n_components = 10
gm = GradientMaps(n_components=n_components,
                  approach='dm',
                  kernel='cosine',
                  random_state=0)
gm.fit(FC_Mat)
fc_grads = gm.gradients_

# Standardize gradients
fc_grads_z = np.zeros(np.shape(fc_grads))
for gi in range(0,np.shape(fc_grads)[1]):
    col = fc_grads[:,gi]
    col_z = (col - np.mean(col))/np.std(col)
    fc_grads_z[:,gi] = col_z
fc_grads = fc_grads_z

# Visualise gradients
text_labels = [None] * grad_max
for g in range(grad_max):
    text_labels[g]='Grad'+str(g+1)
grad = [None] * n_components
for i, g in enumerate(fc_grads.T):
    grad[i] = map_to_labels(g, labeling, mask=mask, fill=np.nan)
plot_hemispheres(surf_lh, surf_rh, array_name=grad[:grad_max], size=(1200, 400), cmap='bwr',
                 color_bar=True, label_text=text_labels, zoom=1.5,
                 embed_nb=False, screenshot=True, filename=Fig_Dir+'/FC_gradients_1-'+str(grad_max)+'_conte69-32k_'+atlas+'.png')

# Correlate gradient with isv
for n in range(n_components):
    [r, p] = spearmanr(SC_ISV, fc_grads[:, n], nan_policy='omit')
    print('Spearman correlation fmri gradient {} and sc_isv: {:.2f} ({:.6f})'.format(n, r, p))

print('Bootstrapping...')
rho, pval, ci = bootstrap_corr(SC_ISV, fc_grads[:, :grad_max], alpha=0.05, n_rep=10000, corr_type='spearman')
print('Done.')
#L/R separate
# rhoL, pvalL, ciL = bootstrap_corr(SC_ISV[:180], fc_grads[:180, :gmax], alpha=0.05, n_rep=10000, corr_type='spearman')
# rhoR, pvalR, ciR = bootstrap_corr(SC_ISV[180:], fc_grads[180:, :gmax], alpha=0.05, n_rep=10000, corr_type='spearman')
ci_err = np.zeros(np.shape(ci))
ci_err[:,0] = rho - ci[:,0]
ci_err[:,1] = ci[:,1] - rho

for n in range(grad_max):
    print('Spearman correlation fmri gradient {} and sc_isv: r={:.5f},p={:.6f},[{:.2f},{:.2f}]'.format(n + 1, rho[n],
                                                                                                       pval[n],
                                                                                                       ci[n, 0],
                                                                                                       ci[n, 1]))
    # print('Spearman correlation Left fmri gradient {} and sc_isv: r={:.5f},p={:.6f},[{:.2f},{:.2f}]'.format(n + 1, rhoL[n],
    #                                                                                                    pvalL[n],
    #                                                                                                    ciL[n, 0],
    #                                                                                                    ciL[n, 1]))
    # print('Spearman correlation Right fmri gradient {} and sc_isv: r={:.5f},p={:.6f},[{:.2f},{:.2f}]'.format(n + 1, rhoR[n],
    #                                                                                                    pvalR[n],
    #                                                                                                    ciR[n, 0],
    #                                                                                           ciR[n, 1]))
# scatter plots
for n in range(grad_max):
    plt.figure(figsize=(10, 7))
    sns.regplot(x=fc_grads[:, n], y=SC_ISV, 
    scatter_kws={'color':'black'}, line_kws={'color':'red'})
    # plt.xlim(-0.05,0.17)
    # plt.ylim(0, 0.27)
    plt.xlabel('Gradient ' + str(n + 1), fontsize=18)
    plt.ylabel('sc-isv', fontsize=18)
    plt.title('Structural variability vs Gradient {}'.format(n + 1), fontsize=18)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.tight_layout()
    # plt.show()
    plt.savefig(Fig_Dir+'/FC_gradient_'+str(n+1)+'_SC_ISV_correlation_'+atlas+'.png')

# Plot eigenvalues with correlation
fig, (ax1,ax2) = plt.subplots(2,1, figsize=(6, 12))
ax1.plot(range(1,1+gm.lambdas_.size), gm.lambdas_,marker='.',markersize=10)
ax1.set(title='Component variance',xlabel='Component',
        ylabel='Eigenvalue')

ax2.errorbar(range(1,grad_max+1),rho,ci_err.T)
ax2.plot([1,grad_max],[0,0],'k--',linewidth=2)
ax2.set(title='Component SC_ISV-FC_grad Correlation',xlabel='Component',
        ylabel='Correlation Coeff.')

# This was to plot cumulative eigenvalues:
# ax2.plot(range(1,1+gm.lambdas_.size), np.cumsum(gm.lambdas_),marker='.',markersize=10)
# ax2.set(title='Cumulative variance',xlabel='Component',
#         ylabel='Eigenvalue')

# plt.show()
plt.savefig(Fig_Dir+'/FC_gradient_eigenvalues_'+atlas+'.png')