# Generate selected figures for SC ISV and FC gradients

import numpy as np
# import argparse
# from correlations import get_mean_grads, make_surrogates
from scipy.stats import spearmanr
import matplotlib.pyplot as plt
from brainsmash.mapgen.base import Base

def corr_grad_scisv(sc_isv,fc_grads,n_grads=5):
    """Calculate spearman correlation coeff. of given SC ISV and functional gradients"""

    corrs_r = np.zeros((n_grads))
    corrs_p = np.zeros((n_grads))
    for grad in range(0,n_grads):
        corrs_r[grad], corrs_p[grad] = spearmanr(sc_isv, fc_grads[:,grad], nan_policy='omit')
    
    return corrs_r, corrs_p

def corr_grad_surrogates_loop_methods(sc_isv,dist_mats,grad_dir,approaches,kernels,thresholds,save_plot=False,space='fsaverage5',parcellation='glasser-360',session='rest',matrix='FC',alignment='procrustes',n_grads=5,n_surrogates=1000):
    """Calculate spearman correlation coeff. of surrogate SC ISV maps and functional gradients calculated from different methods
    
    Parameters
    ----------
    sc_isv : (2,N) np.ndarray
        Structural connectivity inter-subject variability as (lh, rh)
    dist_mats (2,N,N) np.ndarray
        Distance matrices for specified parcellation as (lh,rh)
    grad_dir : str
        Path to directory containing gradient data
    approaches : list of strs
        Embedding approaches to loop over, possible options:
        'dm' or DiffusionMaps: embedding using diffusion maps.
        'le' or LaplacianEigenmaps: embedding using Laplacian eigenmaps.
        'pca' or PCAMaps: embedding using PCA.
    kernels : list of strs
        Kernel function used to build the affinity matrix. Possible options: {'cosine', 'normalized_angle', 'gaussian'}.
    thresholds : list of strs
        Thresholding percentage values (as strings)
    save_plot : str or bool, optional
        Path to output figure or False for no figure. Default is False
    space : str, optional
        Surface mesh name, default is 'fsaverage5'.
    parcellation : str, optional
        Parcellation name, default is 'glasser-360'
    session : {'rest' or 'movie'}, optional
        Session type, default is 'rest'
    matrix : {'FC' or 'timeseries'}, optional
        Matrix type, default is 'FC'
    alignment : {'procrustes' or 'joint'}, optional
        Approach used to align gradients across subjects.
    n_grads : int, optional
        Gradients number to plot up to, default is 5.
    n_surrogates : int, optional
        Number of surrogate SC ISV maps to generate, default is 100

    Returns
    -------
    corrs_r : (N_approaches,N_kernels,N_thresholds,N_surrogates,N_grads) np.ndarray
        Correlation coefficients
    corrs_p : (N_approaches,N_kernels,N_thresholds,N_surrogates,N_grads) np.ndarray
        Correlation p-values
    """

    # Make SC ISV surrogate maps
    sc_isv_lh_surrogates = make_surrogates(sc_isv[0], dist_mats[0], n=n_surrogates, n_jobs=10)
    sc_isv_rh_surrogates = make_surrogates(sc_isv[1], dist_mats[1], n=n_surrogates, n_jobs=10)
    sc_isv_surrogates = np.concatenate((sc_isv_lh_surrogates,sc_isv_rh_surrogates),axis=0)

    # Loop over methods
    corrs_r = np.zeros((len(approaches),len(kernels),len(thresholds),n_surrogates,n_grads))
    corrs_p = np.zeros((len(approaches),len(kernels),len(thresholds),n_surrogates,n_grads))
    for idx_0 in range(0,len(approaches)):
        approach = approaches[idx_0]
        for idx_1 in range(0,len(kernels)):
            kernel = kernels[idx_1]
            for idx_2 in range(0,len(thresholds)):
                threshold = thresholds[idx_2]

                grad_fil = grad_dir + space+'_'+parcellation+'_'+session+'_'+matrix+'_'+approach+'_'+kernel+'_'+alignment+'_'+threshold+'_aligned-grads.npy'
                
                mean_grads = get_mean_grads(grad_fil,standardize=True)

                for idx_3 in range(0,n_surrogates):
                    for idx_4 in range(0,n_grads):
                        [corrs_r[idx_0,idx_1,idx_2,idx_3,idx_4], corrs_p[idx_0,idx_1,idx_2,idx_3,idx_4]] = spearmanr(sc_isv_surrogates[:,idx_3], mean_grads[:,idx_4], nan_policy='omit')
    
    if save_plot:
        fig = plot_correlations(corrs_r,approaches=approaches,kernels=kernels,thresholds=thresholds)

        plt.savefig(save_plot)

    return corrs_r, corrs_p

def plot_correlations(corrs,approaches=None,kernels=None,thresholds=None,ci=95):
    """Plot mean correlation coefficients with 95% CIs from bootstrapping against surrogate maps for a range of methods and thresholding"""

    n_approaches = np.shape(corrs)[0]
    n_kernels = np.shape(corrs)[1]
    n_thresholds = np.shape(corrs)[2]
    n_surrogates = np.shape(corrs)[3]
    n_grads = np.shape(corrs)[4]

    colours = ['#007991','#6BA368','#FF784F','#C33C54']
    perturb = np.linspace(-0.2,0.2,n_thresholds)

    fig, ax_all = plt.subplots(n_approaches,n_kernels, figsize=(6, 6))

    for idx_0 in range(0,n_approaches):
        for idx_1 in range(0,n_kernels):
            ax = ax_all[idx_0,idx_1]
            for idx_2 in range(0,n_thresholds):

                r_means = np.mean(corrs[idx_0,idx_1,idx_2],axis=0)
                r_lower = np.percentile(corrs[idx_0,idx_1,idx_2],(100-ci)/2,axis=0)
                r_upper = np.percentile(corrs[idx_0,idx_1,idx_2],(100+ci)/2,axis=0)
                r_bounds = np.row_stack((r_means-r_lower,r_upper-r_means))

                ax.errorbar(np.asarray(range(1,n_grads+1))+perturb[idx_2],r_means,r_bounds,ls='none',color=colours[idx_2])
                ax.plot(np.asarray(range(1,n_grads+1))+perturb[idx_2],r_means,marker='.',markersize=10,ls='none',color=colours[idx_2])
    
    if approaches==None:
        approaches = [None]*n_approaches
        for n in range(0,n_approaches):
            approaches[n] = 'Approach ' + str(n+1)
    if kernels==None:
        kernels = [None]*n_kernels
        for n in range(0,n_kernels):
            kernels[n] = 'Kernel ' + str(n+1)
    if thresholds==None:
        thresholds = [None]*n_thresholds
        for n in range(0,n_thresholds):
            thresholds[n] = 'Threshold ' + str(n+1)

    for idx_0 in range(0,n_approaches):
        ax_all[idx_0,0].set(ylabel=approaches[idx_0]+"\n Corr. coef.")
        ylim_max = list(ax_all[idx_0,0].get_ylim())
        for idx_1 in range(1,n_kernels):
            ax_all[idx_0,idx_1].set_yticklabels([])
            ax_all[idx_0,idx_1].set_yticks([])
            ylim = ax_all[idx_0,idx_1].get_ylim()
            if ylim[0] < ylim_max[0]:
                ylim_max[0] = ylim[0]
            if ylim[1] > ylim_max[1]:
                ylim_max[1] = ylim[1]
        for idx_1 in range(0,n_kernels):
            ax_all[idx_0,idx_1].set_ylim(ylim_max)

    for idx_1 in range(0,n_kernels):
        ax_all[n_approaches-1,idx_1].set(xlabel='Component')
        ax_all[0,idx_1].set(title=kernels[idx_1])
        for idx_0 in range(0,n_approaches-1):
            ax_all[idx_0,idx_1].set_xticklabels([])
            ax_all[idx_0,idx_1].set_xticks([])

    for idx_0 in range(0,n_approaches):
        for idx_1 in range(0,n_kernels):
            ax_all[idx_0,idx_1].plot([1,n_grads],[0,0],'k--',linewidth=2)

    ax_all[0,0].legend(thresholds)

    return fig, ax

def plot_correlations_single(corrs,ci=95):
    """Plot mean correlation coefficients with 95% CIs from bootstrapping against surrogate maps"""

    n_surrogates = np.shape(corrs)[0]
    n_grads = np.shape(corrs)[1]

    colour = ['#000000']

    fig, ax = plt.plot()

    r_means = np.mean(corrs,axis=0)
    r_lower = np.percentile(corrs,(100-ci)/2,axis=0)
    r_upper = np.percentile(corrs,(100+ci)/2,axis=0)
    r_bounds = np.row_stack((r_means-r_lower,r_upper-r_means))

    ax.errorbar(np.asarray(range(1,n_grads+1)),r_means,r_bounds,ls='none',color=colour)
    ax.plot(np.asarray(range(1,n_grads+1)),r_means,marker='.',markersize=10,ls='none',color=colour)

    ax.set(ylabel="Corr. Coef.")
    ax.set(xlabel="Component")

    return fig, ax

def get_mean_grads(grad_fil,standardize=True):
    """Return mean (and possibly standardized) gradients across subjects"""

    grads = np.load(grad_fil)
    
    if standardize:
        for subj in range(0,np.shape(grads)[0]):
            for grad_i in range(0,np.shape(grads)[2]):
                col = grads[subj,:,grad_i]
                grads[subj,:,grad_i] = (col - np.mean(col)) / np.std(col)
    
    mean_grads = np.mean(grads,0)

    return mean_grads

def make_surrogates(map_data, dist_mat, n=100, resample=True, n_jobs=1):
    """Make surrogates maps using brainsmash to remove spatial autocorrelation"""

    base = Base(x=map_data, D=dist_mat, resample=resample, n_jobs=n_jobs)

    surrogates = np.transpose(base(n=n))

    return surrogates

if __name__ == "__main__":
    # Settings
    PROJ_DIR = '/home/cnorman/Documents/CUBRIC/ALSPAC/CardiffFC/'

    sc_isv_fil = PROJ_DIR + 'data/structural/isv_sc_z_noarea.csv'
    sc_isv_raw = np.loadtxt(sc_isv_fil)
    sc_isv = np.row_stack((sc_isv_raw[:int(len(sc_isv_raw)/2)],sc_isv_raw[int(len(sc_isv_raw)/2):]))

    dist_mats = np.zeros((2,180,180))
    dist_mats[0,:,:] = np.loadtxt(PROJ_DIR + 'parcellations/HCP-MMP1/lh.fsaverage-32k_glasser-360_geodesic_distmat.txt')
    dist_mats[1,:,:] = np.loadtxt(PROJ_DIR + 'parcellations/HCP-MMP1/rh.fsaverage-32k_glasser-360_geodesic_distmat.txt')

    space = 'fsaverage5'
    grad_dir = PROJ_DIR + 'data/gradients/' + space + '/'
    approaches = ['pca','dm','le']
    kernels = ['cosine','normalized_angle','gaussian']
    thresholds = ['75','85','90','95']

    plot_path = PROJ_DIR + 'gradients/correlations.png'

    corrs_r, corrs_p = corr_grad_surrogates_loop_methods(sc_isv,dist_mats,grad_dir,approaches,kernels,thresholds,save_plot=plot_path,n_surrogates=10000)
