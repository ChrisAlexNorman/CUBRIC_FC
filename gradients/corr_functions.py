import numpy as np
from scipy.stats import spearmanr, norm
import matplotlib.pyplot as plt
from brainsmash.mapgen.base import Base
from visualisation import plot_correlations_ci

def corr_single(sc_isv,fc_grads,n_grads=5):
    """Calculate spearman correlation coeff. of given SC ISV and functional gradients"""

    corrs_r = [0] * n_grads
    corrs_p = [0] * n_grads
    for grad in range(0,n_grads):
        corrs_r[grad], corrs_p[grad] = spearmanr(sc_isv, fc_grads[:,grad], nan_policy='omit')
    
    return corrs_r, corrs_p


def corr_surrogates(sc_isv_surr,fc_grads,n_grads=5):
    """Calculate spearman correlation coeff. of given SC ISV surrogates and functional gradients"""

    n_surr = sc_isv_surr.shape[1]

    corrs_r_surr = np.zeros((n_surr,n_grads))
    corrs_p_surr = np.zeros((n_surr,n_grads))
    for surr in range(0,n_surr):
        corrs_r_surr[surr,:], corrs_p_surr[surr,:] = corr_single(sc_isv_surr[:,surr],fc_grads,n_grads=n_grads)
    
    return corrs_r_surr, corrs_p_surr


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


def sig_test(val,dist,test='nonpar'):
    """Return p value against test distribution"""

    if test == 'normal':
        # Assume data is normally distributed
        r_surr_mean = np.mean(dist)
        r_surr_var  = np.var(dist)
        # Two-tailed test against null normal distribution
        if val < r_surr_mean:
            p = 2 * norm(loc=r_surr_mean,scale=r_surr_var**0.5).cdf(val)
        else:
            p = 2 * (1-norm(loc=r_surr_mean,scale=r_surr_var**0.5).cdf(val))
        return p

    elif test == 'nonpar':
        return np.sum(np.abs(dist) > abs(val)) / len(dist)

    else:
        raise Exception("Statistical test '"+test+"' not recognized.")


def sig_test_surrogates(r, r_surr,test='nonpar'):
    """Return p values for each gradient"""

    p = [None] * len(r)

    for grad in range(len(r)):
        p[grad] = sig_test(r[grad],r_surr[:,grad],test=test)

    return p


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
        fig = plot_correlations_ci(corrs_r,approaches=approaches,kernels=kernels,thresholds=thresholds)

        plt.savefig(save_plot)

    return corrs_r, corrs_p


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
