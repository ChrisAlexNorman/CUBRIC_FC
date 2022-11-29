import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr


def basic_one_plot(corr,p=None,corr_surr=None,ci=95,sig=0.01):
    """Template function to be replaced with appropriate figure later"""

    n_grads = len(corr)

    fig = plt.figure(figsize=(4,2))
    ax = plt.axes()

    # Add raw correlation values for each grad
    ax.plot(range(1,n_grads+1),corr,'.',color='k',ms=10)

    # Add errorbars (confidence intervals) from surrogates
    if np.shape(corr_surr):
        corr_surr_mean = np.mean(corr_surr,axis=0)
        corr_surr_lower = np.percentile(corr_surr,(100-ci)/2,axis=0)
        corr_surr_upper = np.percentile(corr_surr,(100+ci)/2,axis=0)
        corr_surr_bounds = np.row_stack((corr_surr_mean-corr_surr_lower,corr_surr_upper-corr_surr_mean))
        ax.errorbar(np.asarray(range(1,n_grads+1)),corr,corr_surr_bounds,ls='none',color='k')

    # Colour statistically significant results red
    if np.shape(p):
        sig_idx = [inx for inx in range(len(p)) if p[inx] < sig]
        ax.plot([idx+1 for idx in sig_idx],corr[sig_idx],'.',color='r',ms=10)
        if np.shape(corr_surr):
            ax.errorbar([idx+1 for idx in sig_idx],corr[sig_idx],corr_surr_bounds[:,sig_idx],ls='none',color='r')

    ax.plot([1,n_grads],[0,0],'k--')
    ax.set_xticks(range(1,n_grads))
    ax.set_xlabel('Component')
    ax.set_ylabel('Corr. Coeff.')

    return fig, ax

