import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import seaborn as sns


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


def plot_distribution(test,dist):
    """Plot test statistic against null distribution from SA-preserving surrogate maps"""

    sac = '#377eb8'  # autocorr-preserving
    # rc = '#e41a1c'  # randomly shuffled
    bins = np.linspace(-1, 1, 51)  # correlation b

    fig = plt.figure(figsize=(3, 3))
    ax = fig.add_axes([0.2, 0.25, 0.6, 0.6])  # autocorr preserving
    # ax2 = ax.twinx()  # randomly shuffled

    # plot the data
    ax.axvline(test, 0, 0.8, color='k', linestyle='dashed', lw=1)
    ax.hist(dist, bins=bins, color=sac, alpha=1, density=True, clip_on=False, zorder=1)
    # ax2.hist(naive_brainmap_corrs, bins=bins, color=rc, alpha=0.7,density=True, clip_on=False, zorder=2)

    # make the plot nice...
    ax.set_xticks(np.arange(-1, 1.1, 0.5))
    ax.spines['left'].set_color(sac)
    ax.tick_params(axis='y', colors=sac)
    # ax2.spines['right'].set_color(rc)
    # ax2.tick_params(axis='y', colors=rc)
    ax.set_ylim(0, 2)
    # ax2.set_ylim(0, 6)
    ax.set_xlim(-1, 1)
    [s.set_visible(False) for s in [
        ax.spines['top'], ax.spines['right']]] # , ax2.spines['top'], ax2.spines['left']]]
    # ax.text(0.97, 1.1, 'SA-independent', ha='right',va='bottom', color=rc, transform=ax.transAxes)
    ax.text(0.97, 1.03, 'SA-preserving', ha='right', va='bottom',
        color=sac, transform=ax.transAxes)
    # ax.text(test, 1.65, "T1w/T2w\nmap", ha='center', va='bottom')
    ax.text(0.5, -0.2, "Corr. Coeff.",
        ha='center', va='top', transform=ax.transAxes)
    ax.text(-0.3, 0.5, "Density", rotation=90, ha='left', va='center', transform=ax.transAxes)

    return fig, ax


def plot_methods_grid(r,p,approaches=None,kernels=None,thresholds=None):
    """Plot mean correlation coefficients, with * for significance, for a range of methods and thresholding"""

    n_approaches = np.shape(r)[0]
    n_kernels = np.shape(r)[1]
    n_thresholds = np.shape(r)[2]
    n_grads = np.shape(r)[3]

    colours = ['#007991','#6BA368','#FF784F','#C33C54']
    perturb = np.linspace(-0.2,0.2,n_thresholds)

    fig, ax_all = plt.subplots(n_approaches,n_kernels, figsize=(6, 6))

    for idx_0 in range(0,n_approaches):
        for idx_1 in range(0,n_kernels):
            ax = ax_all[idx_0,idx_1]
            for idx_2 in range(0,n_thresholds):

                ax.plot(np.asarray(range(1,n_grads+1))+perturb[idx_2],r[idx_0,idx_1,idx_2],ls='none',marker='.',color=colours[idx_2],ms=10)

    for idx_0 in range(0,n_approaches):
        for idx_1 in range(0,n_kernels):
            ax = ax_all[idx_0,idx_1]
            for idx_2 in range(0,n_thresholds):

                idx_crit = np.asarray([idx_grad for idx_grad in range(n_grads) if p[idx_0,idx_1,idx_2,idx_grad] < 0.05])
                if np.size(idx_crit) > 0:
                    ax.plot(idx_crit+1+perturb[idx_2],r[idx_0,idx_1,idx_2,idx_crit],ls='none',marker='x',color='r',ms=5)
                

                # for idx_grad in range(0,n_grads):
                #     if p[idx_0,idx_1,idx_2,idx_grad] < 0.05:
                #         marker = 'x'
                #     else:
                #         marker = '.'
                #     ax.plot(idx_grad+1+perturb[idx_2],r[idx_0,idx_1,idx_2,idx_grad],ls='none',marker=marker,color=colours[idx_2])
    
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
        ax_all[n_approaches-1,idx_1].set_xticks(range(1,n_grads+1))
        ax_all[n_approaches-1,idx_1].set(xlabel='Component')
        ax_all[0,idx_1].set(title=kernels[idx_1])
        for idx_0 in range(0,n_approaches-1):
            ax_all[idx_0,idx_1].set_xticklabels([])
            ax_all[idx_0,idx_1].set_xticks([])

    for idx_0 in range(0,n_approaches):
        for idx_1 in range(0,n_kernels):
            ax_all[idx_0,idx_1].plot([1,n_grads],[0,0],'k--',linewidth=2)

    markers = [None] * n_thresholds
    for idx_2 in range(0,n_thresholds):
        markers[idx_2] = mlines.Line2D([],[],ls='none',marker='.',color=colours[idx_2],ms=10)
    ax_all[0,0].legend(markers,thresholds)

    return fig, ax


###### ARCHIVE ######

def plot_correlations_ci(corrs,approaches=None,kernels=None,thresholds=None,ci=95):
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
