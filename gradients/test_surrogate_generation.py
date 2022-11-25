import numpy as np
from scipy.stats import spearmanr

# Settings
n_surr = int(10**3)
ci = 95
PROJ_DIR = '/home/cnorman/Documents/CUBRIC/ALSPAC/CardiffFC/'
sc_isv_fil = PROJ_DIR + 'data/structural/isv_sc_z_noarea.csv'
sc_isv_surr_fil = PROJ_DIR + 'data/structural/isv_sc_z_noarea_'+str(n_surr)+'_surrogates.csv'

# Load data
sc_isv = np.loadtxt(open(sc_isv_fil,"r"),delimiter=' ')
n_reg = np.shape(sc_isv)[0]
sc_isv_surr = np.loadtxt(open(sc_isv_surr_fil,"r"),delimiter=',')

# Check 'auto'correlations with surrogate maps
corrs_r_auto = np.zeros(n_surr)
corrs_p_auto = np.zeros(n_surr)
for surr in range(0,n_surr):
    corrs_r_auto[surr], corrs_p_auto[surr] = spearmanr(sc_isv,sc_isv_surr[:,surr],nan_policy='omit')
corrs_r_auto_mean = np.mean(corrs_r_auto,axis=0)
corrs_r_auto_lower = np.percentile(corrs_r_auto,(100-ci)/2,axis=0)
corrs_r_auto_upper = np.percentile(corrs_r_auto,(100+ci)/2,axis=0)
corrs_r_auto_bounds = np.row_stack((corrs_r_auto_mean-corrs_r_auto_lower,corrs_r_auto_upper-corrs_r_auto_mean))
corrs_p_auto_mean = np.mean(corrs_p_auto,axis=0)
corrs_p_auto_lower = np.percentile(corrs_p_auto,(100-ci)/2,axis=0)
corrs_p_auto_upper = np.percentile(corrs_p_auto,(100+ci)/2,axis=0)
corrs_p_auto_bounds = np.row_stack((corrs_p_auto_mean-corrs_p_auto_lower,corrs_p_auto_upper-corrs_p_auto_mean))
print('Autocorrelation rho: {} ({},{})'.format(corrs_r_auto_mean,corrs_r_auto_lower,corrs_r_auto_upper))
print('Autocorrelation p: {} ({},{})'.format(corrs_p_auto_mean,corrs_p_auto_lower,corrs_p_auto_upper))
