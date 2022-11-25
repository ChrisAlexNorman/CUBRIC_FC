import numpy as np
from correlations import make_surrogates
from brainsmash.mapgen.eval import base_fit
import matplotlib.pyplot as plt

# Settings
plt.ioff()
n_jobs = 10
n_surrogates = 10**4
PROJ_DIR = '/home/cnorman/Documents/CUBRIC/ALSPAC/CardiffFC/'
sc_isv_fil = PROJ_DIR + 'data/structural/isv_sc_z_noarea.csv'
out_fil_surr = PROJ_DIR+'data/structural/isv_sc_z_noarea_'+str(n_surrogates)+'_surrogates.csv'
out_fil_lh_vario = PROJ_DIR+'gradients/Figures/sc_isv_lh_'+str(n_surrogates)+'_surrogates_variogram.png'
out_fil_rh_vario = PROJ_DIR+'gradients/Figures/sc_isv_rh_'+str(n_surrogates)+'_surrogates_variogram.png'

# Load SC ISV data
sc_isv_raw = np.loadtxt(open(sc_isv_fil,"r"),delimiter=' ')
nreg = np.shape(sc_isv_raw)[0]
sc_isv = np.row_stack((sc_isv_raw[:int(len(sc_isv_raw)/2)],sc_isv_raw[int(len(sc_isv_raw)/2):]))

# Load geodesic distances
dist_mats = np.zeros((2,int(nreg/2),int(nreg/2)))
dist_mats[0,:,:] = np.loadtxt(PROJ_DIR + 'parcellations/HCP-MMP1/lh.fsaverage-32k_glasser-360_geodesic_distmat.txt')
dist_mats[1,:,:] = np.loadtxt(PROJ_DIR + 'parcellations/HCP-MMP1/rh.fsaverage-32k_glasser-360_geodesic_distmat.txt')

# Make & save surrogates
lh_surrogates = make_surrogates(sc_isv[0], dist_mats[0], n=n_surrogates, n_jobs=n_jobs)
rh_surrogates = make_surrogates(sc_isv[1], dist_mats[1], n=n_surrogates, n_jobs=n_jobs)
surrogates = np.concatenate((lh_surrogates,rh_surrogates),axis=0)
np.savetxt(out_fil_surr,surrogates,delimiter=',')

# Make & save variogram figures
# Left
base_fit(sc_isv[0], dist_mats[0],nsurr=n_surrogates)
plt.savefig(out_fil_lh_vario)
# Right
base_fit(sc_isv[1], dist_mats[1],nsurr=n_surrogates)
plt.savefig(out_fil_rh_vario)

