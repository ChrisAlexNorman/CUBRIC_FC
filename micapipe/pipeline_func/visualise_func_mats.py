import os
import sys
sys.path.insert(0,'/nfshome/store02/users/c.sapcn2/.local/lib/python3.7/site-packages')
# sys.path.append('/home/c.sapcn2/.local/bin')
import numpy as np
from nilearn import plotting
import matplotlib as plt

results_dir = '/scratch/scw1648/proj_cn/micapipe/examples/EXAMPLE_mica_outputs'
fig_dir = '/scratch/scw1648/proj_cn/micapipe/examples'


#########################
## FUNCTIONAL CONNECTOMES

subjectIDs=['sub-1011','sub-1011','sub-1011','sub-1028']
atlases=['schaefer-200','schaefer-800','schaefer-200','schaefer-200']
sessions=['rest','rest','movie','rest']

mat_type = 'FC'
for n in range(len(subjectIDs)):
	subjectID = subjectIDs[n]
	atlas = atlases[n]
	session = sessions[n]
	
	cnt_fs = results_dir + '/micapipe/' + subjectID + '/func/' + session + '/surfaces/' + subjectID + '_rsfmri_space-fsnative_atlas-' + atlas + '_desc-' + mat_type + '.txt'
	
	# Load the connectome
	mtx_fs = np.loadtxt(cnt_fs, dtype=np.float, delimiter=' ')
	# Fill the lower triangle of the matrix
	mtx_fcSym = np.triu(mtx_fs,1)+mtx_fs.T
	# Plot the matrix
	corr_plot = plotting.plot_matrix(mtx_fcSym, figure=(10, 10), labels=None, cmap='Reds')
	# Save figure
	plt.pyplot.savefig(fig_dir+'/'+subjectID+'_'+atlas+'_'+session+'_'+mat_type+'.pdf')

mat_type = 'timeseries'
for n in range(len(subjectIDs)):
	subjectID = subjectIDs[n]
	atlas = atlases[n]
	session = sessions[n]
	
	cnt_time = results_dir + '/micapipe/' + subjectID + '/func/' + session + '/surfaces/' + subjectID + '_rsfmri_space-fsnative_atlas-' + atlas + '_desc-' + mat_type + '.txt'
	
	# Load the time series
	mtx_time = np.loadtxt(cnt_time, dtype=np.float, delimiter=' ')
	# Plot as a matrix
	corr_plot = plotting.plot_matrix(mtx_time.T, figure=(12, 5), labels=None, cmap='Reds')
	# Save Figure
	plt.pyplot.savefig(fig_dir+'/'+subjectID+'_'+atlas+'_'+session+'_'+mat_type+'.pdf')

