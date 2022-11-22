import numpy as np
import nibabel as nb

# Label filepaths
labels_dir = './HCP_MMP1/'
label_path_micapipe = labels_dir + 'micapipe/glasser-360_conte69.csv'
label_paths_micapipe_2 = [labels_dir + 'micapipe/glasser-360_conte69_lh.label.gii', labels_dir + 'micapipe/glasser-360_conte69_rh.label.gii']
label_paths_scw = [labels_dir + 'hawk/HCP_MMP0.L.32k_fs_LR.label.gii', labels_dir + 'hawk/HCP_MMP0.R.32k_fs_LR.label.gii']
label_paths_balsa = [labels_dir + 'balsa/Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Final_Final_Areas_Group.32k_fs_LR.label.gii', labels_dir + 'balsa/Q1-Q6_RelatedParcellation210.R.CorticalAreas_dil_Final_Final_Areas_Group.32k_fs_LR.label.gii']
label_paths_balsa_2 = [labels_dir + 'balsa/Q1-Q6_RelatedValidation210.L.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.label.gii', labels_dir + 'balsa/Q1-Q6_RelatedValidation210.R.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.label.gii']
label_paths_balsa_3 = [labels_dir + 'balsa/Q1-Q6_RelatedValidation210.L.CorticalAreas_dil_210P_Orig.32k_fs_LR.label.gii', labels_dir + 'balsa/Q1-Q6_RelatedValidation210.R.CorticalAreas_dil_210P_Orig.32k_fs_LR.label.gii']

# Load labelling data
labels_micapipe = np.array(np.loadtxt(label_path_micapipe)).astype('int')
print('\nmicapipe labelling loaded from:\n' + label_path_micapipe)

labels_micapipe_2_LR = [None] * 2
labels_scw_LR = [None] * 2
labels_balsa_LR = [None] * 2
labels_balsa_2_LR = [None] * 2
labels_balsa_3_LR = [None] * 2
for h in range(2):
    labels_micapipe_2_LR[h] = np.array(nb.load(label_paths_micapipe_2[h]).darrays[0].data).astype('int')
    labels_scw_LR[h] = np.array(nb.load(label_paths_scw[h]).darrays[0].data).astype('int')
    labels_balsa_LR[h] = np.array(nb.load(label_paths_balsa[h]).darrays[0].data).astype('int')
    labels_balsa_2_LR[h] = np.array(nb.load(label_paths_balsa_2[h]).darrays[0].data).astype('int')
    labels_balsa_3_LR[h] = np.array(nb.load(label_paths_balsa_3[h]).darrays[0].data).astype('int')

labels_micapipe_2 = np.concatenate((labels_micapipe_2_LR[0],labels_micapipe_2_LR[1]))
print('\nscw server labelling loaded from:\n' + label_paths_micapipe_2[0] + '\n' + label_paths_micapipe_2[1])
if np.max(labels_micapipe_2) < 360:
    print('NOTE: L and R hemispheres encode the same values')

labels_scw = np.concatenate((labels_scw_LR[0],labels_scw_LR[1]))
print('\nscw server labelling loaded from:\n' + label_paths_scw[0] + '\n' + label_paths_scw[1])
if np.max(labels_scw) < 360:
    print('NOTE: L and R hemispheres encode the same values')

labels_balsa = np.concatenate((labels_balsa_LR[0],labels_balsa_LR[1]))
print('\nbalsa labelling loaded from:\n' + label_paths_balsa[0] + '\n' + label_paths_balsa[1])
if np.max(labels_balsa) < 360:
    print('NOTE: L and R hemispheres encode the same values')

labels_balsa_2 = np.concatenate((labels_balsa_2_LR[0],labels_balsa_2_LR[1]))
print('\nbalsa 2 labelling loaded from:\n' + label_paths_balsa_2[0] + '\n' + label_paths_balsa_2[1])
if np.max(labels_balsa_2) < 360:
    print('NOTE: L and R hemispheres encode the same values')

labels_balsa_3 = np.concatenate((labels_balsa_3_LR[0],labels_balsa_3_LR[1]))
print('\nbalsa 2 labelling loaded from:\n' + label_paths_balsa_3[0] + '\n' + label_paths_balsa_3[1])
if np.max(labels_balsa_3) < 360:
    print('NOTE: L and R hemispheres encode the same values')

# Compare background / medial wall indicies
micapipe_0s = (labels_micapipe == 0)
micapipe_2_0s = (labels_micapipe_2 == 0)
scw_0s = (labels_scw == 0)
balsa_0s = (labels_balsa == 0)
balsa_2_0s = (labels_balsa_2 == 0)
balsa_3_0s = (labels_balsa_3 == 0)

micapipe_micapipe_2_0s_eq = (micapipe_0s == micapipe_2_0s)
if np.all(micapipe_micapipe_2_0s_eq):
    print('\nmicapipe and micapipe 2 background / medial wall vertices match across both hemispheres')
else:
    print('\nmicapipe and micapipe 2 background / medial wall vertices do not match')

micapipe_scw_0s_eq = (micapipe_0s == scw_0s)
if np.all(micapipe_scw_0s_eq):
    print('\nscw and micapipe background / medial wall vertices match across both hemispheres')
else:
    print('\nscw and micapipe background / medial wall vertices do not match')

micapipe_balsa_0s_eq = (micapipe_0s == balsa_0s)
if np.all(micapipe_balsa_0s_eq):
    print('balsa and micapipe background / medial wall vertices match across both hemispheres')
else:
    print('balsa and micapipe background / medial wall vertices do not match')

micapipe_balsa_2_0s_eq = (micapipe_0s == balsa_2_0s)
if np.all(micapipe_balsa_2_0s_eq):
    print('balsa 2 and micapipe background / medial wall vertices match across both hemispheres')
else:
    print('balsa 2 and micapipe background / medial wall vertices do not match')

micapipe_balsa_3_0s_eq = (micapipe_0s == balsa_3_0s)
if np.all(micapipe_balsa_3_0s_eq):
    print('balsa 3 and micapipe background / medial wall vertices match across both hemispheres')
else:
    print('balsa 3 and micapipe background / medial wall vertices do not match')

# Check if all labellings match micapipe's
if np.max(labels_micapipe_2) < 360:
    labels_micapipe_2_LR[1][labels_micapipe_2_LR[1]>0] = labels_micapipe_2_LR[1][labels_micapipe_2_LR[1]>0] + labels_micapipe_2_LR[0].max()
    labels_micapipe_2 = np.concatenate((labels_micapipe_2_LR[0],labels_micapipe_2_LR[1]))
micapipe_micapipe_2_full_eq = (labels_micapipe == labels_micapipe_2)
if np.all(micapipe_micapipe_2_full_eq):
    print('\nAll other micapipe and micapipe 2 labellings match too!')
else:
    print('\nNot all other micapipe and micapipe 2 labellings match though')

if np.max(labels_scw) < 360:
    labels_scw_LR[1][labels_scw_LR[1]>0] = labels_scw_LR[1][labels_scw_LR[1]>0] + labels_scw_LR[0].max()
    labels_scw = np.concatenate((labels_scw_LR[0],labels_scw_LR[1]))
micapipe_scw_full_eq = (labels_micapipe == labels_scw)
if np.all(micapipe_scw_full_eq):
    print('\nAll other scw and micapipe labellings match too!')
else:
    print('\nNot all other scw and micapipe labellings match though')

if np.max(labels_balsa) < 360:
    labels_balsa_LR[1][labels_balsa_LR[1]>0] = labels_balsa_LR[1][labels_balsa_LR[1]>0] + labels_balsa_LR[0].max()
    labels_balsa = np.concatenate((labels_balsa_LR[0],labels_balsa_LR[1]))
micapipe_balsa_full_eq = (labels_micapipe == labels_balsa)
if np.all(micapipe_balsa_full_eq):
    print('All other balsa and micapipe labellings match too!')
else:
    print('Not all other balsa and micapipe labellings match though')

if np.max(labels_balsa_2) < 360:
    labels_balsa_2_LR[1][labels_balsa_2_LR[1]>0] = labels_balsa_2_LR[1][labels_balsa_2_LR[1]>0] + labels_balsa_2_LR[0].max()
    labels_balsa_2 = np.concatenate((labels_balsa_2_LR[0],labels_balsa_2_LR[1]))
micapipe_balsa_2_full_eq = (labels_micapipe == labels_balsa_2)
if np.all(micapipe_balsa_2_full_eq):
    print('All other balsa 2 and micapipe labellings match too!')
else:
    print('Not all other balsa 2 and micapipe labellings match though')

if np.max(labels_balsa_3) < 360:
    labels_balsa_3_LR[1][labels_balsa_3_LR[1]>0] = labels_balsa_3_LR[1][labels_balsa_3_LR[1]>0] + labels_balsa_3_LR[0].max()
    labels_balsa_3 = np.concatenate((labels_balsa_3_LR[0],labels_balsa_3_LR[1]))
micapipe_balsa_3_full_eq = (labels_micapipe == labels_balsa_3)
if np.all(micapipe_balsa_3_full_eq):
    print('All other balsa 3 and micapipe labellings match too!')
else:
    print('Not all other balsa 3 and micapipe labellings match though')

# Check if other labellings match each other
balsa_scw_full_eq = (labels_balsa == labels_scw)
if np.all(balsa_scw_full_eq):
    print('All other balsa and scw labellings match too!')
else:
    print('Not all other balsa and scw labellings match though')

balsa_2_scw_full_eq = (labels_balsa_2 == labels_scw)
if np.all(balsa_2_scw_full_eq):
    print('All other balsa 2 and scw labellings match too!')
else:
    print('Not all other balsa 2 and scw labellings match though')

balsa_balsa_2_full_eq = (labels_balsa == labels_balsa_2)
if np.all(balsa_balsa_2_full_eq):
    print('All other balsa and balsa 2 labellings match too!')
else:
    print('Not all other balsa and balsa 2 labellings match though')

balsa_balsa_3_full_eq = (labels_balsa == labels_balsa_3)
if np.all(balsa_balsa_3_full_eq):
    print('All other balsa and balsa 3 labellings match too!')
else:
    print('Not all other balsa and balsa 3 labellings match though')



