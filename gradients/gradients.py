import numpy as np
import csv
from brainspace.gradient import GradientMaps
import argparse


def get_subj_IDs(indi_file,session='rest'):
    """Get IDs of subjects with processed FC data"""

    ID_indicators = list(csv.reader(open(indi_file,"r"), delimiter = ','))
    indicator = ID_indicators[0].index(session)
    IDs = []
    for subj in range(1,np.shape(ID_indicators)[0]):
        if int(ID_indicators[subj][indicator]):
            IDs.append(ID_indicators[subj][0])
    # Remove problem subject
    IDs.remove('1153')

    return IDs


def loop_grad_methods(subj_dir,indi_file,out_dir,space):
    """ Save functional gradients from a range of methods"""

    # Iterable methods
    approaches = ['pca','dm','le']
    kernels = ['cosine','normalized_angle','gaussian']
    thresholds = [0.75,0.85,0.9,0.95]

    # Get subject IDs
    IDs = get_subj_IDs(indi_file)

    for approach in approaches:
        for kernel in kernels:
            for threshold in thresholds:

                # Calculate gradients
                gm, settings = get_gradients(subj_dir,IDs,space=space,approach=approach,kernel=kernel,threshold=threshold)

                # Save results
                save_path = out_dir + settings["space"] + '_' + settings["atlas"] + '_' + settings["session"] + '_' + settings["matrix"] + '_' + settings["approach"] + '_' + settings["kernel"] + '_' + settings["alignment"] + '_' + str(int(settings["threshold"]*100))
                np.save(save_path + '_lambdas.npy', np.stack(gm.lambdas_))
                np.save(save_path + '_grads.npy', np.stack(gm.gradients_))
                np.save(save_path + '_aligned-grads.npy', np.stack(gm.aligned_))


def single_grad_method(subj_dir,indi_file,out_dir,space):
    """ Save functional gradients from default method"""

    # Get subject IDs
    IDs = get_subj_IDs(indi_file)

    # Calculate gradients
    gm, settings = get_gradients(subj_dir,IDs,space=space)

    # Save results
    save_path = out_dir + settings["space"] + '_' + settings["atlas"] + '_' + settings["session"] + '_' + settings["matrix"] + '_' + settings["approach"] + '_' + settings["kernel"] + '_' + settings["alignment"] + '_' + str(int(settings["threshold"]*100))
    np.save(save_path + '_lambdas.npy', np.stack(gm.lambdas_))
    np.save(save_path + '_grads.npy', np.stack(gm.gradients_))
    np.save(save_path + '_aligned-grads.npy', np.stack(gm.aligned_))
                

def get_gradients(subj_dir, IDs, session='rest', space='conte69-32k', atlas='glasser-360', matrix='FC', n_components=10, approach='dm', kernel='normalized_angle', alignment='procrustes', random_state=0, threshold=0.9):
    """Compute gradients from group data.

    Parameters
    ----------
    subj_dir : str
        Path to subject data parent directory
    IDs : list of strs
        Subject IDs to include in gradient calculation
    session : {'rest' or 'movie'}, optional
        Recording session type, default is 'rest'
    space : {'conte69-32k', 'fsaverage5', 'fsnative'}, optional
        Surface mesh name, default is 'conte69-32k'.
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
    gm : object
        GradientMaps object with fit gradients.
    settings : dict
        Dictionary of input variables used to fit gradients.
    """

    settings = {
        "session":session,
        "space":space,
        "atlas":atlas,
        "matrix":matrix,
        "n_components":n_components,
        "approach":approach,
        "kernel":kernel,
        "alignment":alignment,
        "random_state":random_state,
        "threshold":threshold
        }
    
    # Load FC data
    FC_matrices = [None] * len(IDs)
    for subj in range(0,len(IDs)):
        subj_FC = subj_dir + 'sub-' + IDs[subj] + '/func/' + session + '/surfaces/sub-' + IDs[subj] + '_rsfmri_space-' + space + '_atlas-' + atlas + '_desc-' + matrix + '.txt'
        FC = np.loadtxt(open(subj_FC,"r"), delimiter = ' ')
        # Remove non-cortical regions
        FC = FC[49:,49:]
        # Reflect upper triangle onto lower triangle
        FC = FC + FC.T - np.diag(np.diag(FC))
        FC_matrices[subj] = FC
    
    # Calculate gradients
    gm = GradientMaps(n_components = n_components,
        approach = approach,
        kernel = kernel,
        alignment = alignment,
        random_state = random_state)
    gm.fit(FC_matrices,sparsity = threshold)

    return gm, settings


if __name__ == "__main__":

    PROJ_DIR="/scratch/scw1648/proj_cn/"

    parser = argparse.ArgumentParser()
    parser.add_argument("-d","--subj_dir",default=PROJ_DIR+"data/mica_processed/micapipe/",type=str,help="Subject data directory path")
    parser.add_argument("-i","--indi_file",default=PROJ_DIR+"data/FileIndicators.csv",type=str,help="Subject file indicator reference")
    parser.add_argument("-o","--out_dir",default=PROJ_DIR+"data/gradients/",type=str,help="Results directory")
    parser.add_argument("-m","--mode",default="single",type=str,choices=["single","loop"])
    parser.add_argument("--space",default="conte69-32k",type=str,choices=["conte69-32k","fsaverage5","fsnative"])

    args = parser.parse_args()

    if args.mode=="single":
        single_grad_method(args.subj_dir,args.indi_file,args.out_dir,args.space)
    
    elif args.mode=="loop":
        loop_grad_methods(args.subj_dir,args.indi_file,args.out_dir,args.space)
