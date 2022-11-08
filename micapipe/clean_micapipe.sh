#!/bin/bash
###########
# Author: Chris Norman (NormanC4@cardiff.ac.uk)
##########
# SETTINGS
BIDS_DIR=/scratch/scw1648/proj_cn/data/NIfTI_Data
SAVE_DIR=/scratch/scw1648/proj_cn/data/mica_processed
MICA_IMG=~/Software/Singularity_images/micapipe-v0.1.2.simg
TEMP_DIR=/scratch/scw1648/proj_cn/pipelines/micapipe/pipeline_struct/logs
FREE_LIC=/apps/medical/freesurfer/6.0/license.txt

module load singularity
module load freesurfer/6.0
ID=1153
THREADS=20

######
# MAIN
singularity exec --cleanenv \
    -B $FREE_LIC:/opt/freesurfer-6.0.0/license.txt \
    $MICA_IMG \
    micapipe_cleanup \
        -threads $THREADS \
        -bids $BIDS_DIR \
        -out $SAVE_DIR \
        -tmpDir $TEMP_DIR \
        -sub $ID \
        -proc_structural \
        -proc_freesurfer \
        -post_structural \
	-proc_dwi \
	-SC \
	-proc_rsfmri \
	-MPC \
        -GD \

