#!/bin/bash
###########
# Author: Chris Norman (NormanC4@cardiff.ac.uk)
##########
# SETTINGS
BIDS_DIR=/scratch/scw1648/proj_cn/data/NIfTIs
SAVE_DIR=/scratch/scw1648/proj_cn/data/OUT_DIR
MICA_IMG=~/Software/Singularity_images/micapipe-v0.1.2.simg
TEMP_DIR=/scratch/scw1648/proj_cn/micapipe/pipeline_struct/logs
FREE_LIC=/apps/medical/freesurfer/6.0/license.txt

module load singularity
module load freesurfer/6.0
ID=1153
THREADS=20

######
# MAIN
singularity run --cleanenv \
    -B $FREE_LIC:/opt/freesurfer-6.0.0/license.txt \
    $MICA_IMG \
        -threads $THREADS \
        -bids $BIDS_DIR \
        -out $SAVE_DIR \
        -tmpDir $TEMP_DIR \
        -sub $ID \
        -proc_structural \
        -proc_freesurfer \
        -post_structural \
        -GD \
        -Morphology

