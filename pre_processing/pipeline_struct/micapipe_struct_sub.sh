#!/bin/bash
###########
# Author: Chris Norman (NormanC4@cardiff.ac.uk)
##########
# SETTINGS
BIDS_DIR=/scratch/scw1648/proj_cn/data/NIfTIs
SAVE_DIR=/scratch/scw1648/proj_cn/data/mica_processed
MICA_IMG=~/Software/Singularity_images/micapipe-v0.1.2.simg
TEMP_DIR=/scratch/scw1648/proj_cn/micapipe/pipeline_struct/logs

module load singularity
module load freesurfer/6.0

############
# ID RANGE
ID_START=$1
declare -i ID_END_TEST=$ID_START+5
declare -i ID_END=$(( $ID_END_TEST < 1220 ? $ID_END_TEST : 1220 ))

######
# MAIN
for ID in $(seq $ID_START $ID_END); do
    singularity run --cleanenv \
    -B /apps/medical/freesurfer/6.0/license.txt:/opt/freesurfer-6.0.0/license.txt \
    $MICA_IMG \
        -threads 20 \
        -bids $BIDS_DIR \
        -out $SAVE_DIR \
        -tmpDir $TEMP_DIR \
        -sub $ID \
        -proc_structural \
        -proc_freesurfer \
        -post_structural \
        -GD \
        -Morphology
done

