#!/bin/bash
###########
# Author: Chris Norman (NormanC4@cardiff.ac.uk)
#
# !!!CHANGES TO MICAPIPE DEFAULT FILESYSTEM!!!
# Subdirectories 'rest' and 'movie' created in 
# <outputDirectory>/micapipe/sub-$ID/func (initially empty)
# <outputDirectory>/micapipe/sub-$ID/xfms (contains struct files)
# New files are directed to those subdirectories instead of the default parents.
# Hopefully this doesn't cause any problems down the line...
#
# NOTES:
# freesurfer license path required
# Strangely the .json file for each scan is assumed to be $BIDS_DIR/task-rest_bold.json. This path is re-directed to subject & scan specific file. Destination must be created first.
##########
# SETTINGS
BIDS_DIR=/scratch/scw1648/proj_cn/data/NIfTIs
SAVE_DIR=/scratch/scw1648/proj_cn/data/OUT_DIR
MICA_IMG=~/Software/Singularity_images/micapipe-v0.1.2.simg
TEMP_DIR=/scratch/scw1648/proj_cn/micapipe/pipeline_func/logs
FREE_LIC=/apps/medical/freesurfer/6.0/license.txt

module load singularity
module load freesurfer/6.0
ID=1153
THREADS=20
touch "$BIDS_DIR/task-rest_bold.json"

######
# MAIN

# Check if unwarped files exist & specify scan names
REST_NAME=rest_unwarped
[[ ! -f "$BIDS_DIR/sub-$ID/func/sub-$ID""_$REST_NAME.nii" ]] && REST_NAME=rest
MOVIE_NAME=movie_unwarped
[[ ! -f "$BIDS_DIR/sub-$ID/func/sub-$ID""_$MOVIE_NAME.nii" ]] && MOVIE_NAME=movie

#### rest
# Create subdirectories
[[ ! -d "$SAVE_DIR/micapipe/sub-$ID/func/rest" ]] && mkdir "$SAVE_DIR/micapipe/sub-$ID/func/rest"
[[ ! -d "$SAVE_DIR/micapipe/sub-$ID/xfm/rest" ]] && mkdir "$SAVE_DIR/micapipe/sub-$ID/xfm/rest"

# Run
singularity run --cleanenv \
    -B $FREE_LIC:/opt/freesurfer-6.0.0/license.txt \
    -B "$BIDS_DIR/sub-$ID/func/sub-$ID""_rest.json":"$BIDS_DIR/task-rest_bold.json" \
    -B "$SAVE_DIR/micapipe/sub-$ID/func/rest":"$SAVE_DIR/micapipe/sub-$ID/func" \
    -B "$SAVE_DIR/micapipe/sub-$ID/xfm/rest":"$SAVE_DIR/micapipe/sub-$ID/xfm" \
    $MICA_IMG \
    -threads $THREADS \
    -bids $BIDS_DIR \
    -out $SAVE_DIR \
    -tmpDir $TEMP_DIR \
    -sub $ID \
    -proc_rsfmri \
    -mainScanStr $REST_NAME

#### movie
# Create subdirectories
[[ ! -d "$SAVE_DIR/micapipe/sub-$ID/func/movie" ]] && mkdir "$SAVE_DIR/micapipe/sub-$ID/func/movie"
[[ ! -d "$SAVE_DIR/micapipe/sub-$ID/xfm/movie" ]] && mkdir "$SAVE_DIR/micapipe/sub-$ID/xfm/movie"

# Run
singularity run --cleanenv \
   -B $FREE_LIC:/opt/freesurfer-6.0.0/license.txt \
   -B "$BIDS_DIR/sub-$ID/func/sub-$ID""_movie.json":"$BIDS_DIR/task-rest_bold.json" \
   -B "$SAVE_DIR/micapipe/sub-$ID/func/movie":"$SAVE_DIR/micapipe/sub-$ID/func" \
   -B "$SAVE_DIR/micapipe/sub-$ID/xfm/movie":"$SAVE_DIR/micapipe/sub-$ID/xfm" \
   $MICA_IMG \
   -threads $THREADS \
   -bids $BIDS_DIR \
   -out $SAVE_DIR \
   -tmpDir $TEMP_DIR \
   -sub $ID \
   -proc_rsfmri \
   -mainScanStr $MOVIE_NAME

