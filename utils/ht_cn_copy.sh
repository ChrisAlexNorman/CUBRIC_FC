#!/bin/bash
###########
# Script to copy files between HT and CN project directories
# Author: Chris Norman (NormanC4@cardiff.ac.uk)
##########
# SETTINGS
CN_DIR=/scratch/scw1648/proj_cn/data/NIfTI_Data
HT_DIR=/scratch/scw1648/proj_ht/prac/PR002_topup
INDI_FIL=/scratch/scw1648/proj_cn/data/FileIndicators.csv

######
# MAIN
OLDIFS=$IFS # Keep these safe for now
IFS=','
[ ! -f $INDI_FIL ] && { echo "Filename array: $LINK_FIL file not found"; exit 99; }
while read ID T1 RS RS_PA movie
do
    if [[ "$RS" == 1 ]]; then
        cp "$HT_DIR/sub-$ID/func/sub-$ID""_rest_unwarped.nii" "$CN_DIR/sub-$ID/func/sub-$ID""_rest_unwarped.nii"
    fi

    if [[ "$movie" == 1 ]]; then
        cp "$HT_DIR/sub-$ID/func/sub-$ID""_movie_unwarped.nii" "$CN_DIR/sub-$ID/func/sub-$ID""_movie_unwarped.nii"
    fi
done < $INDI_FIL
IFS=$OLDIFS # Restore default IFS

