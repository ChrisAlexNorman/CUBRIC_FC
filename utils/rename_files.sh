#!/bin/bash
###########
# Simple script to rename a subset of data files
# In the initial case files *T1.* -> *T1w.* the expected format for mica-pipe.
# Author: Chris Norman (NormanC4@cardiff.ac.uk)
##########
# SETTINGS
DATA_DIR=/scratch/scw1648/proj_cn/data/NIfTIs
INDI_FIL=/scratch/scw1648/proj_cn/data/FileIndicators.csv

######
# MAIN
OLDIFS=$IFS # Keep these safe for now
IFS=','
[ ! -f $INDI_FIL ] && { echo "Filename array: $LINK_FIL file not found"; exit 99; }
while read ID T1 RS RS_PA movie
do
    if [[ "$T1" == 1 ]]; then
        mv "$DATA_DIR/sub-$ID/anat/sub-$ID""_T1.nii" "$DATA_DIR/sub-$ID/anat/sub-$ID""_T1w.nii"
        mv "$DATA_DIR/sub-$ID/anat/sub-$ID""_T1.json" "$DATA_DIR/sub-$ID/anat/sub-$ID""_T1w.json"
    fi
done < $INDI_FIL
IFS=$OLDIFS # Restore default IFS

