#!/bin/bash
###########
# Generates NIfTI files (with .json sidecars) from DICOM files.
# Intended for use with ALSPAC dataset only, probably doesn't generalise well.
# Requires array of IDs and file names given duplicates and discrepancies in file names.
# Assumed DICOM path: '<DATA_DIR>/<ID>/Connectom/scans/<SCAN>/DICOM/'
# Creates output file structure in BIDS format (or at least my best interpretation of it).
# Author: Chris Norman (NormanC4@cardiff.ac.uk)
##########
# SETTINGS
DATA_DIR=/scratch/scw1648/full_dataset
LINK_FIL=/scratch/scw1648/proj_cn/data/FileLinks_full_dataset.csv
SAVE_DIR=/scratch/scw1648/proj_cn/data/NIfTIs

######
# MAIN
OLDIFS=$IFS # Keep these safe for now
IFS=','
[ ! -f $LINK_FIL ] && { echo "Filename array: $LINK_FIL file not found"; exit 99; }
while read ID T1 RS RS_PA movie
do
    mkdir {"$SAVE_DIR/sub-$ID","$SAVE_DIR/sub-$ID/anat","$SAVE_DIR/sub-$ID/func"}
    if [ "$T1" != "" ]; then
        dcm2niix -b y -o "$SAVE_DIR/sub-$ID/anat" -f "sub-$ID""_T1" "$DATA_DIR/$ID/Connectom/scans/$T1/DICOM"
    fi
    if [ "$RS" != "" ]; then
        dcm2niix -b y -o "$SAVE_DIR/sub-$ID/func" -f "sub-$ID""_rest" "$DATA_DIR/$ID/Connectom/scans/$RS/DICOM"
    fi
    if [ "$RS_PA" != "" ]; then
        dcm2niix -b y -o "$SAVE_DIR/sub-$ID/func" -f "sub-$ID""_rest_PA" "$DATA_DIR/$ID/Connectom/scans/$RS_PA/DICOM"
    fi
    if [ "$movie" != "" ]; then
        dcm2niix -b y -o "$SAVE_DIR/sub-$ID/func" -f "sub-$ID""_movie" "$DATA_DIR/$ID/Connectom/scans/$movie/DICOM"
    fi
done < $LINK_FIL
IFS=$OLDIFS # Restore default IFS

