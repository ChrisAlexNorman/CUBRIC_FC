#!/bin/bash
###########
# Saves smoothed and spike-regressed timeseries in fsaverage5 space
# Calculates FC and timeseries for specified parcellation in fsaverage5 space
# Author: Chris Norman (NormanC4@cardiff.ac.uk)
##########
# SETTINGS
DATA_DIR=/scratch/scw1648/proj_cn/data/mica_processed/micapipe
LINK_FIL=/scratch/scw1648/proj_cn/data/FileIndicators.csv
PARC_DIR=/scratch/scw1648/proj_cn/atlas_parcs/HCP-MMP1
PARC_FIL=micapipe_fsaverage5_glasser-360.annot
PARC_NAM=glasser-360
PY_FIL=/scratch/scw1648/proj_cn/micapipe/get_fsaverage5_FC.py

TODAY=$(date +"%Y-%m-%d")

source ~/.virtualenvs/ALSPAC-ENV/bin/activate
module load python/3.7.0
######
# MAIN
OLDIFS=$IFS # Keep these safe for now
IFS=','
[ ! -f $LINK_FIL ] && { echo "Filename array: $LINK_FIL file not found"; exit 99; }
{ read
while read ID T1 RS RS_PA movie
do
    if [ "$ID" == "1153" ]; then
	continue
    fi

    SUBJ="sub-$ID"

    if [ "$RS" == "1" ]; then
	SESS=rest
	FUNC_DIR="${DATA_DIR}/${SUBJ}/func/${SESS}"
	LOG_FIL="${DATA_DIR}/${SUBJ}/logs/CN_fsaverage5_${SESS}_${TODAY}.log"
	python3 $PY_FIL "$SUBJ" "$FUNC_DIR" "$PARC_DIR" "$PARC_FIL" "$PARC_NAM" > "$LOG_FIL"
    fi

    if [ "$movie" == "1" ]; then
        SESS=movie
	FUNC_DIR="${DATA_DIR}/${SUBJ}/func/${SESS}"
	LOG_FIL="${DATA_DIR}/${SUBJ}/logs/CN_fsaverage5_${SESS}_${TODAY}.log"
	python3 $PY_FIL "$SUBJ" "$FUNC_DIR" "$PARC_DIR" "$PARC_FIL" "$PARC_NAM" > "$LOG_FIL"
    fi

done
} < $LINK_FIL
IFS=$OLDIFS # Restore default IFS

