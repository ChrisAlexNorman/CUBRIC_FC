#!/bin/bash
###########
# Verifies existence of micapipe functional processing outputs against FileIndicators.csv.
# Author: Chris Norman (NormanC4@cardiff.ac.uk)
##########
# SETTINGS
DATA_DIR=/scratch/scw1648/proj_cn/data/mica_processed/micapipe
LINK_FIL=/scratch/scw1648/proj_cn/data/FileIndicators.csv
OUT_FIL=/scratch/scw1648/proj_cn/micapipe/pipeline_func/logs/func_dir_sizes.csv

######
# MAIN
printf "ID, rest, movie\n" > $OUT_FIL
OLDIFS=$IFS # Keep these safe for now
IFS=','
[ ! -f $LINK_FIL ] && { echo "Filename array: $LINK_FIL file not found"; exit 99; }
{ read
while read ID T1 RS RS_PA movie
do
    if [ "$RS" == "1" ]; then
        REST_SIZE=$(du -s "$DATA_DIR/sub-$ID/func/rest" | awk '{print $1}')
    else
        REST_SIZE="-" 
    fi
    if [ "$movie" == "1" ]; then
        MOVIE_SIZE=$(du -s "$DATA_DIR/sub-$ID/func/movie" | awk '{print $1}')
    else
        MOVIE_SIZE="-" 
    fi
    printf "$ID, $REST_SIZE, $MOVIE_SIZE\n" >> $OUT_FIL

done
} < $LINK_FIL
IFS=$OLDIFS # Restore default IFS

