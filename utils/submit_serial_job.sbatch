#!/bin/bash
#SBATCH --job-name=fsavg5_conv
#SBATCH --output=/scratch/scw1648/proj_cn/micapipe/fsavg5_conv.out
#SBATCH --error=/scratch/scw1648/proj_cn/micapipe/fsavg5_conv.err
#SBATCH -p compute
#SBATCH --time=0-03:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task 1
#SBATCH --mem-per-cpu=5000
#SBATCH --ntasks-per-node=1
#SBATCH --mail-type=END
#SBATCH --mail-user=normanc4@cardiff.ac.uk

/scratch/scw1648/proj_cn/micapipe/get_fsaverage5_FCs.sh
