#!/bin/bash
#SBATCH --job-name=micapipes_func
#SBATCH -p htc
#SBATCH --time=1-00:00
#SBATCH --cpus-per-task 20
#SBATCH --mem-per-cpu=5000
#SBATCH --array=0-36
#SBATCH -o /scratch/scw1648/proj_cn/micapipe/pipeline_func/logs/task-%a.out
#SBATCH -e /scratch/scw1648/proj_cn/micapipe/pipeline_func/logs/task-%a.err

# Submit 37 jobs (max. 40 sumitted jobs allowed by SCW Hawk (HTC node)).
# Each job runs 6 IDs in series. Last job runs remaining 4 IDs

declare -i ID_START=1001+6*$SLURM_ARRAY_TASK_ID

srun /scratch/scw1648/proj_cn/micapipe/pipeline_func/micapipe_func_sub.sh $ID_START

