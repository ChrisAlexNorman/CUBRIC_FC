#!/bin/bash
#SBATCH --job-name=grad_loop_fsaverage5
#SBATCH --output=/scratch/scw1648/proj_cn/analysis/logs/grad_loop_fsaverage5.out
#SBATCH --error=/scratch/scw1648/proj_cn/analysis/logs/grad_loop_fsaverage5.err
#SBATCH -p compute
#SBATCH --time=0-01:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task 1
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=normanc4@cardiff.ac.uk
#SBATCH --mail-type=END

source ~/.virtualenvs/ALSPAC-ENV/bin/activate
module load python/3.7.0

start_time=$SECONDS
python3 ./gradients.py --mode "loop" --space "fsaverage5" --out "/scratch/scw1648/proj_cn/data/gradients/fsaverage5/"
end_time=$SECONDS
elapsed_time=$(( end_time - start_time ))
echo "Time elapsed: $elapsed_time seconds"
