#!/bin/sh
#SBATCH --nodes=1
#SBATCH --time=00:10:00
#SBATCH --job-name=setting_up
#SBATCH --cores-per-socket=1 # Very high number to force unused node instead of shared 
#SBATCH --ntasks=2 # number of tasks to do, node says how to split amongst nodes 
#SBATCH --ntasks-per-node=2     # number of tasks per node
#SBATCH --ntasks-per-core=2
#SBATCH --cpu-freq=HighM1
#SBATCH --mem-per-cpu=1GB 
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-user=<br696985@ucf.edu>
#SBATCH --output=job-%j.out

# ==========================================================
#module load anaconda/anaconda3
#source activate hoomd
#conda activate uber_env
# ==========================================================
python3 Parallel.py
DEPS=$(awk 'BEGIN{ORS=","} {print}' Job_list.txt | sed 's/,$//')
sbatch --dependency=afterany:$DEPS Kill.sh