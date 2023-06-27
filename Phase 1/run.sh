#!/bin/sh
#SBATCH --nodes=8
#SBATCH --time=8:00:00
#SBATCH --job-name=hoomd
#SBATCH --ntasks-per-node=1     # number of tasks per node
#SBATCH --cpus-per-task=8        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-user=<brandon.stine@knights.ucf.edu>
#SBATCH --output=hoomd-%j.out

# ==========================================================
#module load anaconda/anaconda3
#source activate hoomd
# ==========================================================
python3 run_code.py
