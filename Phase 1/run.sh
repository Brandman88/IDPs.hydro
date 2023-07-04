#!/bin/sh
#SBATCH --nodes=3
#SBATCH --time=72:00:00
#SBATCH --job-name=hoomd
#SBATCH --ntasks-per-node=8     # number of tasks per node
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-user=<brandon.stine@knights.ucf.edu>
#SBATCH --output=hoomd-%j.out

# ==========================================================
#module load anaconda/anaconda3
#source activate hoomd
# ==========================================================
python3 run_code.py
