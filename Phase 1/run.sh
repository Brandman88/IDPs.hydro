#!/bin/sh
#SBATCH --nodes=1
#SBATCH --time=100:00:00
#SBATCH --job-name=hoomd
#SBATCH --cores-per-socket=8
#SBATCH --ntasks-per-node=4     # number of tasks per node
#SBATCH --cpu-freq=HighM1[Performance]
#SBATCH --shared
#SBATCH --mem-per-cpu=16GB 
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-user=<brandon.stine@knights.ucf.edu>
#SBATCH --output=hoomd-%j.out

# ==========================================================
#module load anaconda/anaconda3
#source activate hoomd
# ==========================================================
python3 run_code.py
