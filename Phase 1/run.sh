#!/bin/sh
#SBATCH --nodes=4
#SBATCH --time=100:00:00
#SBATCH --job-name=hoomd
#SBATCH --cores-per-socket=14 # Very high number to force unused node instead of shared 
#SBATCH --ntasks-per-node=1     # number of tasks per node
#SBATCH --cpu-freq=HighM1
#SBATCH --mem-per-cpu=16GB 
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-user=<br696985@ucf.edu>
#SBATCH --output=hoomd-%j.out

# ==========================================================
#module load anaconda/anaconda3
#source activate hoomd
# ==========================================================
python3 run_code.py
