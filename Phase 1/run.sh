#!/bin/sh
#SBATCH --nodes=3
#SBATCH --time=80:00:00
#SBATCH --job-name=hoomd
#SBATCH --cores-per-socket=14 # Very high number to force unused node instead of shared 
#SBATCH --ntasks=12 # number of tasks to do, node says how to split amongst nodes 
#SBATCH --ntasks-per-node=4     # number of tasks per node
#SBATCH --ntasks-per-core=2
#SBATCH --cpu-freq=HighM1
#SBATCH --mem-per-cpu=12GB 
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-user=<br696985@ucf.edu>
#SBATCH --output=hoomd-%j.out

# ==========================================================
#module load anaconda/anaconda3
#source activate hoomd
# ==========================================================
python3 run_code.py
