#!/bin/sh
#SBATCH --nodes=1
#SBATCH --time=06:00:00
#SBATCH --job-name=Archival
#SBATCH --exclude=ec1,ec2,ec3,ec4,ec5,ec6,ec7,ec8,ec9,ec10,ec11,ec12,ec13,ec14,ec15,ec16,ec17,ec18,ec19,ec20,ec21,ec22,ec23,ec24,ec25,ec26,ec27,ec28,ec29,ec30,ec31,ec32,ec33,ec34,ec35,ec36,ec37,ec38,ec39,ec40,ec41,ec42,ec43,ec44,ec45,ec46,ec47,ec48,ec49,ec50,ec51,ec52,ec53,ec54,ec55,ec56,ec57,ec58,ec59,ec60,ec61,ec62,ec63,ec64,ec65,ec66,ec67,ec68,ec69,ec70,ec71,ec72,ec73,ec74,ec75,ec76,ec77,ec113,ec114,ec115,ec116,ec117,ec118,ec119,ec120,ec121,ec122,ec123,ec124,ec125,ec126,ec127,ec128,ec129,ec130,ec131,ec132,ec133,ec134,ec135,ec136,ec137,ec138,ec139,ec140,ec141
#SBATCH --cores-per-socket=1 # Very high number to force unused node instead of shared 
#SBATCH --ntasks=2 # number of tasks to do, node says how to split amongst nodes 
#SBATCH --ntasks-per-node=2     # number of tasks per node
#SBATCH --ntasks-per-core=2
#SBATCH --cpu-freq=HighM1
#SBATCH --mem-per-cpu=12GB 
#SBATCH --output=Archival-%j.out

# ==========================================================
#module load anaconda/anaconda3
#source activate hoomd
#conda activate uber_env
# ==========================================================
python3 Parallel.py
rm Job_list.txt