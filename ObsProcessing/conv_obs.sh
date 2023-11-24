#!/bin/bash --login
#SBATCH --job-name=obsconv
#SBATCH --time=12:00:00
#SBATCH --account=n01-nceo
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --ntasks-per-node=64
#SBATCH --ntasks=128
#SBATCH --cpus-per-task=1

export OMP_NUM_THREADS=1
module load cray-python
export PYTHONUSERBASE=/work/n01/n01/ymchen/.local
export PATH=$PYTHONUSERBASE/bin:$PATH
export PYTHONPATH=$PYTHONUSERBASE/lib/python3.8/site-packages:$PYTHONPATH
export LD_LIBRARY_PATH="/work/n01/n01/ymchen/libgeos/lib64/:$LD_LIBRARY_PATH"
srun --ntasks-per-node=64 --ntasks=128 --cpus-per-task=1 python -u getObs.py
