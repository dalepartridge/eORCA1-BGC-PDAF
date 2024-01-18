#!/bin/bash
#SBATCH --job-name=plot
#SBATCH --time=12:00:00
#SBATCH --account=n01-nceo
#SBATCH --partition=serial
#SBATCH --qos=serial
#SBATCH --ntasks=1

# SUBMIT WITH: sbatch --export=year=XXXX,n_ens=XXX cycle_year.slurm

export OMP_NUM_THREADS=1
source ../code/archer2-files/ucx_env
module load cray-python
export PYTHONUSERBASE=/work/n01/n01/ymchen/.local
export PATH=$PYTHONUSERBASE/bin:$PATH
export PYTHONPATH=$PYTHONUSERBASE/lib/python3.8/site-packages:$PYTHONPATH
export LD_LIBRARY_PATH="/work/n01/n01/ymchen/libgeos/lib64/:$LD_LIBRARY_PATH"

srun python -u test.py
