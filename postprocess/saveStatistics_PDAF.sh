#!/bin/bash --login
#SBATCH --job-name=PDAFstat
#SBATCH --time=12:00:00
#SBATCH --account=n01-nceo
#SBATCH --partition=serial
#SBATCH --qos=serial
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
source ../../code/archer2-files/ucx_env

BaseDir=${year}/${icycle}/analysis/

set -e
for fname in $(ls $BaseDir/state_*_001.nc)
do
    if [[ $fname == *201604* ]]; then
        continue
    fi
    if [[ $fname == *201608* ]]; then
        continue
    fi
    if [[ $fname == *201612* ]]; then
        continue
    fi
    start=`date +%s`
    Filename=$(basename $fname)
    filename=${Filename::-7} 
    echo $filename
    cp -v $fname $BaseDir/${filename}.nc
    wait
    export OMP_NUM_THREADS=1
    srun merge_pdaf $BaseDir $filename
    wait
    rm $BaseDir/${filename}_0*.nc

    end=`date +%s`
    runtime=$((end-start))
    echo $runtime
done
