#!/bin/bash --login
#SBATCH --job-name=PDAFstat
#SBATCH --time=12:00:00
#SBATCH --account=n01-nceo
#SBATCH --partition=serial
#SBATCH --qos=serial
#SBATCH --ntasks=1
source ../../code/archer2-files/ucx_env

year=2015
icycle=01

BaseDir=${year}/${icycle}/analysis/

set -e
for fname in $(ls $BaseDir/state_*_001.nc)
do
    if [[ $fname == *201504* ]]; then
        continue
    fi
    if [[ $fname == *201508* ]]; then
        continue
    fi
    if [[ $fname == *201512* ]]; then
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
