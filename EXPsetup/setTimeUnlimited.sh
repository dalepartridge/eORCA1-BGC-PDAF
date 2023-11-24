#!/bin/bash --login
#SBATCH --job-name=setTimeUnlimited
#SBATCH --time=12:00:00
#SBATCH --account=n01-nceo
#SBATCH --partition=serial
#SBATCH --qos=serial
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32

module load nco
export OMP_NUM_THREADS=32

original_path=/work/n01/n01/dapa/NCEO/eORCA1-ERA5
original_path=/work/n01/n01/ymchen/eORCA1-BGC-PDAF/prepostprocess/INPUT
for i in {0..9}
do
    # mkdir INPUT/ERA5_ens_$i
    for f in $original_path/ERA5_ens_$i/ERA5_msr_y2015.nc $original_path/ERA5_ens_$i/ERA5_msdwlwrf_y2015.nc  $original_path/ERA5_ens_$i/ERA5_msdwswrf_y2015.nc $original_path/ERA5_ens_$i/ERA5_mtpr_y2015.nc 
    do
        # echo $f
        # basefilename="$(basename -- $f)"
        # echo $basefilename
        ncks --mk_rec_dmn time -O $f $f # -o /work/n01/n01/ymchen/eORCA1-BGC-PDAF/prepostprocess/INPUT/ERA5_ens_$i/$basefilename
    done
done