#!/bin/bash --login
#SBATCH --job-name=rebuild
#SBATCH --time=12:00:00
#SBATCH --account=n01-nceo
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128

export OMP_NUM_THREADS=1
source /work/n01/n01/ymchen/eORCA1-BGC-PDAF/code/archer2-files/ucx_env

BaseDir=${year}/${icycle}

lfs setstripe -c -1 $BaseDir

cd $BaseDir
for i in  $(seq $nstart $nend)
do
    for fname in $(ls ensemble_$i/eORCA1_1m_*_globmean_*_0000.nc)
    do
        echo $fname
        mv $fname ${fname::-8}.nc
        rm ${fname::-8}_0*.nc
    done    
done

set -e
for i in  $(seq $nstart $nend)
do
    cp /work/n01/n01/ymchen/eORCA1-BGC-PDAF/code/nemo/tools/REBUILD_NEMO/rebuild_nemo.exe ensemble_$i/.
    cp /work/n01/n01/ymchen/eORCA1-BGC-PDAF/code/nemo/tools/REBUILD_NEMO/nam_rebuild ensemble_$i/.    
    for fname in $(ls ensemble_$i/eORCA1_*_0000.nc | grep -v "globmean")
    do
        if [[ $fname == *eORCA1_5d_* ]] || [[ $fname == *eORCA1_1m_* ]] || [[ $fname == *201504-201504* ]] || [[ $fname == *201508-201508* ]] || [[ $fname == *201512-201512* ]]; then
            filename=$(basename "$fname")
            echo ensemble_$i, ${filename::-8}
            cd ensemble_$i
            file_count=$(ls -1 ${filename::-8}_*.nc | wc -l)
            sed -i 's/mesh_mask/'${filename::-8}'/g' nam_rebuild
            sed -i 's/ndomain=36/ndomain='$file_count'/g' nam_rebuild
            sleep 10
            export OMP_NUM_THREADS=128
            srun -n 1 -c 128 rebuild_nemo.exe 
            export OMP_NUM_THREADS=1
            wait
            sed -i 's/'${filename::-8}'/mesh_mask/g' nam_rebuild
            sed -i 's/ndomain='$file_count'/ndomain=36/g' nam_rebuild
            sleep 10
            rm ${filename::-8}_*.nc
            cd ../
        fi
    done
done
