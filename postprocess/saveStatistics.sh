#!/bin/bash --login
#SBATCH --job-name=stat
#SBATCH --time=12:00:00
#SBATCH --account=n01-nceo
#SBATCH --partition=serial
#SBATCH --qos=serial
#SBATCH --ntasks=32
source ../../code/archer2-files/ucx_env

year=2015
icycle=07

BaseDir=${year}/${icycle}

lfs setstripe -c -1 $BaseDir
cp ../../code/nemo/tools/REBUILD_NEMO/rebuild_nemo.exe $BaseDir/.
cp ../../code/nemo/tools/REBUILD_NEMO/nam_rebuild $BaseDir/.

set -e
for fname in $(ls $BaseDir/ensemble_1/eORCA1_1d_*_0000.nc | grep -v "globmean")
do
    if [[ $fname == *201504-201504* ]]; then
        continue
    fi
    if [[ $fname == *201508-201508* ]]; then
        continue
    fi
    if [[ $fname == *201512-201512* ]]; then
        continue
    fi
    start=`date +%s`
    Filename=$(basename $fname)
    filename=${Filename::-8} 
    echo $fname
    file_count=$(ls -1 $BaseDir/ensemble_1/${filename}_*.nc | wc -l)
    cp $BaseDir/ensemble_1/$filename*.nc $BaseDir/
    if [[ $filename == *ptrc1* ]]; then
        for ffname in $(ls $BaseDir/ensemble_1/$filename*.nc); do
            ffnamebase=$(basename $ffname)
            ffnamest=${ffnamebase::-8} 
            ffnameend=${ffnamebase: -8} 
            cp -v ${ffname} $BaseDir/$ffnamest-log$ffnameend
        done
    fi
    wait
    export OMP_NUM_THREADS=1
    if [[ $filename == *ptrc1* ]]; then
        srun -n 32 merge_multi $BaseDir $filename $file_count T
    fi
    srun -n 32 merge_multi $BaseDir $filename $file_count F
    wait
    echo "Current directory: $current_directory"
    rm $BaseDir/ensemble_*/$Filename

    # merge files
    cd $BaseDir
    sed -i 's/mesh_mask/'$filename'/g' nam_rebuild
    sed -i 's/ndomain=36/ndomain='$file_count'/g' nam_rebuild
    export OMP_NUM_THREADS=32
    srun -n 1 -c 32 rebuild_nemo.exe 
    wait
    sed -i 's/'$filename'/mesh_mask/g' nam_rebuild
    sed -i 's/ndomain='$file_count'/ndomain=36/g' nam_rebuild
    if [[ $filename == *ptrc1* ]]; then
        sed -i 's/mesh_mask/'$filename-log'/g' nam_rebuild
        sed -i 's/ndomain=36/ndomain='$file_count'/g' nam_rebuild
        export OMP_NUM_THREADS=32
        srun -n 1 -c 32 rebuild_nemo.exe 
        wait
        sed -i 's/'$filename-log'/mesh_mask/g' nam_rebuild
        sed -i 's/ndomain='$file_count'/ndomain=36/g' nam_rebuild
        rm ${filename}-log_0*.nc
    fi
    rm ${filename}_0*.nc
    cd ../..
    rm $BaseDir/ensemble_*/${filename}_*.nc
    current_directory=$(pwd)
    end=`date +%s`

    runtime=$((end-start))
    echo $runtime
done
