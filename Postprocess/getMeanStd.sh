#!/bin/bash --login
#SBATCH --job-name=mergeoutput
#SBATCH --time=12:00:00
#SBATCH --account=n01-nceo
#SBATCH --partition=serial
#SBATCH --qos=serial
#SBATCH --ntasks=32
source ../../code/archer2-files/ucx_env

BaseDir=${year}/${icycle}

lfs setstripe -c -1 $BaseDir
cp ../../code/nemo/tools/REBUILD_NEMO/rebuild_nemo.exe $BaseDir/.
cp ../../code/nemo/tools/REBUILD_NEMO/nam_rebuild $BaseDir/.


set -e
for fname in $(ls $BaseDir/ensemble_1/eORCA1_1d_*0000.nc | grep -v "globmean")
do
    start=`date +%s`

    filename=$(basename $fname)
    filename=${filename::-8}
    file_count=$(ls -1 $BaseDir/ensemble_1/${filename}_*.nc | wc -l)
    cp $BaseDir/ensemble_1/$Filename $BaseDir/
    wait
    sleep 3
    export OMP_NUM_THREADS=1
    srun -n 32 meanStdForSubdomainFiles $BaseDir $filename $file_count
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
    sleep 3
    sed -i 's/'$filename'/mesh_mask/g' nam_rebuild
    sed -i 's/ndomain='$file_count'/ndomain=36/g' nam_rebuild
    rm ${filename}_0*.nc
    cd ../..
    current_directory=$(pwd)
    end=`date +%s`

    runtime=$((end-start))
    echo $runtime
done
