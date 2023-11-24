n_ens=$1
is_freerun=$2
current_dir=$(pwd)

mkdir -p $RUN_DIR

mkdir -p $RUN_DIR/namelists

INPUTS=$WORK/INPUTS

cd $WORK/RUN/EXP00_MEDUSA/
find -maxdepth 1 -type f -exec cp {} $RUN_DIR/namelists \;

cd ${current_dir}
cp $WORK/code/nemo/cfgs/eORCA1-build/EXP00/nemo $RUN_DIR
cp $WORK/code/xios-build/bin/xios_server.exe $RUN_DIR

echo $starting_iter > $RUN_DIR/current_iter
printf -v iter_start_zero "%08d" $starting_iter


# Duplicate nemo context reference line in iodef file
for i in $(seq 1 $(($n_ens-1)));
do
    sed -i '0,/.*context_nemo.*/{s/.*context_nemo.*/&\n&/}' $RUN_DIR/namelists/iodef.xml
done

# Create numbered context nemo files and update context lines in iodef
for i in $(seq -f "%03g" 1 $n_ens);
do
    cp $RUN_DIR/namelists/context_nemo.xml $RUN_DIR/namelists/context_nemo_$i.xml
    ./setup-update_ensemble_files --con_file $RUN_DIR/namelists/context_nemo_$i.xml \
                        --ens_id $i \
                        --iodef_file $RUN_DIR/namelists/iodef.xml \
                        --ens_file $((10#$i))
done
rm $RUN_DIR/namelists/context_nemo.xml
ln -s $RUN_DIR/namelists/iodef.xml $RUN_DIR/iodef.xml

# ensemble based output file definitions
cd $RUN_DIR/namelists/
tmp_path=$(echo $RUN_DIR | sed 's/\//\\\//g')

for filedef in $(ls file_def*.xml );
do
    echo $filedef
    for j in $(seq 1 $n_ens);
    do
        cp $filedef ${filedef%.xml}_$j.xml
        sed -i 's/name="@expname@_@freq@_@startdate@_@enddate@"/name="'$tmp_path'\/ensemble_'$j'\/@expname@_@freq@_@startdate@_@enddate@"/g' ${filedef%.xml}_$j.xml
    done
done

ens_dirs='-D '
for i in $(seq 1 $n_ens)
do
    EnsRunDir=$RUN_DIR/ensemble_$i
    mkdir -p $EnsRunDir

    # Domain and river files
    ln -s $INPUTS/PHYSICS/DOM/eORCA_R1_zps_domcfg.nc $EnsRunDir/.
    ln -s $INPUTS/PHYSICS/RIV/* $EnsRunDir/.

    # Link restarts
    mkdir -p $EnsRunDir/restarts
    if [ $is_freerun -eq 1 ];
    then
        ln -s $INPUTS/pdaf/ensemble_$i/initPerturb/restart.nc $EnsRunDir/restarts/${NAME}_${iter_start_zero}_restart.nc
        ln -s $INPUTS/PHYSICS/DOM/restart_ice.nc $EnsRunDir/restarts/${NAME}_${iter_start_zero}_restart_ice.nc
        ln -s $INPUTS/MEDUSA/DOM/restart_trc.nc $EnsRunDir/restarts/${NAME}_${iter_start_zero}_restart_trc.nc
    fi
    ln -s $INPUTS/pdaf/ensemble_$i/* $EnsRunDir/restarts/
    

    mkdir -p $EnsRunDir/INPUTS
    ln -s $INPUTS/PHYSICS/DOM/weights* $EnsRunDir/INPUTS/.
    ln -s $INPUTS/PHYSICS/DOM/eddy_viscosity_3D.nc $EnsRunDir/INPUTS/.
    ln -s $INPUTS/PHYSICS/SBC/sss_absolute_salinity_WOA13_decav_Reg1L75_clim.nc $EnsRunDir/INPUTS/.
    ln -s $INPUTS/PHYSICS/BBC/* $EnsRunDir/INPUTS/.
    ln -s $INPUTS/PHYSICS/LBC/* $EnsRunDir/INPUTS/.

    # Link Medusa
    ln -s $INPUTS/MEDUSA/SBC/*nc $EnsRunDir/INPUTS/.
    ln -s $INPUTS/MEDUSA/SBC/pCO2/pCO2a.nc $EnsRunDir/INPUTS/.

    ln -s $INPUTS/PHYSICS/DOM/eddy_viscosity_3D.nc $EnsRunDir/.
    ln -s $INPUTS/PHYSICS/DOM/weights_coreII*nc $EnsRunDir/.

    # Link executable
    ln -s $RUN_DIR/nemo $EnsRunDir/.
    ln -s $RUN_DIR/xios_server.exe $EnsRunDir/.

    # Link namelists
    ln -s $RUN_DIR/namelists/* $EnsRunDir/.

    # perturbed parameter files
    unlink $EnsRunDir/fabm.yaml
    cp -v ${current_dir}/fabm_$i.yaml $EnsRunDir/fabm.yaml

    # Link initial covariance matrix
    cd ${current_dir}
    ln -s $INPUTS/pdaf/cov-log-2D.nc $EnsRunDir/cov.nc

    # Link observations
    ln -s $INPUTS/obs $EnsRunDir/.

    # prepare for slurm scripts
    ens_dirs=${ens_dirs}' ensemble_'$i
    lfs setstripe -c -1 $EnsRunDir
done
cp $EnsRunDir/file_def*.xml $RUN_DIR

# ./pert_params_medusa $RUN_DIR/namelists $n_ens $RUN_DIR
./mkslurm_hetjob_online_ensemble -a n01-nceo -S 16 -s 8 -C 180 -m 16 $ens_dirs -t 05:00:00 -p standard > submit.sh
mv submit.sh $RUN_DIR/
