
mkdir -p $RUN_DIR

cp $WORK/RUN/EXP00_MEDUSA/* $RUN_DIR
cp $WORK/code/nemo/cfgs/eORCA1-build/EXP00/nemo $RUN_DIR
cp $WORK/code/xios-build/bin/xios_server.exe $RUN_DIR
echo $starting_iter > $RUN_DIR/current_iter
printf -v iter_start_zero "%08d" $starting_iter

INPUTS=$WORK/INPUTS

# Domain and river files
ln -s $INPUTS/PHYSICS/DOM/eORCA_R1_zps_domcfg.nc $RUN_DIR
ln -s $INPUTS/PHYSICS/RIV/* $RUN_DIR

# Link restarts
mkdir -p $RUN_DIR/restarts
ln -s $INPUTS/PHYSICS/DOM/restart.nc $RUN_DIR/restarts/${NAME}_${iter_start_zero}_restart.nc
ln -s $INPUTS/PHYSICS/DOM/restart_ice.nc $RUN_DIR/restarts/${NAME}_${iter_start_zero}_restart_ice.nc


mkdir -p $RUN_DIR/INPUTS
ln -s $INPUTS/PHYSICS/DOM/weights* $RUN_DIR/INPUTS/
ln -s $INPUTS/PHYSICS/DOM/eddy_viscosity_3D.nc $RUN_DIR/INPUTS/
ln -s $INPUTS/PHYSICS/SBC/sss_absolute_salinity_WOA13_decav_Reg1L75_clim.nc $RUN_DIR/INPUTS/
ln -s $INPUTS/PHYSICS/BBC/* $RUN_DIR/INPUTS/
ln -s $INPUTS/PHYSICS/LBC/* $RUN_DIR/INPUTS/


# Link Medusa
ln -s $INPUTS/MEDUSA/DOM/restart_trc.nc $RUN_DIR/restarts/${NAME}_${iter_start_zero}_restart_trc.nc
ln -s $INPUTS/MEDUSA/SBC/*nc $RUN_DIR/INPUTS/
ln -s $INPUTS/MEDUSA/SBC/pCO2/pCO2a.nc $RUN_DIR/INPUTS/

ln -s $INPUTS/PHYSICS/DOM/eddy_viscosity_3D.nc $RUN_DIR
ln -s $INPUTS/PHYSICS/DOM/weights_coreII*nc $RUN_DIR
