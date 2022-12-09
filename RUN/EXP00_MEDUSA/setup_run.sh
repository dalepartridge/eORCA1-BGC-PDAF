
year=2000

WORK=/work/n01/n01/dapa/NCEO/eORCA1-BGC-PDAF

#cp $WORK/RUN/EXP00_MEDUSA/* .
cp $WORK/code/nemo/cfgs/eORCA1-build/EXP00/nemo .
cp $WORK/code/xios-build/bin/xios_server.exe .

INPUTS=$WORK/INPUTS

# Domain and river files
ln -s $INPUTS/PHYSICS/DOM/eORCA_R1_zps_domcfg.nc .
ln -s $INPUTS/PHYSICS/RIV/* .

# Link restarts
mkdir restarts
ln -s $INPUTS/PHYSICS/DOM/restart*nc restarts/

mkdir INPUTS
ln -s $INPUTS/PHYSICS/DOM/weights* INPUTS/
ln -s $INPUTS/PHYSICS/DOM/eddy_viscosity_3D.nc INPUTS/
ln -s $INPUTS/PHYSICS/SBC/*y$year.nc INPUTS/
ln -s $INPUTS/PHYSICS/SBC/weights*.nc INPUTS/
ln -s $INPUTS/PHYSICS/SBC/sss_absolute_salinity_WOA13_decav_Reg1L75_clim.nc INPUTS/
ln -s $INPUTS/PHYSICS/BBC/* INPUTS/
ln -s $INPUTS/PHYSICS/LBC/* INPUTS/


# Link Medusa
ln -s $INPUTS/MEDUSA/DOM/* restarts/
ln -s $INPUTS/MEDUSA/SBC/*nc INPUTS/
ln -s $INPUTS/MEDUSA/SBC/Ndep/*y$year.nc INPUTS/
ln -s $INPUTS/MEDUSA/SBC/pCO2/*y$year.nc INPUTS/
ln -s $INPUTS/MEDUSA/SBC/pCO2/pCO2a.nc INPUTS/

ln -s INPUTS/eddy_viscosity_3D.nc .
ln -s INPUTS/weights_coreII*nc .

