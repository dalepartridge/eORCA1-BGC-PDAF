
source ../code/0_set_environment.sh
INPUT_DIR=$WORK/INPUTS

# Source directories
PHYS_DIR=/work/n01/n01/jdha/scratch/eORCA1/nemo/cfgs/eORCA1/INPUTS
MEDUSA_DIR=/work/n01/n01/gle/eORCA1/INPUTS
OBS_DIR=/work/n01/n01/ymchen/obs
COV_DIR=/work/n01/n01/ymchen/cov

##### Link Physics Files #####
mkdir -p $INPUT_DIR/PHYSICS
cd $INPUT_DIR/PHYSICS

#Domain files
mkdir -p ./DOM
ln -s $PHYS_DIR/eORCA_R1_zps_domcfg.nc DOM/
ln -s $INPUT_DIR/../RAW_DATA/restart_20150101.nc DOM/restart.nc
ln -s $INPUT_DIR/../RAW_DATA/restart_ice_20150101.nc DOM/restart_ice.nc
ln -s $PHYS_DIR/weights* DOM/
ln -s $MEDUSA_DIR/../test1/nemo/cfgs/eORCA1/EXP00/eddy_viscosity_3D.nc DOM/

#Surface BCs
mkdir -p ./SBC
ln -s /work/n01/shared/nemo/FORCING/JRA/JRA_v1.5.0_rechunk/*.nc SBC/ # Atmospheric forcing
ln -s $PHYS_DIR/merged_ESACCI_BIOMER4V1R1_CHL_REG05.nc SBC/ #solar radiation
ln -s $PHYS_DIR/sss_absolute_salinity_WOA13_decav_Reg1L75_clim.nc SBC/ # surface salinity restore

#LBCs
mkdir -p ./LBC
ln -s $PHYS_DIR/shlat2d.nc LBC/

#Bottom BCs
mkdir -p ./BBC
ln -s $PHYS_DIR/Goutorbe_ghflux.nc BBC/
ln -s $PHYS_DIR/weights_Goutorbe1_2_eorca1_bilinear.nc BBC/

#River files
mkdir -p ./RIV
ln -s $PHYS_DIR/runoff*nc RIV/

##### Link ERSEM BGC files ######

#mkdir -p $INPUT_DIR/ERSEM
#cd $INPUT_DIR/ERSEM

#mkdir -p ./DOM
#ln -s sdgdhfdgh DOM/restart_trc.nc
#ln -s data_TRC_nomask????

#mkdir -p ./SBC
#ln -s $ERSEM_DIR/eORCA1-CCI-ady-broadband-climatology-1997-2020.nc SBC/
#ln -s $ERSEM_DIR/eORCA1_dust.nc SBC/

##### Link Medusa BGC files ######
mkdir -p $INPUT_DIR/MEDUSA
cd $INPUT_DIR/MEDUSA

mkdir -p ./DOM
ln -s $INPUT_DIR/../RAW_DATA/restart_trc_20000101.nc DOM/restart_trc.nc
#ln -s data_TRC_nomask????

mkdir -p ./SBC
ln -s $MEDUSA_DIR/eORCA1-CCI-ady-broadband-climatology-1997-2020.nc SBC/
ln -s $MEDUSA_DIR/eORCA1_dust.nc SBC/
ln -s $MEDUSA_DIR/eORCA1_Fe_dep_GESAMP.nc SBC/

mkdir -p ./SBC/pCO2
ln -s $MEDUSA_DIR/eORCA1_pCO2a*nc SBC/pCO2/
ln -s $MEDUSA_DIR/pCO2a.nc SBC/pCO2/

mkdir -p ./SBC/Ndep
ln -s $MEDUSA_DIR/ndep*nc SBC/Ndep/

mkdir -p ./RIV
ln -s $PHYS_DIR/runoff*nc RIV/

##### Link PDAF Files #####
cd $INPUT_DIR
mkdir -p ./pdaf
cp $NEMO_PDAF_CLONE/MY_SRC_pdaf/namelist_cfg.pdaf pdaf/
ln -s $COV_DIR/* pdaf/

##### Link Observation Files #####
mkdir -p ./obs
ln -s $OBS_DIR/* obs/