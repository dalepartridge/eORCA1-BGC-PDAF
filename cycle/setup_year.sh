
year=$1
INPUTS=$WORK/INPUTS

# Clean Previous year
yearm=$(($year - 1))
rm -f $RUN_DIR/INPUTS/*y$yearm.nc

# Link current year
ln -s $INPUTS/PHYSICS/SBC/*y$year.nc $RUN_DIR/INPUTS/
ln -s $INPUTS/MEDUSA/SBC/Ndep/*y$year.nc $RUN_DIR/INPUTS/
ln -s $INPUTS/MEDUSA/SBC/pCO2/*y$year.nc $RUN_DIR/INPUTS/
