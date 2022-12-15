
year=$1
INPUTS=$WORK/INPUTS


# Clean Previous year
yearm=$(($year - 1))

for i in $(seq 1 $n_ens)
do
    EnsRunDir=$RUN_DIR/ensemble_$i/
    rm -f $EnsRunDir/INPUTS/*y$yearm.nc

    # Link current year
    ln -s $INPUTS/PHYSICS/SBC/*y$year.nc $EnsRunDir/INPUTS/
    ln -s $INPUTS/MEDUSA/SBC/Ndep/*y$year.nc $EnsRunDir/INPUTS/
    ln -s $INPUTS/MEDUSA/SBC/pCO2/*y$year.nc $EnsRunDir/INPUTS/
done