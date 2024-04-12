
year=$1
n_ens=$2
INPUTS=$WORK/INPUTS

echo $INPUTS
echo $year

# Clean Previous year
yearm=$(($year - 1))
yearadd=$(($year + 1))

for i in $(seq 1 $n_ens)
do
    EnsRunDir=$RUN_DIR/ensemble_$i/
    # Link current year
    ln -s $INPUTS/PHYSICS/SBC/*y$year.nc $EnsRunDir/INPUTS/.
    wait
    ln -s $INPUTS/MEDUSA/SBC/Ndep/*y$year.nc $EnsRunDir/INPUTS/.
    wait
    echo $INPUTS/MEDUSA/SBC/pCO2/*y$year.nc 
    ln -s $INPUTS/MEDUSA/SBC/pCO2/*y$year.nc $EnsRunDir/INPUTS/.
    wait

    ln -s $INPUTS/PHYSICS/SBC/*y$yearadd.nc $EnsRunDir/INPUTS/.
    wait
    ln -s $INPUTS/MEDUSA/SBC/Ndep/*y$yearadd.nc $EnsRunDir/INPUTS/.
    wait
    ln -s $INPUTS/MEDUSA/SBC/pCO2/*y$yearadd.nc $EnsRunDir/INPUTS/.
    wait
done
