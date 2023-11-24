#!/bin/bash
#SBATCH --job-name=eORCA1_ensemble
#SBATCH --time=24:00:00
#SBATCH --account=n01-nceo
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --ntasks=1

# SUBMIT WITH: sbatch --export=year=XXXX,n_ens=XXX cycle_year.slurm

export OMP_NUM_THREADS=1
source ../code/archer2-files/ucx_env
source set_environment.sh
expname=DAChlo
year=2015
n_ens=30
is_free=0

echo year $year
echo n_ens $n_ens

if [ ! -d $RUN_DIR ]; then
    $CYCLE_DIR/setup_run.sh $n_ens $is_free
fi

$CYCLE_DIR/setup_year.sh $year $n_ens

cd $RUN_DIR/namelists/
dt=`grep 'rn_rdt\s*=' namelist_cfg | tr -d '[:space:]' | cut -d'=' -f2 | cut -d'!' -f1` #get time-step
cd $RUN_DIR

set -e 
for month in 01; do
    echo "$SLURM_JOB_ID Submitting year/month" $year/$month >> $RUN_DIR/jobs.log
    # Set namelists
    iter_start=`cat $RUN_DIR/current_iter` #get iteration number
    printf -v iter_start_zero "%08d" $iter_start
    echo $month 
    nday=`cal $month $year | grep -v '[A-Za-z]' | wc -w`
    if [ $year == $ystart ] && [ $month == 01 ];
    then
       nday=$(($nday - 12))
    fi
    echo $nday

    iter_end=$(($iter_start + 86400*$nday/$dt))
    iter_start=$(($iter_start + 1))
    echo $iter_start $iter_end

    $CYCLE_DIR/update_nemo_nl --file $RUN_DIR/namelists/namelist_cfg  \
        --runid $NAME                 \
        --restart true            \
        --next_step $iter_start           \
        --final_step $iter_end          \
        --restart_file ${NAME}_${iter_start_zero}_restart \
        --ice_file $RUN_DIR/namelists/namelist_ice_cfg  \
        --ice_restart_file ${NAME}_${iter_start_zero}_restart_ice \
        --trc_file $RUN_DIR/namelists/namelist_top_cfg  \
        --trc_restart_file ${NAME}_${iter_start_zero}_restart_trc

    bash $CYCLE_DIR/change_stopack_seed.sh $RUN_DIR/namelists
    # for the initial free run, we need to generate
    # state vector perturbation from covariance matrix
    nday=`cal $month $year | grep -v '[A-Za-z]' | wc -w`
    if [ $year != $ystart ] && [ $month != 01 ] && [ $is_free != 1 ];
    then
        sed -i 's/type_ens_init = 3/type_ens_init = 4/g' $RUN_DIR/namelists/namelist_cfg.pdaf
        sed -i 's/ens_restart = .false./ens_restart = .true./g' $RUN_DIR/namelists/namelist_cfg.pdaf
    fi
    # use very large delt_obs to avoid DA in free run
    if [ $is_free == 1 ];
    then
        sed -i 's/delt_obs =/delt_obs = 9999/g' $RUN_DIR/namelists/namelist_cfg.pdaf
    fi
    # Launch Run
    echo "Launching $year $month at $(date +%s) seconds since 1970-01-01 00:00:00"
    sbatch --wait $RUN_DIR/submit.sh
    echo $iter_end > $RUN_DIR/current_iter
    if [ $is_free == 1 ];
    then
        # remove delt_obs in namelist
        sed -i 's/delt_obs = '${step}'/delt_obs =/g' $RUN_DIR/namelists/namelist_cfg.pdaf
    fi

    # Archive results
    for i in $(seq 1 $n_ens); do
        outdir=$OUTPUT_DIR/$expname/$year/$month/ensemble_$i
        mkdir -p $outdir
        mv ensemble_$i/$NAME*nc ensemble_$i/*.output ensemble_$i/run.stat* $outdir
    done
    outdir=$OUTPUT_DIR/$expname/$year/$month
    cp namelists/namelist* $outdir
    mv ensemble_1/state_*.nc $outdir
    cd $OUTPUT_DIR/$expname
    cp ../../prepostprocess/merge_domain.sh .
    # sbatch --export=year=$year,icycle=$month mergef.sh
    # sbatch --export=year=$year,icycle=$month merge_domain.sh
    cd $RUN_DIR
done
