#!/bin/bash
#SBATCH --job-name=eORCA1_ensemble
#SBATCH --time=24:00:00
#SBATCH --account=n01-nceo
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --ntasks=1

# submit job for a deterministic run (only one ensemble member)
# SUBMIT WITH: sbatch --export=year=XXXX,n_ens=XXX cycle_year.slurm

export OMP_NUM_THREADS=1
source ../code/archer2-files/ucx_env
source set_environment-spinup.sh
expname=determine
year=2015
n_ens=1
is_free=1

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
    # run one full year
    nday=365

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

    nday=`cal $month $year | grep -v '[A-Za-z]' | wc -w`
    # use very large delt_obs to avoid DA in free run
    sed -i 's/delt_obs = 1/delt_obs = 999999/g' $RUN_DIR/namelists/namelist_cfg.pdaf
    # setup the task the ensemble size
    sed -i 's/tasks =/tasks = '${n_ens}'/g' $RUN_DIR/namelists/namelist_cfg.pdaf
    # Launch Run
    echo "Launching $year $month at $(date +%s) seconds since 1970-01-01 00:00:00"
    sbatch --wait $RUN_DIR/submit.sh
    echo $iter_end > $RUN_DIR/current_iter
done
