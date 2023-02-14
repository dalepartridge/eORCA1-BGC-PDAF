#!/bin/bash
#SBATCH --job-name=eORCA1_ensemble
#SBATCH --time=24:00:00
#SBATCH --account=n01-nceo
#SBATCH --partition=serial
#SBATCH --qos=serial
#SBATCH --nodes=1
#SBATCH --ntasks-per-core=1

# SUBMIT WITH: sbatch --export=year=XXXX --export=n_ens=XXX cycle_year.slurm

export OMP_NUM_THREADS=1
source ../code/archer2-files/ucx_env
source set_environment.sh

if [ ! -d $RUN_DIR ]; then
    $CYCLE_DIR/setup_run.sh $n_ens
fi

$CYCLE_DIR/setup_year.sh $year $n_ens

cd $RUN_DIR/namelists/
dt=`grep 'rn_rdt\s*=' namelist_cfg | tr -d '[:space:]' | cut -d'=' -f2 | cut -d'!' -f1` #get time-step
cd $RUN_DIR

month=04
for month in 01 04 07 10; do
    echo "$SLURM_JOB_ID Submitting year/month" $year/$month >> $RUN_DIR/jobs.log
    
    # Set namelists
    iter_start=`cat $RUN_DIR/current_iter` #get iteration number
    printf -v iter_start_zero "%08d" $iter_start
    nday=`cal $month $year | grep -v '[A-Za-z]' | wc -w`
    nday1=`cal $(($month+1)) $year | grep -v '[A-Za-z]' | wc -w`
    nday2=`cal $(($month+2)) $year | grep -v '[A-Za-z]' | wc -w`
    nday=$(($nday+$nday1+$nday2))
    echo $nday
    iter_end=$(($iter_start + 86400*$nday/$dt))
    iter_start=$(($iter_start + 1))
    
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

    # update the first assimilation time
    nday=`cal $month $year | grep -v '[A-Za-z]' | wc -w`
    step=$((86400*$nday/2/$dt))
    if [ $year != $ystart ] || [ $month != 01 ];
    then
        sed -i 's/type_ens_init = 3/type_ens_init = 4/g' $RUN_DIR/namelist_cfg.pdaf
        sed -i 's/ens_restart = .false./ens_restart = .true./g' $RUN_DIR/namelist_cfg.pdaf
    fi
    sed -i 's/delt_obs =/delt_obs = '${step}'/g' $RUN_DIR/namelist_cfg.pdaf

    # Launch Run
    echo "Launching $year $month at $(date +%s) seconds since 1970-01-01 00:00:00"
    sbatch --wait $RUN_DIR/submit.sh

    # remove delt_obs in namelist
    sed -i 's/delt_obs = '${step}'/delt_obs =/g' $RUN_DIR/namelist_cfg.pdaf
    s=`grep delt_obs $RUN_DIR/namelist_cfg.pdaf`
    echo $s

    # Archive results
    for i in {1..30}; do
        outdir=$OUTPUT_DIR/$year/$month/ensemble_$i
        mkdir -p $outdir
        mv ensemble_$i/$NAME*nc ensemble_$i/*.output $outdir
    done
    outdir=$OUTPUT_DIR/$year/$month
    cp namelists/namelist* $outdir
    mv ensemble_1/state_*.nc $outdir
    mv ensemble_1/variance_*.nc $outdir
    
    echo $iter_end > $RUN_DIR/current_iter
done

yearp=$(($year + 1))
if [ $yearp -le $yend ]; then
   echo 'Submitting' $yearp '...'
   sbatch --export=year=$yearp --export=n_ens=$n_ens cycle_year.slurm
   echo "Done."
else
   echo "All done."
fi

