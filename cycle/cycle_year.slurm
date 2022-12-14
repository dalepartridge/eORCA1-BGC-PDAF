#!/bin/bash
#SBATCH --job-name=eORCA1
#SBATCH --time=03:20:00
#SBATCH --account=n01-NCEO
#SBATCH --partition=standard
#SBATCH --qos=standard

#SBATCH --nodes=4
#SBATCH --ntasks-per-core=1

# SUBMIT WITH: sbatch --export=year=XXXX cycle_year.slurm

export OMP_NUM_THREADS=1
source ../code/archer2-files/ucx_env
source set_environment.sh

if [ ! -d $RUN_DIR ]; then
    $CYCLE_DIR/setup_run.sh
fi

$CYCLE_DIR/setup_year.sh $year

cd $RUN_DIR
dt=`grep 'rn_rdt\s*=' namelist_cfg | tr -d '[:space:]' | cut -d'=' -f2 | cut -d'!' -f1` #get time-step

for month in 0{1..9} {10..12}; do

    echo "$SLURM_JOB_ID Submitting year/month" $year/$month >> $RUN_DIR/jobs.log
    
    # Set namelists
    iter_start=`cat $RUN_DIR/current_iter` #get iteration number
    printf -v iter_start_zero "%08d" $iter_start
    nday=`cal $month $year | grep -v '[A-Za-z]' | wc -w`
    iter_end=$(($iter_start + 86400*$nday/$dt))
    iter_start=$(($iter_start + 1))
    
    $CYCLE_DIR/update_nemo_nl --file $RUN_DIR/namelist_cfg  \
        --runid $NAME                 \
        --restart true            \
        --next_step $iter_start           \
        --final_step $iter_end          \
        --restart_file ${NAME}_${iter_start_zero}_restart \
        --ice_file $RUN_DIR/namelist_ice_cfg  \
        --ice_restart_file ${NAME}_${iter_start_zero}_restart_ice \
        --trc_file $RUN_DIR/namelist_top_cfg  \
        --trc_restart_file ${NAME}_${iter_start_zero}_restart_trc

    # Launch Run
    echo "Launching $year $month at $(date +%s) seconds since 1970-01-01 00:00:00"
    srun --mem-bind=local \
        --het-group=0 --nodes=3 --ntasks=171 --ntasks-per-node=57 --cpu-bind=v,mask_cpu:0x1,0x10000,0x40000,0x100000,0x400000,0x1000000,0x4000000,0x10000000,0x40000000,0x100000000,0x400000000,0x1000000000,0x4000000000,0x10000000000,0x40000000000,0x100000000000,0x400000000000,0x1000000000000,0x4000000000000,0x10000000000000,0x40000000000000,0x100000000000000,0x400000000000000,0x1000000000000000,0x4000000000000000,0x10000000000000000,0x40000000000000000,0x100000000000000000,0x400000000000000000,0x1000000000000000000,0x4000000000000000000,0x10000000000000000000,0x40000000000000000000,0x100000000000000000000,0x400000000000000000000,0x1000000000000000000000,0x4000000000000000000000,0x10000000000000000000000,0x40000000000000000000000,0x100000000000000000000000,0x400000000000000000000000,0x1000000000000000000000000,0x4000000000000000000000000,0x10000000000000000000000000,0x40000000000000000000000000,0x100000000000000000000000000,0x400000000000000000000000000,0x1000000000000000000000000000,0x4000000000000000000000000000,0x10000000000000000000000000000,0x40000000000000000000000000000,0x100000000000000000000000000000,0x400000000000000000000000000000,0x1000000000000000000000000000000,0x4000000000000000000000000000000,0x10000000000000000000000000000000,0x40000000000000000000000000000000 ./myscript_wrapper.sh \
        : --het-group=1 --nodes=1 --ntasks=9 --ntasks-per-node=9 --cpu-bind=v,mask_cpu:0x1,0x10000,0x40000,0x100000,0x400000,0x1000000,0x4000000,0x10000000,0x40000000 ./myscript_wrapper.sh

    # Archive results
    outdir=$OUTPUT_DIR/$year/$month
    mkdir -p $outdir
    mv $NAME*nc *.output $outdir
    cp namelist* $outdir
    
    echo $iter_end > $RUN_DIR/current_iter
done

yearp=$(($year + 1))
if [ $yearp -le $yend ]; then
   echo 'Submitting' $yearp '...'
   sbatch --export=year=$yearp cycle_year.slurm
   echo "Done."
else
   echo "All done."
fi

