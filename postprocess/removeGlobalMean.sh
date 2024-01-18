#!/bin/bash --login

year=2015
icycle=01
BaseDir=${year}/${icycle}
cd $BaseDir
for i in  {1..30}
do
    for fname in $(ls ensemble_$i/eORCA1_1m_*_globmean_*_0000.nc)
    do
        echo $fname
        mv $fname ${fname::-8}.nc
        rm ${fname::-8}_0*.nc
    done    
done