#!/bin/bash

s1=843052
s2=378174
s3=60998
s4=398305

for i in {1..30};
do
    unlink $1/../ensemble_$i/namelist_cfg
	sed 's/nn_stopack_seed=843052,378174,60998,398305,/nn_stopack_seed='$s1','$s2','$s3','$s4',/g' $1/namelist_cfg > $1/../ensemble_$i/namelist_cfg
	s1=$((s1+1))
	s2=$((s2+1))
	s3=$((s3+1))
	s4=$((s4+1))
done
