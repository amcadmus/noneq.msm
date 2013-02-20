#!/bin/bash

num_of_proc=5
seed_dir=set

stop_num=$(($num_of_proc-1))

for i in `seq 0 $stop_num`
do
    printi=`printf %02d $i`
    if test -d $seed_dir.$printi; then
	echo "existing dir $seed_dir.$printi, exit"
	exit
    fi
    echo "# prepare $seed_dir.$printi"
    cp -a $seed_dir $seed_dir.$printi
    cd $seed_dir.$printi
    sed -e "s/pert_parallel_num_pro=.*/pert_parallel_num_pro=$num_of_proc/g" parameters.sh |
    sed -e "s/pert_parallel_my_id=.*/pert_parallel_my_id=$i/g" > tmp.tmp.tmp
    mv -f tmp.tmp.tmp parameters.sh
    cd ..
done
    