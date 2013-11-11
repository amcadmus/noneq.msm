#!/bin/bash

num_of_proc=16
every=1
seed_dir=set

stop_num=$(($num_of_proc-1))

for i in `seq 0 $stop_num`
do
    printi=`printf %03d $i`
    if test -d $seed_dir.$printi; then
	echo "existing dir $seed_dir.$printi, exit"
	exit
    fi
    tmpi=`echo "$i % $every" | bc`
    test $tmpi -ne 0 && continue
    echo "# prepare $seed_dir.$printi"
    cp -a $seed_dir $seed_dir.$printi
    cd $seed_dir.$printi
    sed -e "s/fht_parallel_num_pro=.*/fht_parallel_num_pro=$num_of_proc/g" parameters.sh |
    sed -e "s/fht_parallel_my_id=.*/fht_parallel_my_id=$i/g" > tmp.tmp.tmp
    mv -f tmp.tmp.tmp parameters.sh
    cd ..
done
    