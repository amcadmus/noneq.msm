#!/bin/bash

num_of_proc=5
seed_dir=set
stop_num=$(($num_of_proc-1))

cwd=`pwd`

rm -f $seed_dir/angle.name
rm -f $seed_dir/gxs.name
rm -f $seed_dir/energy.name

for i in `seq 0 $stop_num`
do
    printi=`printf %02d $i`

    echo "# collect $seed_dir.$printi"
    cd $seed_dir
    test ! -d result.perts && mkdir result.perts
    cd result.perts
    if test -d ../../$seed_dir.$printi/result.perts; then
	ln -sf ../../$seed_dir.$printi/result.perts/pert.*   .
    fi
    cd ..
    test -f ../$seed_dir.$printi/angle.name  && cat ../$seed_dir.$printi/angle.name  >> ./angle.name
    test -f ../$seed_dir.$printi/gxs.name  && cat ../$seed_dir.$printi/gxs.name  >> ./gxs.name
    test -f ../$seed_dir.$printi/energy.name && cat ../$seed_dir.$printi/energy.name >> ./energy.name
    
    cd ..
done
    