#!/bin/bash

make -C tools/average -j8

pert_main_dir=result.perts

targets=`ls $pert_main_dir | grep ^pert`

count=0
for i in $targets;
do
    count=$(($count+1))
done
count=$(($count-1))

rm -f energy.name
count1=0
for i in $targets;
do
    echo doing $pert_main_dir/$i
    cd $pert_main_dir/$i
    echo 14 15 17 19 | g_energy -xvg none -nice 0 &> /dev/null
    cd - &>/dev/null
    echo $pert_main_dir/$i/energy.xvg >> energy.name
    
    count1=$(($count1+1))
    if test $count1 -eq $count; then
	break;
    fi
done

echo doing average
tools/average/average.file -f energy.name -o avg.energy.out
