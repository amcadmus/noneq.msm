#!/bin/bash

make -C tools/average -j8

pert_main_dir=result.perts
backup_dir=../../../set.000/result.perts/`ls ../set.000/result.perts/ | grep backup | head -n 1`

targets=`ls $pert_main_dir | grep ^pert`

count=0
for i in $targets;
do
    count=$(($count+1))
done
count=$(($count-1))

rm -f dipole.name
count1=0
for i in $targets;
do
    echo doing $pert_main_dir/$i
    cd $pert_main_dir/$i
    if test ! -d $backup_dir; then
	echo "no backup dir $backup"
	exit
    fi
    rm -f Mtot.xvg
    echo 2 | g_dipoles -f butane.xtc -s $backup_dir//topol.tpr -nice 0 &> /dev/null
    rm -f epsilon.xvg aver.xvg dipdist.xvg
    cd - &>/dev/null
    echo $pert_main_dir/$i/Mtot.xvg >> dipole.name
    
    count1=$(($count1+1))
    if test $count1 -eq $count; then
	break;
    fi
done

echo doing average
tools/average/average.file -f dipole.name -o avg.dipole.out
