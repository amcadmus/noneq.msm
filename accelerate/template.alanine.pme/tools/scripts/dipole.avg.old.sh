#!/bin/bash

make -C tools/average -j8

pert_main_dir=result.perts
backup_dir=../../../set.00/result.perts/`ls ../set.00/result.perts/ | grep backup | head -n 1`
backup_dir=/home/cocktail/wanghan/study/adress.noneq/alanine.nanma/results.save/ext.mode1.010.Ex.01.00.t1000ps.recheck/set.000/result.perts/backup.pert.000200

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
    echo 16 | g_dipoles -n $backup_dir/index.ndx  -f alanine.xtc -s $backup_dir//topol.tpr -nice 0 &> /dev/null
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
