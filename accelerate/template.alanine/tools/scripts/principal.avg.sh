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

rm -f principal.*name
count1=0
for i in $targets;
do
    echo doing $pert_main_dir/$i
    cd $pert_main_dir/$i
    if test ! -d $backup_dir; then
	echo "no backup dir $backup"
	exit
    fi
    rm -f moi.dat axis1.dat axis2.dat axis3.dat
    echo 16 | g_principal -n $backup_dir/index.ndx  -f alanine.xtc -s $backup_dir//topol.tpr -nice 0 &> /dev/null
    cd - &>/dev/null
    echo $pert_main_dir/$i/moi.dat >> principal.name
    echo $pert_main_dir/$i/axis1.dat >> principal.a1.name
    echo $pert_main_dir/$i/axis2.dat >> principal.a2.name
    echo $pert_main_dir/$i/axis3.dat >> principal.a3.name
    
    count1=$(($count1+1))
    if test $count1 -eq $count; then
	break;
    fi
done

echo doing average
tools/average/average.file -f principal.name -o avg.principal.out
echo doing average a1
tools/average/average.file -f principal.a1.name -o avg.principal.a1.out
echo doing average a2
tools/average/average.file -f principal.a2.name -o avg.principal.a2.out
echo doing average a3
tools/average/average.file -f principal.a3.name -o avg.principal.a3.out
