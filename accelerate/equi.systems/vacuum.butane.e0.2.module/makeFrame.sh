#!/bin/bash

rm -f equi.frame

dt=`grep dt grompp.mdp | cut -d "=" -f 2`
nsteps=`grep nsteps grompp.mdp | cut -d "=" -f 2`
nstout=`grep nstvout grompp.mdp | cut -d "=" -f 2`
begin=100
end=`echo "$dt * $nsteps" | bc -l`
step=`echo "$dt * $nstout" | bc -l`

count=1
for i in `seq $begin $step $end`
do
    echo $count
    pi=`printf %06d $count`
    time=`echo "$i - $step/2.0" | bc -l`
    echo $pi $time >> equi.frame
    count=$(($count+1))
done
