#!/bin/bash

if test $# -ne 3; then
    echo "num of command line parameter is $#"
    echo "usage:"
    echo "$0" start step end
    exit
fi

# echo "$0 $1 $2 $3"

for i in `seq $1 $2 $3`
do
    printi=`printf %03d $i`
    dir=set.$printi
    echo $dir;
    cd $dir
    nohup ./run.sh &
    cd ..
done