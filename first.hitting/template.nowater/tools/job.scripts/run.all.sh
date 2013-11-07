#!/bin/bash

for i in `ls | grep 'set.'`
do
    echo $i;
    cd $i
    nohup ./run.gaussian.sh &
    cd ..
done
