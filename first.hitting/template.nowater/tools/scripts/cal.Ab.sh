#!/bin/bash

make -C /home/mi/wanghan/study/noneq.msm/first.hitting/template.nowater/tools/analysis -j8

source parameters.sh

if test $# -lt 1; then
    echo usage:
    echo $0 gateValue
    exit
fi

rm -f base.info

for ii in $fht_gaussian_base_position;
do
    echo $fht_gaussian_base_max >> base.info
done

/home/mi/wanghan/study/noneq.msm/first.hitting/template.nowater/tools/analysis/cal.Ab.noAvg -g $1
