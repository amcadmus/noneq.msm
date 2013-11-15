#!/bin/bash

make -C /home/mi/wanghan/study/noneq.msm/first.hitting/template.inwater/tools/analysis -j8 &>/dev/null

source parameters.sh

if test $# -lt 1; then
    echo usage:
    echo $0 gateValue
    exit
fi

# rm -f base.info
# for ii in $fht_gaussian_base_position;
# do
#     echo $fht_gaussian_base_max >> base.info
# done

filename=fht.`printf %.3f $1`.out

if test -f $filename; then
    echo "# back up $filename"
    mv $filename $filename.`date +%s`
fi

/home/mi/wanghan/study/noneq.msm/first.hitting/template.inwater/tools/analysis/cal.Ab.noAvg -g $1 > $filename
