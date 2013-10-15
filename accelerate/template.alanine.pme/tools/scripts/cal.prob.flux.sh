#!/bin/bash

make -C /home/mi/wanghan/study/noneq.msm/accelerate/template.alanine.pme/tools/angles/ makedir
make -C /home/mi/wanghan/study/noneq.msm/accelerate/template.alanine.pme/tools/angles/ -j
make -C /home/mi/wanghan/study/noneq.msm/accelerate/template.alanine.pme/tools/average/ makedir
make -C /home/mi/wanghan/study/noneq.msm/accelerate/template.alanine.pme/tools/average/ -j

source parameters.sh
lag_time=`echo "1.0 / $pert_rescale" | bc -l`
echo calculate with actual lag time $lag_time

/home/mi/wanghan/study/noneq.msm/accelerate/template.alanine.pme/tools/angles/average.traj.corr -c $lag_time &> /dev/null
/home/mi/wanghan/study/noneq.msm/accelerate/template.alanine.pme/tools/average/integrate.file -f meta.flux.out -o meta.flux.inte.out
