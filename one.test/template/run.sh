#!/bin/bash

source parameters.sh

make -C tools/simulator makedir &> make.log
make -C tools/simulator -j8 &>> make.log

command="`pwd`/tools/simulator/timeCorr \
--dt $dt --nst $nst --nstprint 100000 --gamma $gamma --temperature $temperature \
--double-well-k $double_well_k --double-well-a $double_well_a \
--branch-feq $branch_feq \
--noneq-check-feq $noneq_check_feq --noneq-time $noneq_time \
--quench-temperature $quench_temperature --quench-time $quench_time \
--corr-time $corr_time --corr-step $corr_step \
--warm-time $warm_time --seed $seed --pert-strength $pert_strength"

echo "# command is $command"
$command > output.timeCorr &
echo "# pid is $!"

command="`pwd`/tools/simulator/noneq \
--dt $dt --nst $nst --nstprint 100000 --gamma $gamma --temperature $temperature \
--double-well-k $double_well_k --double-well-a $double_well_a \
--branch-feq $branch_feq \
--noneq-check-feq $noneq_check_feq --noneq-time $noneq_time \
--quench-temperature $quench_temperature --quench-time $quench_time \
--warm-time $warm_time --seed $seed --pert-strength $pert_strength"

echo "# command is $command"
$command > output.noneq &
echo "# pid is $!"

