#!/bin/bash

source parameters.sh

make -C tools/simulator makedir &> make.log
make -C tools/simulator -j8 &>> make.log

rm -f run.noneq.log

command="`pwd`/tools/simulator/$project_name.noneq \
--dt $dt --nst $nst -p 100000 --gamma $gamma --temperature $temperature \
$command_line_param_print \
--branch-feq $branch_feq \
--noneq-check-feq $noneq_check_feq --noneq-time $noneq_time \
--quench-temperature $quench_temperature --quench-time $quench_time \
--x-low $x0 --x-up $x1 --v-low $v0 --v-up $v1 --x-grid $nx --v-grid $nv \
--warm-time $warm_time --seed $seed --pert-strength $pert_strength"

echo "# run on `hostname`" &>> run.noneq.log
echo "# command is $command" &>> run.noneq.log
$command > output.noneq 2> error.noneq &
echo "# pid is $!" &>> run.noneq.log

