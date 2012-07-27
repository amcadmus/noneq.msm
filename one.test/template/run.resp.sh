#!/bin/bash

source parameters.sh

make -C tools/simulator makedir &> make.log
make -C tools/simulator -j8 &>> make.log

rm -f run.resp.log

command="`pwd`/tools/simulator/$project_name.noneqResp \
--dt $dt --nst $nst -p 100000 --gamma $gamma --temperature $temperature \
$command_line_param_print \
--branch-feq $branch_feq \
--noneq-check-feq $noneq_check_feq --noneq-time $noneq_time \
--quench-temperature $quench_temperature --quench-time $quench_time \
--x-low $x0 --x-up $x1 --v-low $v0 --v-up $v1 --x-grid $nx --v-grid $nv \
--warm-time $warm_time --seed $seed \
--pert-strength0 $refe_strength --pert-strength1 $pert_strength"

echo "# run on `hostname`" &>> run.resp.log
echo "# command is $command" &>> run.resp.log
$command > output.resp 2> error.resp &
echo "# pid is $!" &>> run.resp.log

