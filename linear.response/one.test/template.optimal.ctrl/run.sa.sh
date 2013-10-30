#!/bin/bash

source parameters.sh

make -C tools/simulator makedir &> make.log
make -C tools/simulator -j8 &>> make.log

rm -f run.resp.log

print_noneq_step=`printf %.2f $noneq_check_feq`
print_noneq_time=`printf %.1f $noneq_time`
print_quench_time=`printf %.1f $quench_time`
print_x0=`printf %.1f $x0`
print_x1=`printf %.1f $x1`
print_v0=`printf %.1f $v0`
print_v1=`printf %.1f $v1`
print_nx=`printf %d $nx`
print_nv=`printf %d $nv`

if test -f $init_ctr; then
    init_ctr_opt=" --init $init_ctr"
else
    init_ctr_opt=" "
fi

command="`pwd`/tools/simulator/$project_name.ctr.sa
--dt $dt \
    --nst $nst \
    -p 1000000 \
    --gamma $gamma \
    --temperature $temperature \
    $command_line_param_print \
    --branch-feq $branch_feq \
    --time-resolution $time_resolution \
    --beta $beta \
    --noneq-check-feq $noneq_check_feq \
    --noneq-time $noneq_time \
    --seed $seed \
    --order $resp_order \
    $init_ctr_opt \
    --saNst $saNst \
    --saTmax $saTmax \
    --saSigma $saSigma \
    --saChangeMin $saChangeMin \
    --pert-strength0 $pert_strength"

echo "# run on `hostname`" &>> run.ctr.log
echo "# command is $command" &>> run.ctr.log
mpiexec -np 8 $command > output.ctr 2> error.ctr 
echo "# pid is $!" &>> run.ctr.log
