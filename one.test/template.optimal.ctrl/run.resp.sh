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

equi_time=`echo "$dt * $nst" | octave | grep ans | cut -d '=' -f 2`
print_equi_time=`printf %.2f $equi_time`
print_quench_time=`printf %.2f $quench_time`
print_refStr=`printf %.2f $refe_strength`

if test ! -d $saved_resp_dir; then
    echo "# make dir $saved_resp_dir"
    mkdir -p $saved_resp_dir
fi

if echo $load_saved_resp | grep yes &> /dev/null; then
    saved_resp_name=$project_name.$save_corr_param_note.save.resp.refStr$print_refStr.step$print_noneq_step.time$print_noneq_time.quench$print_quench_time.equiT$print_equi_time.grid$print_x0.$print_x1.$print_nx.$print_v0.$print_v1.$print_nv
    if test -f $saved_resp_dir/$saved_resp_name; then
	echo "# existing saved resp $saved_resp_dir/$saved_resp_name, load it"
	command="`pwd`/tools/simulator/$project_name.noneqResp \
--dt $dt --nst $nst -p 1000000 --gamma $gamma --temperature $temperature \
$command_line_param_print \
--branch-feq $branch_feq \
--noneq-check-feq $noneq_check_feq --noneq-time $noneq_time \
--quench-temperature $quench_temperature --quench-time $quench_time \
--load-corr $saved_resp_dir/$saved_resp_name \
--x-low $x0 --x-up $x1 --v-low $v0 --v-up $v1 --x-grid $nx --v-grid $nv \
--warm-time $warm_time --seed $seed \
--order $resp_order \
--pert-strength0 $refe_strength --pert-strength1 $pert_strength"
    else
	echo "# no saved resp $saved_resp_dir/$saved_resp_name, exit"
	exit
    fi
else
    saved_resp_name=$project_name.$save_corr_param_note.save.resp.refStr$print_refStr.step$print_noneq_step.time$print_noneq_time.quench$print_quench_time.equiT$print_equi_time.grid$print_x0.$print_x1.$print_nx.$print_v0.$print_v1.$print_nv
    echo "# save resp to $saved_resp_dir/$saved_resp_name"
    command="`pwd`/tools/simulator/$project_name.noneqResp \
--dt $dt --nst $nst -p 1000000 --gamma $gamma --temperature $temperature \
$command_line_param_print \
--branch-feq $branch_feq \
--noneq-check-feq $noneq_check_feq --noneq-time $noneq_time \
--quench-temperature $quench_temperature --quench-time $quench_time \
--save-corr $saved_resp_dir/$saved_resp_name \
--x-low $x0 --x-up $x1 --v-low $v0 --v-up $v1 --x-grid $nx --v-grid $nv \
--warm-time $warm_time --seed $seed \
--order $resp_order \
--pert-strength0 $refe_strength --pert-strength1 $pert_strength"
fi

# command="`pwd`/tools/simulator/$project_name.noneqResp \
# --dt $dt --nst $nst -p 100000 --gamma $gamma --temperature $temperature \
# $command_line_param_print \
# --branch-feq $branch_feq \
# --noneq-check-feq $noneq_check_feq --noneq-time $noneq_time \
# --quench-temperature $quench_temperature --quench-time $quench_time \
# --x-low $x0 --x-up $x1 --v-low $v0 --v-up $v1 --x-grid $nx --v-grid $nv \
# --warm-time $warm_time --seed $seed \
# --pert-strength0 $refe_strength --pert-strength1 $pert_strength"

echo "# run on `hostname`" &>> run.resp.log
echo "# command is $command" &>> run.resp.log
$command > output.resp 2> error.resp &
echo "# pid is $!" &>> run.resp.log

