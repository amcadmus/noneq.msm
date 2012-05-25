#!/bin/bash

source parameters.sh

make -C tools/simulator makedir &> make.log
make -C tools/simulator -j8 &>> make.log

rm -f run.log

print_corr_step=`printf %.3f $corr_step`
print_corr_time=`printf %.1f $corr_time`
if test ! -d $saved_corr_dir; then
    echo "# make dir $saved_corr_dir"
    mkdir -p $saved_corr_dir
fi

if echo $load_saved_corr | grep yes &> /dev/null; then
    print_equi_time=`printf %.2f $saved_corr_equi_time`
    saved_corr_name=save.corr.step$print_corr_step.time$print_corr_time.equiT$print_equi_time
    if test -f $saved_corr_dir/$saved_corr_name; then
	echo "# existing saved corr $saved_corr_dir/$saved_corr_name, load it"
	command="`pwd`/tools/simulator/timeCorr \
--dt $dt --nst $nst -p 100000 --gamma $gamma --temperature $temperature \
--double-well-k $double_well_k --double-well-a $double_well_a \
--branch-feq $branch_feq \
--noneq-check-feq $noneq_check_feq --noneq-time $noneq_time \
--quench-temperature $quench_temperature --quench-time $quench_time \
--corr-time $corr_time --corr-step $corr_step \
--load-corr $saved_corr_dir/$saved_corr_name \
--warm-time $warm_time --seed $seed --pert-strength $pert_strength"
    else
	echo "# no saved corr $saved_corr_dir/$saved_corr_name, exit"
	exit
    fi
else
    equi_time=`echo "$dt * $nst" | octave | grep ans | cut -d '=' -f 2`
    print_equi_time=`printf %.2f $equi_time`
    saved_corr_name=save.corr.step$print_corr_step.time$print_corr_time.equiT$print_equi_time
    echo "# save corr to $saved_corr_dir/$saved_corr_name"
    command="`pwd`/tools/simulator/timeCorr \
--dt $dt --nst $nst -p 100000 --gamma $gamma --temperature $temperature \
--double-well-k $double_well_k --double-well-a $double_well_a \
--branch-feq $branch_feq \
--noneq-check-feq $noneq_check_feq --noneq-time $noneq_time \
--quench-temperature $quench_temperature --quench-time $quench_time \
--corr-time $corr_time --corr-step $corr_step \
--save-corr $saved_corr_dir/$saved_corr_name \
--warm-time $warm_time --seed $seed --pert-strength $pert_strength"
fi

echo "# command is $command" &>> run.log
$command > output.timeCorr &
echo "# pid is $!" &>> run.log

command="`pwd`/tools/simulator/noneq \
--dt $dt --nst $nst -p 100000 --gamma $gamma --temperature $temperature \
--double-well-k $double_well_k --double-well-a $double_well_a \
--branch-feq $branch_feq \
--noneq-check-feq $noneq_check_feq --noneq-time $noneq_time \
--quench-temperature $quench_temperature --quench-time $quench_time \
--warm-time $warm_time --seed $seed --pert-strength $pert_strength"

echo "# command is $command" &>> run.log
$command > output.noneq &
echo "# pid is $!" &>> run.log

