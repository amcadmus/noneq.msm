#!/bin/bash

source parameters.sh
source functions.sh

make -C tools/angles/ makedir
make -C tools/angles/ -j4

cwd=`pwd`
# prepare files
if test ! -d $pert_equi_result ; then
    echo "no dir $pert_equi_result, exit"
    exit
fi
if test ! -f $pert_equi_result/equi.frame ; then
    echo "no file $pert_equi_result/equi.frame, exit"
    exit
fi
if test ! -d $gro_dir; then
    echo "no dir $gro_dir, exit"
    exit
fi
#rm -f angle.name
#rm -f gxs.name

targets=`awk '{print $1}' $pert_equi_result/equi.frame | head -n $pert_num_conf_use`

pert_main_dir=result.perts
for i in $targets;
do
    count=$i
    runid=`echo "$count % $pert_parallel_num_pro" | bc`
    test $runid -ne $pert_parallel_my_id && continue
    test ! -d $pert_main_dir && mkdir -p $pert_main_dir
    my_dir=$pert_main_dir/pert.$count
    if test -d $my_dir; then
	echo "existing dir $my_dir, continue"
	continue
	rm -fr $my_dir
    fi
    cp -a $gro_dir $my_dir
    cd $my_dir
    rm -f run.log
    set_parameters_pert grompp.mdp
    start_time=`grep $count $cwd/$pert_equi_result/equi.frame | awk '{print $2}'`
    echo "# run with command `which grompp`" &>> run.log
    $grompp_command -t $cwd/$pert_equi_result/traj.trr -time $start_time &>> run.log
    echo "# run with command `which mdrun`" &>> run.log
    $mdrun_command &>> run.log

    echo 3 14 | trjconv -center -pbc whole
    mv -f trajout.xtc alanine.xtc
    $cwd/tools/angles/evolve -f alanine.xtc -s angle.dat &> angle.log
    rm -f traj.xtc traj.trr state*.cpt topol.tpr conf.gro index.ndx 
    
    cd $cwd
    
    echo "$my_dir/angle.dat" >> angle.name    
    echo "$my_dir/gxs.out" >> gxs.name    
done

