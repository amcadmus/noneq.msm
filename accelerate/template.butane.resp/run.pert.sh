#!/bin/bash

source env.sh
source parameters.sh
source functions.sh

make -C tools/angles/ makedir
# make -C tools/angles/ clean
make -C tools/angles/ -j8
make -C tools/parse.top/ makedir
# make -C tools/parse.top/ clean
make -C tools/parse.top/ -j8
make -C tools/potential/ makedir
# make -C tools/potential/ clean
make -C tools/potential/ -j8

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
#rm -f angle.name
#rm -f gxs.name
#rm -f success.dir.name
touch success.dir.name

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
	if grep $my_dir success.dir.name &> /dev/null; then
	    echo "# existing successful dir $my_dir, continue"
	    continue
	else
	    echo "# existing dir $my_dir, but failed, remrove"
	    rm -fr $my_dir
	fi
    fi
    echo "# doing in dir $my_dir"
    mkdir $my_dir
    cp $pert_equi_result/conf.gro	$my_dir
    cp $pert_equi_result/grompp.mdp	$my_dir
    cp $pert_equi_result/angle.ndx	$my_dir
    cp $pert_equi_result/topol.top	$my_dir
    cd $my_dir    
    rm -f run.log
    ../../tools/potential/tab.coulomb --scale $pert_rescale --scale-ele-start 1000 --scale-ele-end 1001 --scale-ele-scale-0 1.0 --scale-vdw-scale-0 $pert_vdw_rescale_scale0 --scale-vdw-start $pert_vdw_rescale_start --scale-vdw-end $pert_vdw_rescale_end --rc $pert_rcut --r-smooth $pert_rsmooth --erf $pert_erf -o table.xvg
    ../../tools/potential/tab.coulomb --scale $pert_rescale --scale-ele-start 1000 --scale-ele-end 1001 --scale-ele-scale-0 1.0 --scale-vdw-scale-0 $pert_vdw_rescale_scale0 --scale-vdw-start $pert_vdw_rescale_start --scale-vdw-end $pert_vdw_rescale_end --rc $pert_rcut --r-smooth $pert_rcut --erf $pert_erf -o tablep.xvg
    if [ $? -ne 0 ]; then
	echo "failed at tab.rf, exit"; exit
    fi
    ../../tools/parse.top/rescale.param -f topol.top -b $pert_bond_rescale -a $pert_angle_rescale --scal-angle-bond $pert_angle_bond_rescale -d $pert_rescale --scal-cmap $pert_rescale -o out.top
    mv -f out.top topol.top
    set_parameters_pert grompp.mdp
    start_time=`grep ^$count $pert_equi_result/equi.frame | awk '{print $2}'`
    echo "# run with command `which grompp`"
    $grompp_command -t $pert_equi_result/traj.trr -time $start_time &> run.log
    if [ $? -ne 0 ]; then
	echo "failed at grompp exit"; exit
    fi
    # echo "# run with command `which g_pme_error`"
    # g_pme_error -tune -self 0.0001 &> run.log
    # if [ $? -ne 0 ]; then
    # 	echo "failed at g_pme_error exit"; exit
    # fi
    # mv -f tuned.tpr topol.tpr
    echo "# run with command `which mdrun`"
    $mdrun_command &> run.log
    if [ $? -ne 0 ]; then
	echo "failed at mdrun exit"; exit
    fi

    echo 2 2 | trjconv -center -pbc whole
    mv -f trajout.xtc butane.xtc
    # $cwd/tools/angles/evolve -f alanine.xtc -s angle.dat &> angle.log
    # if [ $? -ne 0 ]; then
    # 	echo "failed at evolve exit"; exit
    # fi
    g_angle -type dihedral -od angdist.xvg  -ov angaver.xvg -xvg none &> run.log

    tmpid=`echo "$count - $pert_parallel_num_pro" | bc -l`
    echo "tmpid is $tmpid"
    if [ $tmpid -lt $pert_parallel_num_pro ]; then
	cp -a ..//pert.$count ..//backup.pert.$count
    fi
    rm -f traj.xtc traj.trr state*.cpt topol.tpr conf.gro index.ndx angle.log md.log genbox.log mdout.mdp protein.gro run.log tablep.xvg table.xvg grompp.mdp topol.top angdist.xvg 
    
    cd $cwd
    echo "$my_dir/angle.dat" >> angle.name
    echo "$my_dir/gxs.out" >> gxs.name
    echo "$my_dir" >> success.dir.name
#    sync
done

