#!/bin/bash

source env.sh
source parameters.sh
source functions.sh

make -C tools/analysis/ makedir
make -C tools/analysis/ -j8
make -C tools/dihedral.table/ makedir
make -C tools/dihedral.table/ -j8

cwd=`pwd`
# prepare files
if test ! -d $fht_equi_dir ; then
    echo "no dir $fht_equi_dir, exit"
    exit
fi
if test ! -f $fht_equi_frame_name ; then
    echo "no file $fht_equi_frame_name, exit"
    exit
fi
#rm -f angle.name
#rm -f gxs.name
#rm -f success.dir.name
touch success.dir.name
touch last.time.record

make_tables
rm -f $cwd/coreset*dat
cp $fht_coreset_data $cwd/
#targets=`awk '{print $1}' $fht_equi_dir/$fht_equi_frame_name | head -n $fht_num_conf_use`
targets=`seq 0 $(($fht_num_conf_use-1))`

fht_main_dir=result.traj
for i in $targets;
do
    count=`printf %08d $i`
    runid=`echo "$count % $fht_parallel_num_pro" | bc`
    test $runid -ne $fht_parallel_my_id && continue
    test ! -d $fht_main_dir && mkdir -p $fht_main_dir
    my_dir=$fht_main_dir/fht.$count
    if grep $my_dir success.dir.name &> /dev/null; then
	echo "# existing successful dir $my_dir, continue"
	continue
    fi
    if test -d $my_dir; then
	echo "# existing dir $my_dir, but failed, remrove"
	rm -fr $my_dir
    fi
    echo "# doing in dir $my_dir"
    mkdir -p $my_dir
    cp -L $fht_equi_dir/conf.gro	 $my_dir
    cp -L $fht_equi_dir/grompp.mdp	 $my_dir
    cp -L $fht_equi_dir/angle.ndx	 $my_dir
    cp -L $fht_equi_dir/alanine.itp	 $my_dir
    cp -L $cwd/$fht_base_phi_k_file	 $my_dir
    cp -L $cwd/$fht_base_psi_k_file	 $my_dir
    cp -L $fht_equi_dir/topol.top	 $my_dir
    cp $fht_coreset_data $my_dir/coreset.dat
    cd $my_dir
    make_top
#    make_tables
    test -f $cwd/table_d0.xvg && cp $cwd/table_d*xvg .
    set_parameters_fht grompp.mdp

    nframe_equi=`cat $fht_equi_frame_count`
    equi_code=`echo "($count % $nframe_equi) + 1" | bc`
    start_time=`echo "$equi_code - 0.5" | bc -l`
    echo "# run with command `which grompp`, starting time $start_time"
    $grompp_command -t $fht_equi_frame_traj -time $start_time &> run.log
    if [ $? -ne 0 ]; then
	echo "failed at grompp exit"; exit
    fi
    echo "# run with command `which mdrun`"
    $mdrun_command &> run.log
    if [ $? -ne 0 ]; then
	echo "failed at mdrun exit"; exit
    fi

    # echo 2 2 | trjconv -center -pbc whole &> run.log
    # if [ $? -ne 0 ]; then
    # 	echo "failed at trjconv exit"; exit
    # fi
    # mv -f trajout.xtc butane.xtc
    trajfile=traj.xtc
    if test -f traj.trr; then
	trajfile=traj.trr;
    fi
    echo 0 | g_angle -f $trajfile -n angle.ndx -type dihedral -od angdist.xvg  -ov angaver.xvg -xvg none &> run.log
    if [ $? -ne 0 ]; then
	echo "failed at g_angle exit"; exit
    fi
    mv -f angdist.xvg angdist.phi.xvg
    mv -f angaver.xvg angaver.phi.xvg
    echo 1 | g_angle -f $trajfile -n angle.ndx -type dihedral -od angdist.xvg  -ov angaver.xvg -xvg none &> run.log
    if [ $? -ne 0 ]; then
	echo "failed at g_angle exit"; exit
    fi
    mv -f angdist.xvg angdist.psi.xvg
    mv -f angaver.xvg angaver.psi.xvg
    
    traj_last_time=`tail -n 1 angaver.phi.xvg | awk '{print $1}'`
    traj_last_angle_phi=`tail -n 1 angaver.phi.xvg | awk '{print $2}'`
    traj_last_angle_psi=`tail -n 1 angaver.psi.xvg | awk '{print $2}'`
    last_angle_ind=`$cwd/tools/analysis/angle.ind --input-angle-phi $traj_last_angle_phi --input-angle-psi $traj_last_angle_psi --input-coreset-file $fht_coreset_data`
    bool_hit_set=`echo "$last_angle_ind == $fht_coreset_target" | bc `
    bool_time_out=`echo "$last_angle_ind == $fht_coreset_notin" | bc `

    if [ $bool_time_out -eq 1 ]; then
	echo "WARNING: find traj out of time"
    fi

    tmpid=`echo "$count - $fht_parallel_num_pro" | bc -l`
    echo "tmpid is $tmpid"
    if [ $tmpid -lt 0 ]; then
	if [ $count -eq 0 ] ; then
	    cp -a ..//fht.$count ..//backup.fht.$count
	fi
    fi

    # hit meta, remove useless files
    if test $bool_hit_set -eq 1; then
	rm -f traj.xtc traj.trr state*.cpt topol.tpr conf.gro index.ndx angle.log md.log genbox.log mdout.mdp protein.gro run.log tablep.xvg table.xvg grompp.mdp topol.top angdist.*.xvg angle.ndx gxs.out table_d*xvg ener.edr cos.k.in confout.gro coreset.dat  base.k.phi base.k.psi alanine.itp &
    fi
    
    cd $cwd

    # does not hit, remove the dir!
    if test $bool_hit_set -eq 0; then
	rm -fr $my_dir &
    fi
	
    # echo "$my_dir/angle.dat" >> angle.name
    # echo "$my_dir/gxs.out" >> gxs.name
    echo "$my_dir" >> success.dir.name
    echo "$my_dir $traj_last_time $traj_last_angle_phi $traj_last_angle_psi $last_angle_ind" >> last.time.record
#    sync
done

