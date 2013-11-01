#!/bin/bash

#!/bin/bash

source env.sh
source parameters.sh
source functions.sh

make -C tools/dihedral.table/ makedir
make -C tools/dihedral.table/ -j8

cwd=`pwd`
# prepare files
if test ! -d $equi_result ; then
    echo "no dir $equi_result, exit"
    exit
fi
if test ! -f $equi_result/equi.frame ; then
    echo "no file $equi_result/equi.frame, exit"
    exit
fi
#rm -f angle.name
#rm -f gxs.name
#rm -f success.dir.name
touch success.dir.name

targets=`awk '{print $1}' $equi_result/equi.frame | head -n $fht_num_conf_use`

fht_main_dir=result.fhts
for i in $targets;
do
    count=$i
    runid=`echo "$count % $fht_parallel_num_pro" | bc`
    test $runid -ne $fht_parallel_my_id && continue
    test ! -d $fht_main_dir && mkdir -p $fht_main_dir
    my_dir=$fht_main_dir/fht.$count
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
    cp $equi_result/conf.gro	$my_dir
    cp $equi_result/grompp.mdp	$my_dir
    cp $equi_result/angle.ndx	$my_dir
    cp $equi_result/topol.top	$my_dir
    cd $my_dir
    make_top
    count=0
    for ii in $fht_gaussian_base_position;
    do
	./tools/dihedral.table/gaussian --max $fht_gaussian_base_max --mu $ii --sigma $fht_gaussian_base_sigma -o table_d${count}.xvg
	count=$(($count+1))
    done
    set_parameters_fht grompp.mdp
    
    start_time=`grep ^$count $equi_result/equi.frame | awk '{print $2}'`
    echo "# run with command `which grompp`"
    $grompp_command -t $equi_result/traj.trr -time $start_time &> run.log
    if [ $? -ne 0 ]; then
	echo "failed at grompp exit"; exit
    fi
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

    tmpid=`echo "$count - $fht_parallel_num_pro" | bc -l`
    echo "tmpid is $tmpid"
    if [ $tmpid -lt $fht_parallel_num_pro ]; then
	cp -a ..//fht.$count ..//backup.fht.$count
    fi
    rm -f traj.xtc traj.trr state*.cpt topol.tpr conf.gro index.ndx angle.log md.log genbox.log mdout.mdp protein.gro run.log tablep.xvg table.xvg grompp.mdp topol.top angdist.xvg 
    
    cd $cwd
    # echo "$my_dir/angle.dat" >> angle.name
    # echo "$my_dir/gxs.out" >> gxs.name
    echo "$my_dir" >> success.dir.name
#    sync
done

