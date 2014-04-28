#!/bin/bash

function make_tables () {
    if test $fht_base_phi_number -gt `grep -v \# $fht_base_phi_k_file | wc -l | awk '{print $1}'`; then
	echo "input number of base $fht_base_phi_number is larger than the number of base in file $fht_base_phi_k_file, do nothing";
	exit
    fi
    rm -f base.phi.info
    table_count=0
    for ii in `seq 1 $fht_base_phi_number`;
    do
	make_top_line=`grep -v \# $fht_base_phi_k_file | head -n $ii | tail -n 1`
	make_top_type=`echo $make_top_line | awk '{print $1}'`
	make_top_gamma=`echo $make_top_line | awk '{print $2}'`
	make_top_kk=`echo $make_top_line | awk '{print $3}'`
	if echo $make_top_type | grep cos &> /dev/null; then
	    # cos base
	    $cwd/tools/dihedral.table/cos --gamma $make_top_gamma -o table_d${table_count}.xvg > /dev/null
	else if echo $make_top_type | grep sin &> /dev/null; then
	    # cos base
	    $cwd/tools/dihedral.table/sin --gamma $make_top_gamma -o table_d${table_count}.xvg > /dev/null
	else
	    echo "unknow type: $make_top_type"
	fi
	fi
	echo $make_top_kk >> base.phi.info
	table_count=$(($table_count+1))
    done

    if test $fht_base_psi_number -gt `grep -v \# $fht_base_psi_k_file | wc -l | awk '{print $1}'`; then
	echo "input number of base $fht_base_psi_number is larger than the number of base in file $fht_base_psi_k_file, do nothing";
	exit
    fi
    rm -f base.psi.info
    for ii in `seq 1 $fht_base_psi_number`;
    do
	make_top_line=`grep -v \# $fht_base_psi_k_file | head -n $ii | tail -n 1`
	make_top_type=`echo $make_top_line | awk '{print $1}'`
	make_top_gamma=`echo $make_top_line | awk '{print $2}'`
	make_top_kk=`echo $make_top_line | awk '{print $3}'`
	if echo $make_top_type | grep cos &> /dev/null; then
	    # cos base
	    $cwd/tools/dihedral.table/cos --gamma $make_top_gamma -o table_d${table_count}.xvg > /dev/null
	else if echo $make_top_type | grep sin &> /dev/null; then
	    # cos base
	    $cwd/tools/dihedral.table/sin --gamma $make_top_gamma -o table_d${table_count}.xvg > /dev/null
	else
	    echo "unknow type: $make_top_type"
	fi
	fi
	echo $make_top_kk >> base.psi.info
	table_count=$(($table_count+1))
    done
}


function make_top () {
    echo '; alanine dipeptide'							>  topol.top
    echo '#include "amber99sb-ildn.ff/forcefield.itp"'				>> topol.top
    echo ''									>> topol.top
    echo '#include "alanine.itp"'						>> topol.top
    echo ''									>> topol.top
    echo '#include "amber99sb-ildn.ff/tip3p.itp"'				>> topol.top
    echo ''									>> topol.top
    echo '[ system ]'								>> topol.top
    echo 'Alanine in vacumm'							>> topol.top
    echo ''									>> topol.top
    echo '[ molecules ]'							>> topol.top
    echo 'Protein  1'								>> topol.top
    echo 'SOL      0'								>> topol.top
    echo ''									>> topol.top
    echo ''									>> topol.top
    
    echo '[ dihedrals ]'							>> alanine.itp
    echo ''									>> alanine.itp

    top_item_count=0
    for ii in `seq 1 $fht_base_phi_number`;
    do
	make_top_line=`grep -v \# $fht_base_phi_k_file | head -n $ii | tail -n 1`
	make_top_kk=`echo $make_top_line | awk '{print $3}'`
	echo "5      7       9      15       8	$top_item_count	$make_top_kk"	>> alanine.itp
	top_item_count=$(($top_item_count+1))
    done    
    for ii in `seq 1 $fht_base_psi_number`;
    do
	make_top_line=`grep -v \# $fht_base_psi_k_file | head -n $ii | tail -n 1`
	make_top_kk=`echo $make_top_line | awk '{print $3}'`
	echo "7      9      15      17       8	$top_item_count	$make_top_kk"	>> alanine.itp
	top_item_count=$(($top_item_count+1))
    done    

    echo ''									>> alanine.itp
}


function set_parameters_fht () {
    file=$1
    if ! `grep userreal5 $file &>/dev/null`; then
	echo "userreal5 = 1.0" >> $file
    fi
    if ! `grep userreal6 $file &>/dev/null`; then
	echo "userreal6 = 0.0" >> $file
    fi
    if ! `grep userreal7 $file &>/dev/null`; then
	echo "userreal7 = 0.0" >> $file
    fi
    if ! `grep nstcalcenergy $file &>/dev/null`; then
	echo "nstcalcenergy = 1.0" >> $file
    fi
    fht_time=`echo "$fht_stop_time + $fht_frame_feq" | bc -l`
    fht_nstep=`echo "$fht_time / $fht_dt" | bc -l | cut -d '.' -f 1`
    fht_xtcout_feq=`echo "$fht_frame_feq / $fht_dt" | bc -l | cut -d '.' -f 1`
    fht_xvout_feq=$fht_xtcout_feq
    fht_xvout_feq=0
    fht_energy_feq_r=`echo "$fht_energy_feq / $fht_dt" | bc -l | cut -d '.' -f 1`
    fht_seed=`date +%s`
    fht_seed=`echo "$fht_seed + $fht_num_conf_use * 10 * $fht_parallel_my_id" | bc`
    sed -e "/^dt/s/=.*/= $fht_dt/g" $file |\
    sed -e "/^integrator/s/=.*/= $fht_integrator/g" |\
    sed -e "/^Tcoupl /s/=.*/= no/g" |\
    sed -e "/^Pcoupl /s/=.*/= $fht_barostat/g" |\
    sed -e "/^ref_t/s/=.*/= $fht_T/g" |\
    sed -e "/^tau_t/s/=.*/= $fht_taut/g" |\
    sed -e "/^tau_p/s/=.*/= $fht_taup/g" |\
    sed -e "/^nstep/s/=.*/= $fht_nstep/g" |\
    sed -e "/^nstxout/s/=.*/= $fht_xvout_feq/g" |\
    sed -e "/^nstvout/s/=.*/= $fht_xvout_feq/g" |\
    sed -e "/^nstfout/s/=.*/= 0/g" |\
    sed -e "/^nstenergy/s/=.*/= $fht_energy_feq_r/g" |\
    sed -e "/^nstcalcenergy/s/=.*/= $fht_energy_feq_r/g" |\
    # sed -e "/^userreal1/s/=.*/= $fht_noSdRange/g" |\
    # sed -e "/^userreal2/s/=.*/= $fht_ele_rescale_start/g" |\
    # sed -e "/^userreal3/s/=.*/= $fht_ele_rescale_end/g" |\
    # sed -e "/^userreal4/s/=.*/= $fht_rescale/g" |\
    # sed -e "/^userreal5/s/=.*/= $fht_ele_rescale_scale0/g" |\
    sed -e "/^userreal6/s/=.*/= $fht_meta_low/g" |\
    sed -e "/^userreal7/s/=.*/= $fht_meta_up/g" |\
    sed -e "/^E-x /s/=.*/= 0 0.0 0.0/g" |\
    sed -e "/^ld-seed/s/=.*/= $fht_seed/g" |\
    sed -e "/^gen_vel /s/=.*/= no/g" |\
    sed -e "/^gen-vel /s/=.*/= no/g" |\
    sed -e "/^nstxtcout/s/=.*/= $fht_xtcout_feq/g" > tmptmptmp.mdp
    mv -f tmptmptmp.mdp $file
}

