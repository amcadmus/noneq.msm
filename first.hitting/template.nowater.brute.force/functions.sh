function make_gaussian_tables () {
    table_count=0
    table_list=""
    table_scale=""
    rm -f base.info
    for ii in $fht_gaussian_base_position;
    do
	$cwd/tools/dihedral.table/gaussian --max 1.0 --mu $ii --sigma $fht_gaussian_base_sigma -o table_d${table_count}.xvg
	table_list="$table_list $cwd/table_d${table_count}.xvg"
	table_scale="$table_scale 1.0"
	table_count=$(($table_count+1))
	echo $fht_gaussian_base_max >> base.info
    done
    # if test $table_count -ne 0; then
    # 	tmp_s1=" \" ${table_list} \" "
    # 	tmp_s2=" \" ${table_scale} \" "	
    # 	# command="$cwd/tools/dihedral.table/combine.table.1 --input-list ' ${table_list} ' --scale-list ' ${table_scale} ' -o $cwd/$my_dir/effective.table.xvg"
    # 	# echo command is  $command
    # 	$cwd/tools/dihedral.table/combine.table.1 --input-list $tmp_s1 --scale-list $tmp_s2
    # fi
}

function make_cos_tables () {
    table_count=0
    if test $fht_cos_base_number -gt `grep -v \# $fht_cos_base_k_file | wc -l | awk '{print $1}'`; then
	echo "input number of base $fht_cos_base_number is larger than the number of base in file $fht_cos_base_k_file, do nothing";
	exit
    fi
    rm -f base.info
    for ii in `seq 1 $fht_cos_base_number`;
    do
	make_cos_top_line=`grep -v \# $fht_cos_base_k_file | head -n $ii | tail -n 1`
	make_cos_top_gamma=`echo $make_cos_top_line | awk '{print $1}'`
	make_cos_top_kk=`echo $make_cos_top_line | awk '{print $2}'`
	$cwd/tools/dihedral.table/cos --gamma $make_cos_top_gamma -o table_d${table_count}.xvg > /dev/null
	echo $make_cos_top_kk >> base.info
	table_count=$(($table_count+1))
    done
}


function make_gaussian_top () {
    echo '; n-butane'							>  topol.top
    echo '; by Han Wang (han.wang@fu-berlin.de) on 29.3.2013'		>> topol.top
    echo '#include "gromos45a3.ff/forcefield.itp"'				>> topol.top
    echo ''	>> topol.top
    echo '[ moleculetype ]'>> topol.top
    echo '; Name  nrexcl'>> topol.top
    echo 'butane  3'>> topol.top
    echo ''>> topol.top
    echo '[ atoms ]'>> topol.top
    echo '; nr    type    resdnr  resd    atom    cgnr    charge  mass'>> topol.top
    # echo '1       CH3     1       C4      CH3     1       0.0     15.035'>> topol.top
    # echo '2       CH2     1       C4      CH2     1       0.0     14.027'>> topol.top
    # echo '3       CH2     1       C4      CH2     1       0.0     14.027'>> topol.top
    # echo '4       CH3     1       C4      CH3     1       0.0     15.035'>> topol.top
    echo '1       CH3     1       C4      CH3     1       0.0     1.'>> topol.top
    echo '2       CH2     1       C4      CH2     1       0.0     1.'>> topol.top
    echo '3       CH2     1       C4      CH2     1       0.0     1.'>> topol.top
    echo '4       CH3     1       C4      CH3     1       0.0     1.'>> topol.top
    echo ''>> topol.top
    echo '[ bonds ]'>> topol.top
    echo '; ai    aj      funct   param'>> topol.top
    echo '1       2       2       gb_26'>> topol.top
    echo '2       3       2       gb_26'>> topol.top
    echo '3       4       2       gb_26'>> topol.top
    echo ''>> topol.top
    echo '[ pairs ]'>> topol.top
    echo '; ai    aj      funct'>> topol.top
    echo ';1      4       1'>> topol.top
    echo ''>> topol.top
    echo '[ exclusions ]'>> topol.top
    echo '4 1'>> topol.top
    echo ''>> topol.top
    echo '[ angles ]'>> topol.top
    echo '; ai    aj      ak      funct   param'>> topol.top
    echo '1       2       3       2       ga_14'>> topol.top
    echo '2       3       4       2       ga_14'>> topol.top
    echo ''>> topol.top
    echo '[ dihedrals ] '>> topol.top
    echo '; ai    aj      ak      al      funct   tab_num k'>> topol.top
    echo '1       2       3       4       1       gd_17'>> topol.top
    top_item_count=0
    for i in `echo $fht_gaussian_base_position`
    do
	echo "1      2       3       4       8	$top_item_count	$fht_gaussian_base_max" >> topol.top
	top_item_count=$(($top_item_count+1))
    done    
    echo ''>> topol.top
    echo '#include "gromos45a3.ff/spce.itp"'>> topol.top
    echo ''>> topol.top
    echo '[ system ]'>> topol.top
    echo 'Butane in vacumm'>> topol.top
    echo ''>> topol.top
    echo '[ molecules ]'>> topol.top
    echo 'butane  1'>> topol.top
    echo 'SOL     0 '>> topol.top
    echo ''>> topol.top
    echo ''>> topol.top
}


function make_cos_top () {
    echo '; n-butane'							>  topol.top
    echo '; by Han Wang (han.wang@fu-berlin.de) on 29.3.2013'		>> topol.top
    echo '#include "gromos45a3.ff/forcefield.itp"'				>> topol.top
    echo ''	>> topol.top
    echo '[ moleculetype ]'>> topol.top
    echo '; Name  nrexcl'>> topol.top
    echo 'butane  3'>> topol.top
    echo ''>> topol.top
    echo '[ atoms ]'>> topol.top
    echo '; nr    type    resdnr  resd    atom    cgnr    charge  mass'>> topol.top
    # echo '1       CH3     1       C4      CH3     1       0.0     15.035'>> topol.top
    # echo '2       CH2     1       C4      CH2     1       0.0     14.027'>> topol.top
    # echo '3       CH2     1       C4      CH2     1       0.0     14.027'>> topol.top
    # echo '4       CH3     1       C4      CH3     1       0.0     15.035'>> topol.top
    echo '1       CH3     1       C4      CH3     1       0.0     1.'>> topol.top
    echo '2       CH2     1       C4      CH2     1       0.0     1.'>> topol.top
    echo '3       CH2     1       C4      CH2     1       0.0     1.'>> topol.top
    echo '4       CH3     1       C4      CH3     1       0.0     1.'>> topol.top
    echo ''>> topol.top
    echo '[ bonds ]'>> topol.top
    echo '; ai    aj      funct   param'>> topol.top
    echo '1       2       2       gb_26'>> topol.top
    echo '2       3       2       gb_26'>> topol.top
    echo '3       4       2       gb_26'>> topol.top
    echo ''>> topol.top
    echo '[ pairs ]'>> topol.top
    echo '; ai    aj      funct'>> topol.top
    echo ';1      4       1'>> topol.top
    echo ''>> topol.top
    echo '[ exclusions ]'>> topol.top
    echo '4 1'>> topol.top
    echo ''>> topol.top
    echo '[ angles ]'>> topol.top
    echo '; ai    aj      ak      funct   param'>> topol.top
    echo '1       2       3       2       ga_14'>> topol.top
    echo '2       3       4       2       ga_14'>> topol.top
    echo ''>> topol.top
    echo '[ dihedrals ] '>> topol.top
    echo '; ai    aj      ak      al      funct   tab_num k'>> topol.top
    echo '1       2       3       4       1       gd_17'>> topol.top
    top_item_count=0
    for ii in `seq 1 $fht_cos_base_number`;
    do
	make_cos_top_line=`grep -v \# $fht_cos_base_k_file | head -n $ii | tail -n 1`
	make_cos_top_gamma=`echo $make_cos_top_line | awk '{print $1}'`
	make_cos_top_kk=`echo $make_cos_top_line | awk '{print $2}'`
	echo "1      2       3       4       8	$top_item_count	$make_cos_top_kk" >> topol.top
	top_item_count=$(($top_item_count+1))
    done    
    echo ''>> topol.top
    echo '#include "gromos45a3.ff/spce.itp"'>> topol.top
    echo ''>> topol.top
    echo '[ system ]'>> topol.top
    echo 'Butane in vacumm'>> topol.top
    echo ''>> topol.top
    echo '[ molecules ]'>> topol.top
    echo 'butane  1'>> topol.top
    echo 'SOL     0 '>> topol.top
    echo ''>> topol.top
    echo ''>> topol.top
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
    fht_energy_feq_r=`echo "$fht_energy_feq / $fht_dt" | bc -l | cut -d '.' -f 1`
    fht_seed=`date +%s`
    fht_seed=`echo "$fht_seed + $fht_num_conf_use * 10 * $fht_parallel_my_id" | bc`
    sed -e "/^dt/s/=.*/= $fht_dt/g" $file |\
    sed -e "/^integrator/s/=.*/= $fht_integrator/g" |\
    sed -e "/^Tcoupl /s/=.*/= no/g" |\
    sed -e "/^Pcoupl /s/=.*/= $fht_barostat/g" |\
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

