#!/bin/bash

function set_parameters_equi () {
    file=$1
    equi_nstep_warm=`echo "$equi_warm_time/$equi_dt" | bc -l | cut -d '.' -f 1`
    equi_nstep_count=`echo "$equi_frame_feq * $equi_num_frame / $equi_dt" | bc -l | cut -d '.' -f 1`
    equi_nstep=`echo "$equi_nstep_warm + $equi_nstep_count" | bc`
    equi_xvout_feq=`echo "$equi_frame_feq / $equi_dt" | bc -l | cut -d '.' -f 1`
    sed -e "/^integrator/s/=.*/= sd/g" $file |\
    sed -e "/^dt/s/=.*/= $equi_dt/g" |\
    sed -e "/^Tcoupl/s/=.*/= no/g" |\
    sed -e "/^nstep/s/=.*/= $equi_nstep/g" |\
    sed -e "/^nstxout/s/=.*/= $equi_xvout_feq/g" |\
    sed -e "/^nstvout/s/=.*/= $equi_xvout_feq/g" |\
    sed -e "/^ld-seed/s/=.*/= $equi_seed/g" |\
    sed -e "/^gen_vel/s/=.*/= no/g" |\
    sed -e "/^tau_t/s/=.*/= $equi_taut/g" |\
    sed -e "/^nstxtcout/s/=.*/= $equi_xvout_feq/g" > tmptmptmp.mdp
    mv -f tmptmptmp.mdp $file
}

# function set_parameters_long_equi () {
#     file=$1
#     long_equi_nstep_warm=`echo "$long_equi_warm_time/$long_equi_dt" | bc -l | cut -d '.' -f 1`
#     long_equi_nstep_count=`echo "$long_equi_frame_feq * $long_equi_num_frame / $long_equi_dt" | bc -l | cut -d '.' -f 1`
#     long_equi_nstep=`echo "$long_equi_nstep_warm + $long_equi_nstep_count" | bc`
#     long_equi_xvout_feq=`echo "$long_equi_frame_feq / $long_equi_dt" | bc -l | cut -d '.' -f 1`
#     sed -e "/^integrator/s/=.*/= sd/g" $file |\
#     sed -e "/^dt/s/=.*/= $long_equi_dt/g" |\
#     sed -e "/^nstep/s/=.*/= $long_equi_nstep/g" |\
#     sed -e "/^nstxout/s/=.*/= 0/g" |\
#     sed -e "/^nstvout/s/=.*/= 0/g" |\
#     sed -e "/^ld-seed/s/=.*/= $long_equi_seed/g" |\
#     sed -e "/^tau_t/s/=.*/= $long_equi_taut/g" |\
#     sed -e "/^nstxtcout/s/=.*/= $long_equi_xvout_feq/g" > tmptmptmp.mdp
#     mv -f tmptmptmp.mdp $file
# }

function set_parameters_pert () {
    file=$1
    pert_nstep=`echo "$pert_time / $pert_dt" | bc -l | cut -d '.' -f 1`
    pert_xtcout_feq=`echo "$pert_frame_feq / $pert_dt" | bc -l | cut -d '.' -f 1`
    pert_xvout_feq=$pert_xtcout_feq
    pert_energy_feq=$pert_xtcout_feq
    pert_strength_r=`echo "$pert_strength * $pert_rescale" | bc -l`
    pert_warm_time_r=`echo "$pert_warm_time / $pert_rescale" | bc -l`
#    pert_epsilon_r=`echo "1. / $pert_rescale" | bc -l`
    pert_taut_r=`echo "$pert_taut / $pert_rescale" | bc -l`
    sed -e "/^dt/s/=.*/= $pert_dt/g" $file |\
    sed -e "/^integrator/s/=.*/= $pert_integrator/g" |\
    sed -e "/^Tcoupl /s/=.*/= no/g" |\
    sed -e "/^Pcoupl /s/=.*/= $pert_barostat/g" |\
    sed -e "/^tau_t/s/=.*/= $pert_taut_r/g" |\
    sed -e "/^tau_p/s/=.*/= $pert_taup/g" |\
    sed -e "/^nstep/s/=.*/= $pert_nstep/g" |\
    sed -e "/^nstxout/s/=.*/= $pert_xvout_feq/g" |\
    sed -e "/^nstvout/s/=.*/= $pert_xvout_feq/g" |\
    sed -e "/^nstfout/s/=.*/= 0/g" |\
    sed -e "/^nstenergy/s/=.*/= $pert_energy_feq/g" |\
    sed -e "/^userreal1/s/=.*/= $pert_noSdRange/g" |\
    sed -e "/^userreal2/s/=.*/= $pert_ele_rescale_start/g" |\
    sed -e "/^userreal3/s/=.*/= $pert_ele_rescale_end/g" |\
    sed -e "/^userreal4/s/=.*/= $pert_rescale/g" |\
    sed -e "/^rlist /s/=.*/= $pert_rcut/g" |\
    sed -e "/^rcoulomb /s/=.*/= $pert_rcut/g" |\
    sed -e "/^rvdw /s/=.*/= $pert_rcut/g" |\
#    sed -e "/^epsilon_r /s/=.*/= $pert_epsilon_r/g" |\
    sed -e "/^epsilon_rf /s/=.*/= $pert_erf/g" |\
    sed -e "/^E-x /s/=.*/= 1 $pert_strength_r 0.0/g" |\
    sed -e "/^ld-seed/s/=.*/= `date +%s`/g" |\
    sed -e "/^gen_vel /s/=.*/= no/g" |\
    sed -e "/^gen-vel /s/=.*/= no/g" |\
    sed -e "/^nstxtcout/s/=.*/= $pert_xtcout_feq/g" > tmptmptmp.mdp
    mv -f tmptmptmp.mdp $file
    if test $pert_mode -eq 1; then
	sed -e "/^E-xt /s/=.*/= 1 $pert_warm_time_r 0.0/g" $file > tmptmptmp.mdp
    else if test $pert_mode -eq 2; then
	pi=3.14159265359
	my_phi=`echo "$pert_phi / 180.0 * $pi" | bc -l`
	my_T=`echo "2. * $pi / $pert_warm_time_r" | bc -l`
	sed -e "/^E-xt /s/=.*/= 2 $my_T $my_phi $pert_shift 0.0/g" $file > tmptmptmp.mdp
    else
	echo "wrong pert mode: $pert_mode"
	exit
    fi
    fi
    mv -f tmptmptmp.mdp $file
}

function set_interaction () {
    file=$1    
    pert_rescale_sqrt=`echo "sqrt($pert_rescale)" | bc -l`
    param_nb_ch2_c06_r=`echo "$param_nb_ch2_c06 * $pert_rescale" | bc -l`
    param_nb_ch3_c06_r=`echo "$param_nb_ch3_c06 * $pert_rescale" | bc -l`
    param_nb_ow_c06_r=`echo "$param_nb_ow_c06 * $pert_rescale" | bc -l`
    param_nb_ch2_ow_c06_r=`echo "$param_nb_ch2_ow_c06 * $pert_rescale" | bc -l`
    param_nb_ch3_ow_c06_r=`echo "$param_nb_ch3_ow_c06 * $pert_rescale" | bc -l`
    param_nb_ch3_ch2_c06_r=`echo "$param_nb_ch3_ch2_c06 * $pert_rescale" | bc -l`
    param_pair_ow_ow_c06_r=`echo "$param_pair_ow_ow_c06 * $pert_rescale" | bc -l`
    param_pair_ch2_ch2_c06_r=`echo "$param_pair_ch2_ch2_c06 * $pert_rescale" | bc -l`
    param_pair_ch3_ch3_c06_r=`echo "$param_pair_ch3_ch3_c06 * $pert_rescale" | bc -l`
    param_pair_ch2_ow_c06_r=`echo "$param_pair_ch2_ow_c06 * $pert_rescale" | bc -l`
    param_pair_ch3_ow_c06_r=`echo "$param_pair_ch3_ow_c06 * $pert_rescale" | bc -l`
    param_pair_ch3_ch2_c06_r=`echo "$param_pair_ch3_ch2_c06 * $pert_rescale" | bc -l`
    param_charge_ch2_r=`echo "$param_charge_ch2 * $pert_rescale_sqrt" | bc -l`
    param_charge_ch3a_r=`echo "$param_charge_ch3a * $pert_rescale_sqrt" | bc -l`
    param_charge_ch3b_r=`echo "$param_charge_ch3b * $pert_rescale_sqrt" | bc -l`
    param_charge_ow_r=`echo "$param_charge_ow * $pert_rescale_sqrt" | bc -l`
    param_charge_h_r=`echo "$param_charge_h * $pert_rescale_sqrt" | bc -l`
    param_dihedral_c_c_c_c_k=`echo "$param_dihedral_c_c_c_c_k * $pert_rescale" | bc -l`
    param_bond_c_c_k_r=`echo "$param_bond_c_c_k * $pert_bond_rescale" | bc -l`
    param_bond_ow_h_k_r=`echo "$param_bond_ow_h_k * $pert_bond_rescale" | bc -l`
    param_angle_h_ow_h_k_r=`echo "$param_angle_h_ow_h_k * $pert_angle_rescale" | bc -l`
    param_angle_c_c_c_k_r=`echo "$param_angle_c_c_c_k * $pert_rescale" | bc -l`
    
    sed -e "/^param_angle_c_c_c_k/s/=.*/=$param_angle_c_c_c_k_r/g" $file |\
    sed -e "/^param_angle_h_ow_h_k/s/=.*/=$param_angle_h_ow_h_k_r/g" |\
    sed -e "/^param_bond_c_c_k/s/=.*/=$param_bond_c_c_k_r/g" |\
    sed -e "/^param_bond_ow_h_k/s/=.*/=$param_bond_ow_h_k_r/g" |\
    sed -e "/^param_dihedral_c_c_c_c_k/s/=.*/=$param_dihedral_c_c_c_c_k/g" > tmp.tmp.tmp

    # sed -e "/^param_nb_ch2_c06/s/=.*/=$param_nb_ch2_c06_r/g" $file |\
    # sed -e "/^param_nb_ch3_c06/s/=.*/=$param_nb_ch3_c06_r/g" |\
    # sed -e "/^param_nb_ow_c06/s/=.*/=$param_nb_ow_c06_r/g" |\
    # sed -e "/^param_nb_ch2_ow_c06/s/=.*/=$param_nb_ch2_ow_c06_r/g" |\
    # sed -e "/^param_nb_ch3_ow_c06/s/=.*/=$param_nb_ch3_ow_c06_r/g" |\
    # sed -e "/^param_nb_ch3_ch2_c06/s/=.*/=$param_nb_ch3_ch2_c06_r/g" |\
    # sed -e "/^param_pair_ow_ow_c06/s/=.*/=$param_pair_ow_ow_c06_r/g" |\
    # sed -e "/^param_pair_ch2_ch2_c06/s/=.*/=$param_pair_ch2_ch2_c06_r/g" |\
    # sed -e "/^param_pair_ch3_ch3_c06/s/=.*/=$param_pair_ch3_ch3_c06_r/g" |\
    # sed -e "/^param_pair_ch2_ow_c06/s/=.*/=$param_pair_ch2_ow_c06_r/g" |\
    # sed -e "/^param_pair_ch3_ow_c06/s/=.*/=$param_pair_ch3_ow_c06_r/g" |\
    # sed -e "/^param_pair_ch3_ch2_c06/s/=.*/=$param_pair_ch3_ch2_c06_r/g" |\
    # sed -e "/^param_charge_ch2/s/=.*/=$param_charge_ch2_r/g" |\
    # sed -e "/^param_charge_ch3a/s/=.*/=$param_charge_ch3a_r/g" |\
    # sed -e "/^param_charge_ch3b/s/=.*/=$param_charge_ch3b_r/g" |\
    # sed -e "/^param_charge_ow/s/=.*/=$param_charge_ow_r/g" |\
    # sed -e "/^param_charge_h/s/=.*/=$param_charge_h_r/g" |\
    # sed -e "/^param_angle_c_c_c_k/s/=.*/=$param_angle_c_c_c_k_r/g" $file |\

    mv -f tmp.tmp.tmp $file
}

function split_trr () {
    echo 2 0 | trjconv -f traj.trr -o out.gro -center
    equi_xvout_feq=`echo "$equi_frame_feq / $equi_dt" | bc -l | cut -d '.' -f 1`
    equi_warm_num_frame=`echo "$equi_warm_time / $equi_frame_feq" | bc -l | cut -d '.' -f 1`
    equi_nline_gro=`wc -l conf.gro | awk '{print $1}'`
    tmp_nline_del=`echo "$equi_nline_gro * ($equi_warm_num_frame+1)" | bc -l`
    sed "1,$tmp_nline_del d" out.gro > tmptmp.gro
    split -l $equi_nline_gro -a 5 -d tmptmp.gro
    mkdir -p equiConfs
    for i in `ls | grep "^x"`
    do
	tmp_number=`echo $i | cut -d 'x' -f 2 | cut -d 'c' -f 1`
	mv -f $i conf.$tmp_number.gro
	mv conf.$tmp_number.gro equiConfs
    done
    rm -f out.gro tmptmp.gro
}

function print_equi_frame_time () {
    rm -f equi.frame
    equi_frame_start_time=`echo "$equi_warm_time + $equi_frame_feq" | bc -l`
    equi_frame_end_time=`  echo "$equi_warm_time + $equi_frame_feq * $equi_num_frame" | bc -l`
    count=1
    for i in `seq $equi_frame_start_time $equi_frame_feq $equi_frame_end_time`
    do
	myprint=`echo "$i - $equi_frame_feq * 0.5" | bc -l`
	mycount=`printf %06d $count`
	echo "$mycount $myprint" >> equi.frame
	count=$(($count+1))
    done
}
