#!/bin/bash

source env.sh
source parameters.sh

make -C $tool_dir/make.equi -j8

function set_parameters_init () {
    file=$1
    init_nstep=`echo "$gmx_time / $gmx_dt" | bc -l | cut -f 1 -d '.'`
    init_xtc_feq=`echo "$gmx_frame_feq / $gmx_dt" | bc -l | cut -f 1 -d '.'`
    init_trr_feq=`echo "$gmx_frame_feq / $gmx_dt" | bc -l | cut -f 1 -d '.'`
    init_log_feq=`echo "$gmx_log_feq / $gmx_dt" | bc -l | cut -f 1 -d '.'`
    init_energy_feq=`echo "$gmx_log_feq / $gmx_dt" | bc -l | cut -f 1 -d '.'`
#    echo "init_xtc_feq is $init_xtc_feq"
    init_seed=`date +%s`
    sed -e "/^dt/s/=.*/= $gmx_dt/g" $file |\
    sed -e "/^nstep/s/=.*/= $init_nstep/g" |\
    sed -e "/^ld-seed/s/=.*/= $init_seed/g" |\
    sed -e "/^gen_vel/s/=.*/= no/g" |\
    sed -e "/^nstlog/s/=.*/= $init_log_feq/g" |\
    sed -e "/^nstenergy/s/=.*/= $init_energy_feq/g" |\
    sed -e "/^nstxout/s/=.*/= $init_trr_feq/g" |\
    sed -e "/^nstvout/s/=.*/= $init_trr_feq/g" |\
    sed -e "/^nstxtcout/s/=.*/= $init_xtc_feq/g" > tmptmptmp.mdp
    mv -f tmptmptmp.mdp $file
}

rm -f run.log

set_parameters_init grompp.mdp
grompp &> run.log
mdrun &> run.log

g_angle -type dihedral -od angdist.xvg  -ov angaver.xvg -xvg none -nice 0

./$tool_dir/make.equi/main -f angaver.xvg --lower $cis_start --upper $cis_end --start $equi_time
