#!/bin/bash

# non-equilibrium settings
pert_equi_result="$HOME/study/noneq.msm/accelerate/equi.systems/water.butane.e0.5.charge1-4.module.rf.tab/"
pert_num_conf_use=10	#
pert_mode=2			# 1: relax, 2: cos wave
pert_strength=1.0		# nm/ps velocity
pert_warm_time=40		# ps 1: warm time, 2: periodicity
pert_shift=0			# unitless
pert_phi=270			# deg. shift of the phase
pert_time=40			# ps
pert_frame_feq=0.5		# ps
pert_dt=0.0005			# ps
pert_rcut=1.3			# nm
pert_rsmooth=1.2		# nm
pert_erf=78			# 

pert_integrator=sd-baoab
pert_taut=5.0			# ps
pert_barostat=no		# Parrinello-Rahman or no
pert_taup=2.0			# ps

pert_rescale=5.0		# ps
pert_bond_rescale=1.0		# 
pert_angle_rescale=1.0		# 
pert_ele_rescale_start=0.00	# nm
pert_ele_rescale_end=0.2	# nm
pert_vdw_rescale_start=0.00	# nm
pert_vdw_rescale_end=0.2	# nm

pert_parallel_num_pro=1		# n
pert_parallel_my_id=0		# n


grompp_command="grompp"
mdrun_command="mdrun -v -tablep table.xvg"

pert_rescale_sqrt=`echo "sqrt($pert_rescale)" | bc -l`
