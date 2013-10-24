#!/bin/bash

# non-equilibrium settings
pert_equi_result="$HOME/study/noneq.msm/accelerate/equi.systems/alanine.charmm.pme.tab/"
pert_num_conf_use=100000	#
pert_mode=2			# 1: relax, 2: cos wave
pert_strength=1.0		# nm/ps velocity
pert_warm_time=40		# ps 1: warm time, 2: periodicity
pert_shift=0			# unitless
pert_phi=270			# deg. shift of the phase
pert_time=100			# ps
pert_frame_feq=0.5		# ps
pert_dt=0.0005			# ps
pert_rlist=1.33			# nm
pert_rcut=1.0			# nm
pert_rsmooth=0.9		# nm
pert_erf=78			# 
pert_fourierspacing=0.12	# nm
pert_pme_order=4		#

pert_integrator=sd-baoab
pert_taut=0.2			# ps
pert_noSdRange=1.0		# nm
pert_barostat=no		# Parrinello-Rahman or no
pert_taup=2.0			# ps

pert_rescale=1.0		# ps
pert_bond_T_cap=8.631997e-15	# s taken from water
pert_angle_T_cap=6.955613e-15	# s taken from water
pert_ele_rescale_scale0=1.0	# 
pert_vdw_rescale_scale0=1.0	# 

pert_ele_rescale_start=0.00	# nm
pert_ele_rescale_end=0.01	# nm
pert_vdw_rescale_start=0.00	# nm
pert_vdw_rescale_end=0.01	# nm

pert_parallel_num_pro=1		# n
pert_parallel_my_id=0		# n


grompp_command="grompp"
mdrun_command="mdrun -v"

pert_rescale_sqrt=`echo "sqrt($pert_rescale)" | bc -l`
