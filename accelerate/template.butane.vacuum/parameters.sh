#!/bin/bash

# non-equilibrium settings
pert_equi_result="$HOME/study/adress.noneq/alanine.nanma/results/run.result.equi.box02.70/result.equi/"
pert_num_conf_use=100000	#
pert_mode=2			# 1: relax, 2: cos wave
pert_strength=1.0		# nm/ps velocity
pert_warm_time=10		# ps 1: warm time, 2: periodicity
pert_shift=1			# unitless
pert_phi=180			# deg. shift of the phase
pert_time=20			# ps
pert_frame_feq=0.5		# ps
pert_dt=0.002			# ps
pert_taut=0.1			# ps
pert_noSdRange=1.0		# nm
pert_barostat=Parrinello-Rahman # Parrinello-Rahman or no
pert_taup=2.0			# ps
pert_parallel_num_pro=1		# n
pert_parallel_my_id=0		# n

grompp_command="grompp -n index.ndx"
mdrun_command="mdrun"

