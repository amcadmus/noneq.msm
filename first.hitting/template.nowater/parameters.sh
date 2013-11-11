#!/bin/bash

fht_equi_dir=/home/mi/wanghan/study/noneq.msm/first.hitting/equi.systems/nowater.gromos.45a3.new/
fht_num_conf_use=200

fht_integrator=sd1
fht_dt=0.001			# ps
fht_stop_time=0.5		# ps
fht_frame_feq=0.005		# ps
fht_energy_feq=10.0		# ps

fht_noSdRange=0.0		# nm
fht_T=300			# K
fht_taut=0.1			# ps
fht_barostat=no			# Parrinello-Rahman or no
fht_taup=2.0			# ps

fht_meta_low=150		# deg.
fht_meta_up=-150		# deg.

# parameters for gaussian base
fht_gaussian_base_max=1.0
fht_gaussian_base_sigma=25
fht_gaussian_base_position="120 -120"

# parameters for cos base
fht_cos_base_number=6
fht_cos_base_k_file=cos.k.in

fht_parallel_num_pro=1		# n
fht_parallel_my_id=0		# n

grompp_command="grompp"
mdrun_command="mdrun -nt 1"

