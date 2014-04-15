#!/bin/bash

fht_equi_dir=/home/mi/wanghan/study/noneq.msm/committor/equi.system/nowater.gromos.45a3.confs
fht_equi_frame_name=equi.frame.080-085
fht_coreset_data=`pwd`/coreset.dat
fht_num_conf_use=5000

fht_integrator=sd1
fht_dt=0.0005			# ps
fht_stop_time=1.0		# ps
fht_frame_feq=0.005		# ps
fht_energy_feq=10.0		# ps

fht_noSdRange=0.0		# nm
fht_T=300			# K
fht_taut=0.1			# ps
fht_barostat=no			# Parrinello-Rahman or no
fht_taup=2.0			# ps

fht_meta_low=150		# deg.
fht_meta_up=-150		# deg.

# parameters for cos base
fht_cos_base_number=6
fht_cos_base_k_file=cos.k.in

fht_parallel_num_pro=1		# n
fht_parallel_my_id=0		# n

grompp_command="grompp"
mdrun_command="mdrun -nt 1"

