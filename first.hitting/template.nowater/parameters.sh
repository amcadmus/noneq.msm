#!/bin/bash

fht_equi_dir=/home/mi/wanghan/study/noneq.msm/first.hitting/equi.systems/nowater.gromos.45a3.new/
fht_num_conf_use=20

fht_integrator=sd1
fht_dt=0.001			# ps
fht_time=1000			# ps
fht_frame_feq=0.02		# ps
fht_energy_feq=10.0		# ps

fht_noSdRange=0.0		# nm
fht_taut=1.0			# ps
fht_barostat=no			# Parrinello-Rahman or no
fht_taup=2.0			# ps

fht_meta_low=150		# deg.
fht_meta_up=-150		# deg.

# parameters for gaussian base
fht_gaussian_base_max=10
fht_gaussian_base_sigma=25
fht_gaussian_base_position=

fht_parallel_num_pro=1		# n
fht_parallel_my_id=0		# n

grompp_command="grompp"
mdrun_command="mdrun -nt 1"

