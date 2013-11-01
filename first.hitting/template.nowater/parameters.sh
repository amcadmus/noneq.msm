#!/bin/bash

fht_equi_dir=
fht_num_conf_use=100000

fht_integrator=sd1
fht_dt=0.001			# ps
fht_time=1000			# ps
fht_frame_feq=0.1		# ps
fht_energy_feq=0.1		# ps

fht_noSdRange=0.0		# nm
fht_taut=1.0			# ps
fht_barostat=no			# Parrinello-Rahman or no
fht_taup=2.0			# ps

# parameters for gaussian base
fht_gaussian_base_max=10
fht_gaussian_base_sigma=25
fht_gaussian_base_position=


fht_parallel_num_pro=1		# n
fht_parallel_my_id=0		# n

grompp_command="grompp"
mdrun_command="mdrun -nt 1"

