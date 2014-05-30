#!/bin/bash

fht_equi_dir=/home/mi/wanghan/study/noneq.msm/first.hitting/equi.systems/inwater.alanine.amber99sb-ildn/nvt.1000ns.confs
fht_equi_frame_traj=$fht_equi_dir/alphaR.trr
fht_equi_frame_count=$fht_equi_dir/alphaR.cnt
fht_coreset_data=$fht_equi_dir/coreset.beta.dat
fht_coreset_target=2
fht_coreset_notin=0
fht_num_conf_use=20000

fht_integrator=sd1
fht_dt=0.0005			# ps
fht_stop_time=5.0		# ps
fht_frame_feq=0.1		# ps
fht_energy_feq=10.0		# ps

fht_noSdRange=0.0		# nm
fht_T=300			# K
fht_taut=0.1			# ps
fht_barostat=no			# Parrinello-Rahman or no
fht_taup=2.0			# ps

fht_meta_low=150		# deg.
fht_meta_up=-150		# deg.

# parameters for bases 
fht_base_phi_number=6
fht_base_phi_k_file=base.k.phi
fht_base_psi_number=6
fht_base_psi_k_file=base.k.psi

fht_parallel_num_pro=1		# n
fht_parallel_my_id=0		# n

grompp_command="grompp"
mdrun_command="mdrun -nt 1"

