#!/bin/bash

gmx_dt=0.001
gmx_time=2000000
gmx_frame_feq=1.0
gmx_log_feq=100.0

equi_time=100.0
cis_start=-100.0
cis_end=100.0
tool_dir=./tools

grompp_command="grompp"
mdrun_command="mdrun"

