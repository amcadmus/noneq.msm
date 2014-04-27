#!/bin/bash

equi_time=100
bin_size=15
threshold=1e-2
command=./ram.table
angle_phi=./angaver.phi.xvg
angle_psi=./angaver.psi.xvg
meta_dat=meta.set.dat

phi_l=`grep -v \# $meta_dat | head -n 1 | cut -f 1 -d ',' | cut -f 1 -d ':'`
phi_u=`grep -v \# $meta_dat | head -n 1 | cut -f 1 -d ',' | cut -f 2 -d ':'`
psi_l=`grep -v \# $meta_dat | head -n 1 | cut -f 2 -d ',' | cut -f 1 -d ':'`
psi_u=`grep -v \# $meta_dat | head -n 1 | cut -f 2 -d ',' | cut -f 2 -d ':'`

mycommand="$command --begin $equi_time --threshold $threshold --print-indicator 3 --phi-low $phi_l --phi-up $phi_u --psi-low $psi_l --psi-up $psi_u --input-angle-phi $angle_phi --input-angle-psi $angle_psi --bin-size-phi $bin_size --bin-size-psi $bin_size --output coreset.alphaL.dat"
echo "# excute $mycommand"
$mycommand

phi_l=`grep -v \# $meta_dat | tail -n 1 | cut -f 1 -d ',' | cut -f 1 -d ':'`
phi_u=`grep -v \# $meta_dat | tail -n 1 | cut -f 1 -d ',' | cut -f 2 -d ':'`
psi_l=`grep -v \# $meta_dat | tail -n 1 | cut -f 2 -d ',' | cut -f 1 -d ':'`
psi_u=`grep -v \# $meta_dat | tail -n 1 | cut -f 2 -d ',' | cut -f 2 -d ':'`

mycommand="$command --begin $equi_time --threshold $threshold --print-indicator 2 --phi-low $phi_l --phi-up $phi_u --psi-low $psi_l --psi-up $psi_u --input-angle-phi $angle_phi --input-angle-psi $angle_psi --bin-size-phi $bin_size --bin-size-psi $bin_size --output coreset.beta.dat"
echo "# excute $mycommand"
$mycommand

phi_l=`grep -v \# $meta_dat | head -n 2 | tail -n 1 | cut -f 1 -d ',' | cut -f 1 -d ':'`
phi_u=`grep -v \# $meta_dat | head -n 2 | tail -n 1 | cut -f 1 -d ',' | cut -f 2 -d ':'`
psi_l=`grep -v \# $meta_dat | head -n 2 | tail -n 1 | cut -f 2 -d ',' | cut -f 1 -d ':'`
psi_u=`grep -v \# $meta_dat | head -n 2 | tail -n 1 | cut -f 2 -d ',' | cut -f 2 -d ':'`

mycommand="$command --begin $equi_time --threshold $threshold --print-indicator 1 --phi-low $phi_l --phi-up $phi_u --psi-low $psi_l --psi-up $psi_u --input-angle-phi $angle_phi --input-angle-psi $angle_psi --bin-size-phi $bin_size --bin-size-psi $bin_size --output coreset.alphaR.dat"
echo "# excute $mycommand"
$mycommand


