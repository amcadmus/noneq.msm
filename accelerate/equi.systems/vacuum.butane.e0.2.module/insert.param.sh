#!/bin/bash

function insert_param () {
    file=$1
    outfile=$2
    sed -e "s/param_nb_ch2_c06/$param_nb_ch2_c06/g" $file |\
    sed -e "s/param_nb_ch2_c12/$param_nb_ch2_c12/g" |\
    sed -e "s/param_nb_ch3_c06/$param_nb_ch3_c06/g" |\
    sed -e "s/param_nb_ch3_c12/$param_nb_ch3_c12/g" |\
    sed -e "s/param_nb_ow_c06/$param_nb_ow_c06/g" |\
    sed -e "s/param_nb_ow_c12/$param_nb_ow_c12/g" |\
    sed -e "s/param_nb_ch2_ow_c06/$param_nb_ch2_ow_c06/g" |\
    sed -e "s/param_nb_ch2_ow_c12/$param_nb_ch2_ow_c12/g" |\
    sed -e "s/param_nb_ch3_ow_c06/$param_nb_ch3_ow_c06/g" |\
    sed -e "s/param_nb_ch3_ow_c12/$param_nb_ch3_ow_c12/g" |\
    sed -e "s/param_nb_ch3_ch2_c06/$param_nb_ch3_ch2_c06/g" |\
    sed -e "s/param_nb_ch3_ch2_c12/$param_nb_ch3_ch2_c12/g" |\
    sed -e "s/param_pair_ow_ow_c06/$param_pair_ow_ow_c06/g" |\
    sed -e "s/param_pair_ow_ow_c12/$param_pair_ow_ow_c12/g" |\
    sed -e "s/param_pair_ch2_ch2_c06/$param_pair_ch2_ch2_c06/g" |\
    sed -e "s/param_pair_ch2_ch2_c12/$param_pair_ch2_ch2_c12/g" |\
    sed -e "s/param_pair_ch3_ch3_c06/$param_pair_ch3_ch3_c06/g" |\
    sed -e "s/param_pair_ch3_ch3_c12/$param_pair_ch3_ch3_c12/g" |\
    sed -e "s/param_pair_ch2_ow_c06/$param_pair_ch2_ow_c06/g" |\
    sed -e "s/param_pair_ch2_ow_c12/$param_pair_ch2_ow_c12/g" |\
    sed -e "s/param_pair_ch3_ow_c06/$param_pair_ch3_ow_c06/g" |\
    sed -e "s/param_pair_ch3_ow_c12/$param_pair_ch3_ow_c12/g" |\
    sed -e "s/param_pair_ch3_ch2_c06/$param_pair_ch3_ch2_c06/g" |\
    sed -e "s/param_pair_ch3_ch2_c12/$param_pair_ch3_ch2_c12/g" |\
    sed -e "s/param_charge_ch2/$param_charge_ch2/g" |\
    sed -e "s/param_charge_ch3a/$param_charge_ch3a/g" |\
    sed -e "s/param_charge_ch3b/$param_charge_ch3b/g" |\
    sed -e "s/param_charge_ow/$param_charge_ow/g" |\
    sed -e "s/param_charge_h/$param_charge_h/g" |\
    sed -e "s/param_bond_c_c_k/$param_bond_c_c_k/g" |\
    sed -e "s/param_angle_c_c_c_k/$param_angle_c_c_c_k/g" |\
    sed -e "s/param_dihedral_c_c_c_c_k/$param_dihedral_c_c_c_c_k/g" > $outfile
}