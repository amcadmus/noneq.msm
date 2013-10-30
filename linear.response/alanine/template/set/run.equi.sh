#!/bin/bash

source parameters.sh
source functions.sh
rm -f run.log

if test -d result.equi; then
    echo "existing dir result.equi, mv and backup"
    mv result.equi result.equi.`date +%s`
fi
if test ! -d $gro_dir; then
    echo "no dir $gro_dir, exit"
    exit
fi
cp -a $gro_dir result.equi

cd result.equi
set_parameters_equi grompp.mdp

echo "# run with command `which grompp`"
$grompp_command &>> run.log
echo "# run with command `which mdrun`"
$mdrun_command &>> run.log

split_trr
cd ..
