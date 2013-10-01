#!/bin/bash

dipole_name=dipole.name
angle_name=angle.name
out_dipole_name=consist.dipole.name
out_angle_name=consist.angle.name
bin_dir=/home/cocktail/wanghan/study/adress.noneq/alanine.nanma/template/tools/angles/
bin=$bin_dir/conf.dipole
make -C $bin_dir/ -j 8

echo "# firstly make the angle name and dipole name consistent"

rm -f $out_angle_name $out_dipole_name

for i in `grep -v \# $dipole_name`;
do
    pert_dir=`echo $i | sed 's/.*\(pert\.[0-9]*\).*/\1/g'`
#    echo $pert_dir
    grep $pert_dir $angle_name >> $out_angle_name
done
sort $out_angle_name > tmp.tmp
mv -f tmp.tmp $out_angle_name

for i in `grep -v \# $angle_name`;
do
    pert_dir=`echo $i | sed 's/.*\(pert\.[0-9]*\).*/\1/g'`
#    echo $pert_dir
    grep $pert_dir $dipole_name >> $out_dipole_name
done
sort $out_dipole_name > tmp.tmp
mv -f tmp.tmp $out_dipole_name


echo "# do the calculation"
$bin --input-angle $out_angle_name --input-dipole $out_dipole_name


