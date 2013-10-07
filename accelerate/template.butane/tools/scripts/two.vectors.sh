#!/bin/bash

bin_dir=/home/cocktail/wanghan/study/adress.noneq/alanine.nanma/template/tools/angles/
bin=$bin_dir/two.vectors
make -C $bin_dir/ -j 8


echo "# make xtc name file from angle name file"

rm -f xtc.name

for i in `cat angle.name`;
do
    echo $i | sed 's/angle.dat/alanine\.xtc/g' >> xtc.name
#    echo $i | sed 's/\(.*\)\/angle.dat/\1\//g' 
done

$bin --input-name xtc.name

