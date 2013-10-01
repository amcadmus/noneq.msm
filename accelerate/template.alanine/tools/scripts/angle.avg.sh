#!/bin/bash

make -C ./tools/angles/ -j

cwd=`pwd`

for i in `cat angle.name`; 
do
    if test ! -f $i; then
	echo "# no file $i, make it"
	cd $(dirname $i)
	$cwd/tools/angles/evolve -f alanine.xtc -s angle.dat > /dev/null
	cd $cwd
    fi
done

./tools/angles/average.traj -f angle.name -o avg.angle --refh 5

test ! -d dir.avg.angle && mkdir dir.avg.angle

mv -f avg.angle* dir.avg.angle
