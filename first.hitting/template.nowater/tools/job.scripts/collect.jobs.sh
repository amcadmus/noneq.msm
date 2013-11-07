#!/bin/bash

num_of_proc=16
seed_dir=set
stop_num=$(($num_of_proc-1))

cwd=`pwd`

rm -f $seed_dir/angle.name
rm -f $seed_dir/gxs.name
rm -f $seed_dir/energy.name
rm -f $seed_dir/success.dir.name

for i in `seq 0 $stop_num`
do
    printi=`printf %03d $i`

    echo "# collect $seed_dir.$printi"
    cd $seed_dir
    if test ! -f base.info; then
	echo "# copy base.info from $seed_dir.$printi"
	cp ../$seed_dir.$printi/base.info .
    fi

    test -f ../$seed_dir.$printi/angle.name  && cat ../$seed_dir.$printi/angle.name  >> ./angle.name
    test -f ../$seed_dir.$printi/gxs.name    && cat ../$seed_dir.$printi/gxs.name    >> ./gxs.name
    test -f ../$seed_dir.$printi/energy.name && cat ../$seed_dir.$printi/energy.name >> ./energy.name
    test -f ../$seed_dir.$printi/success.dir.name && cat ../$seed_dir.$printi/success.dir.name >> ./success.dir.name
    
    test ! -d result.fhts && mkdir result.fhts
    cd result.fhts
    if test -d ../../$seed_dir.$printi/result.fhts; then
	ln -sf ../../$seed_dir.$printi/result.fhts/fht.*   .
    fi
    cd ..
    
    cd ..
done


for i in `seq 0 $stop_num`
do
    printi=`printf %03d $i`

    echo "# collect $seed_dir.$printi"
    cd $seed_dir
    test ! -d backup.result.fhts && mkdir backup.result.fhts
    cd backup.result.fhts
    if test -d ../../$seed_dir.$printi/result.fhts; then
	ln -sf ../../$seed_dir.$printi/result.fhts/backup.fht.*   .
    fi
    cd ..    
    cd ..
done
    