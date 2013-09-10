#!/bin/bash

begin=0
end=1000
step=1
rdfcut=1.2
target_dir=backup.result.perts

if test $# -ge 1; then
    begin=$1
fi
if test $# -ge 2; then
    step=$2
fi
if test $# -ge 3; then
    end=$3
fi
if test $# -ge 4; then
    target_dir=$4
fi

avg_bin="$HOME/study/adress.noneq/alanine.nanma/template/tools/average/average.file"
if test ! -f $avg_bin; then
    echo "# no file $avg_bin"
fi

echo "# cal rdf in time interval [ $begin : $step : $end ]"
if test ! -d $target_dir; then
    echo "# no dir $target_dir, exit"
    exit
fi

if test -d dir.avg.rdf; then
    echo "# using existing dir.avg.rdf"
#    echo "# existing dir.avg.rdf, remove"
#    rm -fr dir.avg.rdf
fi
mkdir -p dir.avg.rdf

echo "# make rdf.sol.name"
rm -f rdf.sol.name
for i in `ls $target_dir/`;
do
    echo "$target_dir/$i/rdf.xvg" >> rdf.sol.name
done
    

for ii in `seq $begin $step $end`
do
    echo "# doing time $ii"
    cwd=`pwd`
    for jj in `ls $target_dir/`
    do
	cd $target_dir/$jj
	rm -f rdf.xvg 
	echo 3 15 | g_rdf -n index.ndx -b $ii -e $ii -rdf mol_com -bin 0.01 &> /dev/null
	cat rdf.xvg | head -n 140 > tmp.xvg
	mv -f tmp.xvg rdf.xvg
	cd $cwd
    done
    pii=`printf "%.2f" $ii`
    piii=`echo $pii | cut -d '.' -f 1`
    piii=`printf "%05d" $piii`
    piid=`echo $pii | cut -d '.' -f 2`
    rm -f dir.avg.rdf/rdf.$piii.$piid.xvg
    $avg_bin -f rdf.sol.name -o dir.avg.rdf/rdf.$piii.$piid.xvg
done

