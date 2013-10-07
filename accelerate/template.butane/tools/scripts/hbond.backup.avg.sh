#!/bin/bash

avg_bin="$HOME/study/adress.noneq/alanine.nanma/template/tools/average/average.file"
if test ! -f $avg_bin; then
    echo "# no file $avg_bin"
fi

target_dir=backup.result.perts
if test $# -ge 1; then
    target_dir=$1
fi
angle_cut=20
dist_cut=0.3

rm -f hbond.o.name
rm -f hbond.n.name
cwd=`pwd`

source parameters.sh

for ii in `ls $target_dir | grep pert`
do
    echo "# doing $target_dir/$ii"
    cd $target_dir/$ii
    echo "a O" > command.tmp
    echo "a N and HN" >> command.tmp
    echo "q" >> command.tmp
    cat command.tmp | make_ndx -f confout.gro -o tmp.ndx &> /dev/null
    rm -f command.tmp
    echo 13 15 | g_hbond -n tmp.ndx -nonitacc -a $angle_cut -r $dist_cut -nice 0 &> /dev/null
    grep -v '@' hbnum.xvg > hbnum.o.xvg
    rm -f hbnum.xvg
    echo 13 16 | g_hbond -n tmp.ndx -nonitacc -a $angle_cut -r $dist_cut -nice 0 &> /dev/null
    grep -v '@' hbnum.xvg > hbnum.n.xvg
    rm -f hbnum.xvg
    echo 13 13 | g_hbond -n tmp.ndx -nonitacc -a $angle_cut -r $dist_cut -nice 0 &> /dev/null
    grep -v '@' hbnum.xvg > hbnum.sol.xvg
    rm -f hbnum.xvg
    echo $target_dir/$ii/hbnum.o.xvg >> $cwd/hbond.o.name
    echo $target_dir/$ii/hbnum.n.xvg >> $cwd/hbond.n.name
    echo $target_dir/$ii/hbnum.sol.xvg >> $cwd/hbond.sol.name
    cd $cwd
done

$avg_bin -f hbond.o.name -o avg.hbond.o.out
$avg_bin -f hbond.n.name -o avg.hbond.n.out
$avg_bin -f hbond.sol.name -o avg.hbond.sol.out
