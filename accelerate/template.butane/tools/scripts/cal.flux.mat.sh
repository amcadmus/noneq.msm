#!/bin/bash

file="meta.flux.inte.out"
up=1000
lo=600
if test $# -ge 1; then
    lo=$1
fi
if test $# -ge 2; then
    up=$2
fi
if test $# -ge 3; then
    file=$3
fi

echo "# doing linear fitting for $file in [$lo:$up]"

echo "set xrange [$lo:$up]"		>  tmp.gp
echo "f(x) = a*x + b"			>> tmp.gp
echo "fit f(x) '$file' u 1:3 via a,b"	>> tmp.gp

for ii in `seq 0 4`;
do
    line=""
    for jj in `seq 0 4`;
    do
#	echo $ii $jj
	if test $ii -eq $jj; then
	    line="$line   0.000"
	else
	    index=`echo "$ii * 5 + $jj + 2" | bc`
	    sed "s/u 1:[0-9]* /u 1:$index /g" tmp.gp > tmp.tmp
	    mv -f tmp.tmp tmp.gp
	    gnuplot tmp.gp &> output.tmp
	    a=`grep a output.tmp | grep = | tail -n 1 | awk '{print $3}'`
	    ta=`echo "$a * 1000" | octave | grep = | awk '{print $3}'`
	    pa=`printf %.3f $ta`
#	    echo $a $ta $pa
#	    echo $line  $pa
	    line="$line   $pa"
	fi
    done
    echo $line
done


