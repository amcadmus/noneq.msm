#!/bin/bash

target=`seq 0 0.2 30`
rm -f error.2d.gif.log

for i in $target;
do
    echo "# processing $i"
    printedf=`printf %.2f $i`
    int_i=`echo $printedf | cut -d '.' -f 1`
    printedi=`printf %05d $int_i`
    printednum=$printedi.`echo $printedf | cut -d '.' -f 2`
    out=fig-2d-noneq-$printednum.gif
    out_linear=fig-2d-linear-$printednum.gif
    file=indicator.vx.$printednum.out
    file_linear=indicator.linear.corr.vx.$printednum.out
    if test ! -f $file; then
	echo "# no file $file"
	continue
    fi
    if test ! -f $file_linear; then
	echo "# no file $file_linear"
	continue
    fi
    sed -e "/^set out/s/'.*'/'$out'/g" pl.2d.gif.one.gp |\
    sed -e "/^spl/s/'.*'/'$file'/g" > tmp.gp
    gnuplot tmp.gp 2>> error.2d.gif.log
    name_line_1="$name_line_1 $out"
    sed -e "/^set out/s/'.*'/'$out_linear'/g" pl.2d.gif.one.gp |\
    sed -e "/^spl/s/'.*'/'$file_linear'/g" > tmp.gp
    gnuplot tmp.gp 2>> error.2d.gif.log
    name_line_2="$name_line_2 $out_linear"

    sed -e "s/00000.00/$printednum/g" pl.2d.gif.tog.one.gp > tmp.gp
    gnuplot tmp.gp 2>> error.2d.gif.log
    name_line_3="$name_line_3 fig-2d-tog-$printednum.gif"
done

gifsicle --delay=10 $name_line_1 > fig-2d-noneq.gif
gifsicle --delay=10 $name_line_2 > fig-2d-linear.gif
gifsicle --delay=10 $name_line_3 > fig-2d-tog.gif
