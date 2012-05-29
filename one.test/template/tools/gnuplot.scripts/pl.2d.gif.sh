#!/bin/bash

target=`seq 0 1 150`
rm -f error.2d.gif.log

for i in $target;
do
    echo "# processing $i"
    printedi=`printf %05d $i`
    out=fig-2d-noneq-$printedi.gif
    file=indicator.vx.$printedi.00.out
    if test ! -f $file; then
	echo "# no file $file"
	continue
    fi
    sed -e "/^set out/s/'.*'/'$out'/g" pl.2d.gif.one.gp |\
    sed -e "/^spl/s/'.*'/'$file'/g" > tmp.gp
    gnuplot tmp.gp 2>> error.2d.gif.log
    name_line_1="$name_line_1 $out"
    out=fig-2d-linear-$printedi.gif
    file=indicator.linear.vx.$printedi.00.out
    sed -e "/^set out/s/'.*'/'$out'/g" pl.2d.gif.one.gp |\
    sed -e "/^spl/s/'.*'/'$file'/g" > tmp.gp
    gnuplot tmp.gp 2>> error.2d.gif.log
    name_line_2="$name_line_2 $out"

    sed -e "s/00000/$printedi/g" pl.2d.gif.tog.one.gp > tmp.gp
    gnuplot tmp.gp 2>> error.2d.gif.log
    name_line_3="$name_line_3 fig-2d-tog-$printedi.gif"
done

gifsicle --delay=5 $name_line_1 > fig-2d-noneq.gif
gifsicle --delay=5 $name_line_2 > fig-2d-linear.gif
gifsicle --delay=5 $name_line_3 > fig-2d-tog.gif
