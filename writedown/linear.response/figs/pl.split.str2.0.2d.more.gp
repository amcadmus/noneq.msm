set term post eps enh color solid size 17cm,14cm
set out 'fig-split-str2-2d-more.eps'

#set border 4095 front linetype -1 linewidth 1.000
set view map
set samples 50, 50
set isosamples 50, 50
unset surface
set style data pm3d
set style function pm3d
set ticslevel 0
set xlabel "q [nm]" 
set ylabel "p [nm/ps]"
set xtics 1.0
set ytics 2.0
set pm3d implicit at b
#set size ratio 0.88
set palette rgbformulae 22,13,-31
set size ratio 1.0

# set yrange [-8:8]
set xrange [-1.5:1.5]
set mxtics 2
set yrange [-6:6]

set cbrange [-0.04:0.1]
set format cb "%.2f"

#set size 1.4,1.0
set origin 0, -0.015

set multiplot 

# set ratio 1
set size 0.24,0.3

set origin 0.0,0.69
unset xlabel
set title 't = 0 ps'
spl 'split.warm020.str2.0.nst1e09.smallgrid/indicator.vx.00000.00.out' u 1:2:(-$3) not

unset title
set origin 0.0, 0.46
spl 'split.warm020.str2.0.nst1e09.smallgrid/indicator.linear.vx.00000.00.out' u 1:2:(-$3) not

set origin 0.0, 0.23
spl 'split.warm020.str2.0.nst1e09.ref0.0.resp.order2/indicator.resp.vx.00000.00.out' u 1:2:(-$3) not

set origin 0.0, 0.0
set xlabel "q [nm]" 
spl 'split.warm020.str2.0.nst1e09.ref1.0.resp.order1/indicator.resp.vx.00000.00.out' u 1:2:(-$3) not


unset ylabel


set origin 0.25, 0.69
unset xlabel
set title 't = 5 ps'
spl 'split.warm020.str2.0.nst1e09.smallgrid/indicator.vx.00005.00.out' u 1:2:(-$3) not

set origin 0.5, 0.69
set title 't = 10 ps'
spl 'split.warm020.str2.0.nst1e09.smallgrid/indicator.vx.00010.00.out' u 1:2:(-$3) not

set origin 0.75, 0.69
set title 't = 20 ps'
spl 'split.warm020.str2.0.nst1e09.smallgrid/indicator.vx.00020.00.out' u 1:2:(-$3) not


unset title

set origin 0.25, 0.46
spl 'split.warm020.str2.0.nst1e09.smallgrid/indicator.linear.corr.vx.00005.00.out' u 1:2:(-$3) not

set origin 0.5, 0.46
spl 'split.warm020.str2.0.nst1e09.smallgrid/indicator.linear.corr.vx.00010.00.out' u 1:2:(-$3) not

set origin 0.75, 0.46
spl 'split.warm020.str2.0.nst1e09.smallgrid/indicator.linear.corr.vx.00020.00.out' u 1:2:(-$3) not


set origin 0.25, 0.23
spl 'split.warm020.str2.0.nst1e09.ref0.0.resp.order2/indicator.resp.vx.00005.00.out' u 1:2:(-$3) not

set origin 0.5, 0.23
spl 'split.warm020.str2.0.nst1e09.ref0.0.resp.order2/indicator.resp.vx.00010.00.out' u 1:2:(-$3) not

set origin 0.75, 0.23
spl 'split.warm020.str2.0.nst1e09.ref0.0.resp.order2/indicator.resp.vx.00020.00.out' u 1:2:(-$3) not


set xlabel "q [nm]"

set origin 0.25, 0.0
spl 'split.warm020.str2.0.nst1e09.ref1.0.resp.order1/indicator.resp.vx.00005.00.out' u 1:2:(-$3) not

set origin 0.5, 0.0
spl 'split.warm020.str2.0.nst1e09.ref1.0.resp.order1/indicator.resp.vx.00010.00.out' u 1:2:(-$3) not

set origin 0.75, 0.0
spl 'split.warm020.str2.0.nst1e09.ref1.0.resp.order1/indicator.resp.vx.00020.00.out' u 1:2:(-$3) not


unset multiplot
