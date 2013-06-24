set term post eps enh color solid size 15cm,4cm
set out 'fig-tilt-str2-simple.eps'

#set border 4095 front linetype -1 linewidth 1.000
set view map
set samples 50, 50
set isosamples 50, 50
unset surface
set style data pm3d
set style function pm3d
set ticslevel 0
# set xlabel "q [nm]" 
# set ylabel "p [nm/ps]"
set xrange [-2.0:2.0]
set yrange [-8:8]
set xtics 1.0
set ytics 4.0
set mxtics 2
set mytics 2
set pm3d implicit at b
#set size ratio 0.88
# set palette rgbformulae 22,13,-31
set palette rgbformulae 21,22,23
set size ratio 1.0

# set yrange [-8:8]
# set cbrange [-0.02:0.08]
set cbrange [-0.06:0.06]
set cbtics 0.02
set format cb "%.2f"

#set size 1.4,1.0
set origin 0, 0.0

set multiplot 

# set ratio 1
set size 0.24,1

set origin 0.0,0.0
# set xlabel 'q [nm]'
set xlabel '(a)'
set title 't=0 ps, EX'
# set title 'linear resp., t=00000 ps' font "tt,16"
# unset label
# set label "(a)" at -2.8, 11
# spl 'tilt.warm020.str2.0.nst1e09.smallgrid/indicator.vx.00000.00.out' u 1:2:(-$3) not
spl 'tilt.warm020.str2.0.nst1e10.noneq/indicator.vx.00000.00.out' u 1:2:(-$3) not

unset ylabel

set origin 0.25, 0.0
set xlabel '(b)'
set title 't=20 ps, EX'
# unset label
# set label "(b)" at -2.8, 11
# spl 'tilt.warm020.str2.0.nst1e09.smallgrid/indicator.vx.00020.00.out' u 1:2:(-$3) not
spl 'tilt.warm020.str2.0.nst1e10.noneq/indicator.vx.00020.00.out' u 1:2:(-$3) not

set origin 0.5, 0.0
set xlabel '(c)'
set title 't=20 ps, resp. n=1'
# unset label
# set label "(c)" at -2.8, 11
spl 'tilt.warm020.str2.0.nst1e09.ref0.0.resp.order1/indicator.resp.vx.00020.00.out' u 1:2:(-$3) not

set origin 0.75, 0.0
set xlabel '(d)'
set title '  t = 20 ps, resp. n=2'
# unset label
# set label "(d)" at -2.8, 11
spl 'tilt.warm020.str2.0.nst1e09.ref0.0.resp.order2/indicator.resp.vx.00020.00.out' u 1:2:(-$3) not


unset multiplot
