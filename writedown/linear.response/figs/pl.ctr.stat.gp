set term post eps enh color solid font 16 size 8.5cm,6cm
set out 'fig-ctr-stat.eps'

set style line 01 lt 0 lc 0 lw 3 pt 3
set style line 10 lt 1 lc 1 lw 3 pt 3
set style line 11 lt 0 lc 1 lw 3 pt 3
set style line 20 lt 1 lc 2 lw 3 pt 3
set style line 21 lt 0 lc 2 lw 3 pt 3
set style line 3 lt 1 lc 3 lw 3 pt 3
set style line 4 lt 1 lc 4 lw 3 pt 3

set xlabel "t [ps]"

set lmargin 8

set multiplot

set origin 0,0
set size 1, 0.58
set xrange [0:20]
set yrange [-0.7:-0.5]
set xtics 5
set mxtics 5
set ytics 0.05
set mytics 5
set ylabel "O[F_e](t)"
set format y "%.2f"
set key left bottom

pl -0.5 ls 01 not, \
'ctrl.beta0.01.res1.00/state.step030.out' u 1:2 w l ls 11 not '{/Symbol b} = 0.01',\
'ctrl.beta0.01.res1.00/state.step030.out' u 1:4 w l ls 10 not '{/Symbol b} = 0.01',\
'ctrl.beta0.05.res1.00/state.step030.out' u 1:2 w l ls 21 not '{/Symbol b} = 0.05',\
'ctrl.beta0.05.res1.00/state.step030.out' u 1:4 w l ls 20 not '{/Symbol b} = 0.05'

set origin 0,0.52
set size 1, 0.5
set xrange [0:20]
set yrange [0:4]
set xtics 5
set mxtics 5
set format x  ""
set ytics 1
set mytics 2
unset xlabel
set ylabel "F_e(t)"
set key left top

pl 0.0 ls 01 not, \
'ctrl.beta0.01.res1.00/ctr.step030.out' u 1:2 w l ls 10 not '{/Symbol b} = 0.01',\
'ctrl.beta0.05.res1.00/ctr.step030.out' u 1:2 w l ls 20 not '{/Symbol b} = 0.05'
unset multiplot