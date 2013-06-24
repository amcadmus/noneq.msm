set term post eps enh color font 16 size 8.5cm,6cm
set out 'fig-ctr-stat.eps'

set style line 01 lt 2 lc 0 lw 3 pt 3
set style line 02 lt 1 lc 0 lw 3 pt 3
set style line 10 lt 1 lc 1 lw 3 pt 3
set style line 11 lt 2 lc 1 lw 3 pt 3
set style line 20 lt 1 lc 2 lw 3 pt 3
set style line 21 lt 2 lc 2 lw 3 pt 3
set style line 30 lt 1 lc 3 lw 3 pt 3
set style line 4 lt 1 lc 4 lw 3 pt 3

set xlabel "t [ps]"

set lmargin 8
set rmargin 8

set multiplot
set grid

set origin 0,0
set size 1, 0.58
set xrange [0:20]
set yrange [0.5:0.8]
set xtics 5
set mxtics 5
set ytics 0.1
set mytics 5
set ylabel "O[F_e](t)"
set format y "%.2f"
set key left bottom

pl -0.5 ls 01 not, \
'ctrl.beta0.01.res1.00.1/state.step030.out' u 1:(-$2) w l ls 11 not '{/Symbol b} = 0.01',\
'ctrl.beta0.01.res1.00.1/state.step030.out' u 1:(-$4) w l ls 10 not '{/Symbol b} = 0.01',\
'ctrl.beta0.05.res1.00.1/state.step030.out' u 1:(-$2) w l ls 21 not '{/Symbol b} = 0.05',\
'ctrl.beta0.05.res1.00.1/state.step030.out' u 1:(-$4) w l ls 20 not '{/Symbol b} = 0.05'

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

set style arrow 1 head filled size screen 0.015,20,20 ls 02
set arrow from 3,1 to 3,0.2 as 1
set arrow from 12,1.4 to 12,0.7 as 1
set arrow from 16.5, 3.2 to 18.8, 3.7 as 1
set arrow from 22, 3.2 to 20.2, 2.5 as 1

pl 0.0 ls 01 not, \
'ctrl.beta0.01.res1.00.1/ctr.step030.out' u 1:2 w l ls 10 not '{/Symbol b} = 0.01',\
'ctrl.beta0.05.res1.00.1/ctr.step030.out' u 1:2 w l ls 20 not '{/Symbol b} = 0.05'

set lmargin 0
set rmargin 0
set tmargin 0
set bmargin 0
set size 0.13 , 0.2
set format x  ""
set format y  ""
unset xlabel
unset ylabel
unset grid
set xtics 2
set ytics 10
set mxtics 1
set mytics 1
set noxtic
set noytic
set noborder
unset arrow 

set origin 0.2, 0.55
pl [-2:2][-10:7] 0.5 * 8 * (x*x - 1*1)**2 + -x * 0 w l ls 30 not

set origin 0.5, 0.58
pl [-2:2][-10:7] 0.5 * 8 * (x*x - 1*1)**2 + -x * 0.54 w l ls 30 not

set origin 0.6, 0.75
pl [-2:2][-10:7] 0.5 * 8 * (x*x - 1*1)**2 + -x * 3.73 w l ls 30 not

set origin 0.85, 0.75
pl [-2:2][-10:7] 0.5 * 8 * (x*x - 1*1)**2 + -x * 2.52 w l ls 30 not

unset multiplot