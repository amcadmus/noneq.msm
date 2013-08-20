set term post eps enh color font 16 size 8.5cm,6cm
set out 'fig-ctr-stat-2.eps'

set style line 01 lt 2 lc 0 lw 3 pt 3
set style line 02 lt 1 lc 0 lw 3 pt 3
set style line 10 lt 1 lc 1 lw 3 pt 3
set style line 11 lt 2 lc 1 lw 3 pt 3
set style line 12 lt 4 lc 1 lw 3 pt 3
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
set xrange [0:9]
set yrange [0.5:0.8]
set xtics 1
set mxtics 2
set ytics 0.1
set mytics 5
# set ylabel "{/Symbol h}I and P(q_t)"
set ylabel "I and P"
set format y "%.2f"
set key at 4.5, 0.77 

pl -0.5 ls 01 not, \
'ctrl.beta0.01.res1.00.2.ref.fdiff0.2.more/state.step011.out' u 1:(-$2) w l ls 21 not 'ref {/Symbol b} = 0.01',\
'ctrl.beta0.01.res1.00.2.ref.fdiff0.2.more/state.step011.out' u 1:(-$4) w l ls 20 t 'finite diff.',\
'ctrl.beta0.01.res1.00.2.more/state.step022.out' u 1:(-$2) w l ls 11 not '{/Symbol b} = 0.01',\
'ctrl.beta0.01.res1.00.2.more/state.step022.out' u 1:(-$4) w l ls 10 t 'linear resp.'

# 'ctrl.beta0.05.res1.00.2/state.step012.out' u 1:(-$2) w l ls 21 not '{/Symbol b} = 0.05',\
# 'ctrl.beta0.05.res1.00.2/state.step012.out' u 1:(-$4) w l ls 20 not '{/Symbol b} = 0.05'

set origin 0,0.52
set size 1, 0.5
set xrange [0:9]
set yrange [0:5]
set format x  ""
set ytics 1
set mytics 2
unset xlabel
set ylabel "u(t)"
set key at -0.5,3.5

set style arrow 1 head filled size screen 0.015,20,20 ls 02
set arrow from 2,1.5 to 2,0.7 as 1
set arrow from 5.5,2.5 to 5.95,1.6 as 1
set arrow from 7.2, 4.2 to 7.9, 4.2 as 1
set arrow from 9.8, 4.0 to 9.1, 3.1 as 1

pl 0.0 ls 01 not, \
'ctrl.beta0.01.res1.00.2.ref.fdiff0.2.more/ctr.step011.out' u 1:2 w l ls 20 not 'finite diff.',\
'ctrl.beta0.01.res1.00.2.more/ctr.step022.out' u 1:2 w l ls 10 not 'linear resp.'
# 'ctrl.beta0.05.res1.00.2/ctr.step012.out' u 1:2 w l ls 20 not '{/Symbol b} = 0.05'

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

set origin 0.25, 0.57
pl [-2:2][-10:7] 0.5 * 8 * (x*x - 1*1)**2 + -x * 0.493 w l ls 30 not

set origin 0.45, 0.63
pl [-2:2][-10:7] 0.5 * 8 * (x*x - 1*1)**2 + -x * 1.489 w l ls 30 not

set origin 0.58, 0.75
pl [-2:2][-10:7] 0.5 * 8 * (x*x - 1*1)**2 + -x * 4.176 w l ls 30 not

set origin 0.85, 0.75
pl [-2:2][-10:7] 0.5 * 8 * (x*x - 1*1)**2 + -x * 2.982 w l ls 30 not

unset multiplot