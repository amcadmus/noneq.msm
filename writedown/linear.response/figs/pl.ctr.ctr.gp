set term post eps enh color solid font 16 size 8.5cm,6cm
set out 'fig-ctr-ctr.eps'

set style line 10 lt 0 lc 0 lw 3 pt 3
set style line 1 lt 1 lc 1 lw 3 pt 3
set style line 2 lt 1 lc 2 lw 3 pt 3
set style line 3 lt 1 lc 3 lw 3 pt 3
set style line 4 lt 1 lc 4 lw 3 pt 3

set xlabel "t [ps]"
set ylabel "F_e(t)"
set key left top

pl 0 ls 0 not, \
'ctrl.beta0.01.res1.00/ctr.step030.out' w lp ls 1 t '{/Symbol b} = 0.01',\
'ctrl.beta0.05.res1.00/ctr.step030.out' w lp ls 2 t '{/Symbol b} = 0.05'

