set term post eps enh color solid font 16 size 8.5cm,6cm
set out 'fig-1d.eps'

set style line 10 lt 1 lc 0 lw 3 pt 3
set style line 1 lt 1 lc 1 lw 3 pt 3
set style line 2 lt 1 lc 2 lw 3 pt 3
set style line 3 lt 1 lc 3 lw 3 pt 3
set style line 4 lt 1 lc 4 lw 3 pt 3

set key top right

# pl 0 ls 10 not, \
#    'indicator.x.00000.00.out' w l ls 1 t 't=0', \
#    'indicator.linear.x.00000.00.out' w p ls 1 not,\
#    'indicator.x.00030.00.out' w l ls 2 t 't=30', \
#    'indicator.linear.x.00030.00.out' w p ls 2 not, \
#    'indicator.x.00100.00.out' w l ls 3 t 't=100', \
#    'indicator.linear.x.00100.00.out' w p ls 3 not,\
#    'indicator.x.00150.00.out' w l ls 4 t 't=150', \
#    'indicator.linear.x.00150.00.out' w p ls 4 not

pl 0 ls 10 not, \
'./indicator.x.00030.00.out' w l ls 1, \
'./indicator.linear.corr.x.00030.00.out' w l ls 2,\
'./indicator.linear.x.00030.00.out' w l ls 3,\
'./quenchCorr.x.out' u 1:($2/1) w l ls 4
