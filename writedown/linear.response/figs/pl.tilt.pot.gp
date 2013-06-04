set term post eps enh color solid font 16 size 8.5cm,6cm
set out 'fig-tilt-pot.eps'

set style line 10 lt 0 lc 0 lw 3 pt 3
set style line 1 lt 1 lc 1 lw 3 pt 3
set style line 2 lt 1 lc 2 lw 3 pt 3
set style line 3 lt 1 lc 3 lw 3 pt 3
set style line 4 lt 1 lc 4 lw 3 pt 3

set xlabel "x [nm]"
set ylabel "U(x) + {/Symbol e} F_e(t) V(x) [kJ/mol]"
set mxtics 5
set mytics 5


pl[-2:2][-10:20] 0 ls 0 not,\
0.5 * 8 * (x*x - 1*1)**2 + -x * 0 ls 1 t '{/Symbol e} = 0',\
0.5 * 8 * (x*x - 1*1)**2 + -x * 2 ls 2 t '{/Symbol e} = 2',\
0.5 * 8 * (x*x - 1*1)**2 + -x * 4 ls 3 t '{/Symbol e} = 4'

