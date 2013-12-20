set term post eps enh color size 8.5cm,6cm

set style line 11 lc 1 lw 3 pt 2 lt 1
set style line 12 lc 2 lw 3 pt 2 lt 1
set style line 13 lc 3 lw 3 pt 2 lt 1
set style line 14 lc 4 lw 3 pt 2 lt 1
set style line 15 lc 5 lw 3 pt 2 lt 1

set style line 21 lc 1 lw 3 pt 2 lt 2
set style line 22 lc 2 lw 3 pt 2 lt 2
set style line 23 lc 3 lw 3 pt 2 lt 2
set style line 24 lc 4 lw 3 pt 2 lt 2
set style line 25 lc 5 lw 3 pt 2 lt 2

set style line 31 lc 1 lw 3 pt 2 lt 3
set style line 32 lc 2 lw 3 pt 2 lt 3
set style line 33 lc 3 lw 3 pt 2 lt 3
set style line 34 lc 4 lw 3 pt 2 lt 3
set style line 35 lc 5 lw 3 pt 2 lt 3

set style line 41 lc 1 lw 3 pt 2 lt 4
set style line 42 lc 2 lw 3 pt 2 lt 4
set style line 43 lc 3 lw 3 pt 2 lt 4
set style line 44 lc 4 lw 3 pt 2 lt 4
set style line 45 lc 5 lw 3 pt 2 lt 4

set style line 20 lc 0 lw 3 pt 2 lt 2


set out 'fig-pot.eps'

set xlabel '{/Symbol f} [ deg. ]'
#set ylabel 'V_{{/Symbol f}}({/Symbol f}) + V_{ctrl}({/Symbol f})'
set ylabel 'Dihedral angle energy'
set xtics 20
set mxtics 2
set mytics 5
set grid

pl [40:150] 0 ls 20 not,\
   'table_orig.xvg' u 1:($2 - 5.920000e+00) w l ls 11 t 'V_{{/Symbol f}}({/Symbol f})',\
   'ctrl.1.00.xvg' u 1:($2 - 2.349671e+00) w l ls 12 t 'V_{{/Symbol f}}({/Symbol f}) + V_{ctrl}({/Symbol f}), {/Symbol l}=1.00 ps',\
   'ctrl.0.50.xvg' u 1:($2 - 8.613114e-01) w l ls 13 t 'V_{{/Symbol f}}({/Symbol f}) + V_{ctrl}({/Symbol f}), {/Symbol l}=0.50 ps',\
   'ctrl.0.20.xvg' u 1:($2 - -2.410373e+00) w l ls 14 t 'V_{{/Symbol f}}({/Symbol f}) + V_{ctrl}({/Symbol f}), {/Symbol l}=0.20 ps',\
   'ctrl.0.10.xvg' u 1:($2 - -1.102120e+01) w l ls 15 t 'V_{{/Symbol f}}({/Symbol f}) + V_{ctrl}({/Symbol f}), {/Symbol l}=0.10 ps'