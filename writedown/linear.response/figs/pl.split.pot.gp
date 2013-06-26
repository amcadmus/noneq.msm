set term post eps enh color solid font 16 size 8.5cm,6cm
set out 'fig-split-pot.eps'

set style line 10 lt 0 lc 0 lw 3 pt 3
set style line 1 lt 1 lc 1 lw 3 pt 3
set style line 2 lt 1 lc 2 lw 3 pt 3
set style line 3 lt 1 lc 3 lw 3 pt 3
set style line 4 lt 1 lc 4 lw 3 pt 3

set xlabel "q [nm]"
set ylabel "Energy [kJ/mol]"
set mxtics 5
set mytics 2

sigma2=0.16**2

set key top left
set key spacing 1.5

pl[-1.5:1.5][-1:8] 0 ls 0 not, \
0.0 * 1./sqrt(2*pi*sigma2) * exp(-x**2/(2.*sigma2)) + 4*x**2 ls 1 t'U(q)',\
1.0 * 1./sqrt(2*pi*sigma2) * exp(-x**2/(2.*sigma2)) + 4*x**2 ls 2 t'U(q) + F_e(t_c)V(q)',\
2.0 * 1./sqrt(2*pi*sigma2) * exp(-x**2/(2.*sigma2)) + 4*x**2 ls 3 t'U(q) + [F_e(t_c) + {/Symbol e}{/Symbol D}F_e(t_c)]V(q)'


# pl[-2:2][-10:20] 0 ls 0 not,\
# 0.5 * 8 * (x*x - 1*1)**2 + -x * 0 ls 1 t 'F_e = 0',\
# 0.5 * 8 * (x*x - 1*1)**2 + -x * 2 ls 2 t 'F_e = 2',\
# 0.5 * 8 * (x*x - 1*1)**2 + -x * 4 ls 3 t 'F_e = 4'

