set term post eps enh color font 16 size 8.5cm,6cm

# color at http://html-color-codes.info/
set style line 11 lt 1 lw 4 pt 7 linecolor rgb "#F5A9A9"
set style line 12 lt 1 lw 4 pt 7 linecolor rgb "#ff0000"
set style line 13 lt 4 lc 1 lw 4 pt 7
# set style line 21 lt 1 lw 4 pt 7 linecolor rgb "#bbffbb"
# set style line 22 lt 1 lw 4 pt 7 linecolor rgb "#00ff66"
set style line 21 lt 1 lw 4 pt 7 linecolor rgb "#A9F5A9"
set style line 22 lt 1 lw 4 pt 7 linecolor rgb "#04B404"
set style line 23 lt 4 lc 2 lw 4 pt 7
set style line 31 lt 1 lw 4 pt 7 linecolor rgb "#8181F7"
set style line 32 lt 1 lw 4 pt 7 linecolor rgb "#0000ff"
set style line 33 lt 4 lc 3 lw 4 pt 7
set style line 41 lt 1 lw 4 pt 7 linecolor rgb "#F781F3"
set style line 42 lt 1 lw 4 pt 7 linecolor rgb "#B404AE"
set style line 51 lt 1 lw 4 pt 7 linecolor rgb "#A9F5F2"
set style line 52 lt 1 lw 4 pt 7 linecolor rgb "#01DFD7"

# set style line 11 lt 1 lw 4 pt 7 linecolor rgb "#ffbbbb"
# set style line 12 lt 1 lw 4 pt 7 linecolor rgb "#ff0000"
# set style line 13 lt 4 lc 1 lw 4 pt 7
# # set style line 21 lt 1 lw 4 pt 7 linecolor rgb "#bbffbb"
# # set style line 22 lt 1 lw 4 pt 7 linecolor rgb "#00ff66"
# set style line 21 lt 1 lw 4 pt 7 linecolor rgb "#99ff99"
# set style line 22 lt 1 lw 4 pt 7 linecolor rgb "#33cc33"
# set style line 23 lt 4 lc 2 lw 4 pt 7
# set style line 31 lt 1 lw 4 pt 7 linecolor rgb "#bbbbff"
# set style line 32 lt 1 lw 4 pt 7 linecolor rgb "#0000ff"
# set style line 33 lt 4 lc 3 lw 4 pt 7
# set style line 41 lt 1 lw 4 pt 7 linecolor rgb "#ffccff"
# set style line 42 lt 1 lw 4 pt 7 linecolor rgb "#cc00ff"
# set style line 51 lt 1 lw 4 pt 7 linecolor rgb "#bbffff"
# set style line 52 lt 1 lw 4 pt 7 linecolor rgb "#00ffff"

set xlabel 't [ ps ]'
set ylabel 'Probability of conformations'
# set yrange [0:.6]
set mxtics 5
set ytics .1
set mytics 5
set grid

set out 'fig-meta02.eps'
pl \
   'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale01.0/set/metastable.out' u 1:2 w l ls 11 not,\
   '' u 1:3 w l ls 21 not,\
   '' u 1:4 w l ls 31 not,\
   '' u 1:5 w l ls 41 not,\
   '' u 1:6 w l ls 51 not,\
   'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale02.0/set/metastable.out' u ($1*2):2  w l ls 12 not,\
   '' u ($1*2):3  w l ls 22 not,\
   '' u ($1*2):4  w l ls 32 not,\
   '' u ($1*2):5  w l ls 42 not,\
   '' u ($1*2):6  w l ls 52 not

set out 'fig-meta04.eps'
pl \
   'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale01.0/set/metastable.out' u 1:2 w l ls 11 not,\
   '' u 1:3 w l ls 21 not,\
   '' u 1:4 w l ls 31 not,\
   '' u 1:5 w l ls 41 not,\
   '' u 1:6 w l ls 51 not,\
   'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale04.0/set/metastable.out' u ($1*4):2  w l ls 12 not,\
   '' u ($1*4):3  w l ls 22 not,\
   '' u ($1*4):4  w l ls 32 not,\
   '' u ($1*4):5  w l ls 42 not,\
   '' u ($1*4):6  w l ls 52 not

set out 'fig-meta04-1.eps'
pl \
   'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale01.0/set/metastable.out' u 1:($2+$3) w l ls 11 not,\
   '' u 1:($4+$5) w l ls 31 not,\
   '' u 1:6 w l ls 51 not,\
   'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale04.0/set/metastable.out' u ($1*4):($2+$3)  w l ls 12 not,\
   '' u ($1*4):($4+$5)  w l ls 32 not,\
   '' u ($1*4):6  w l ls 52 not

   
set out 'fig-meta05.eps'
pl \
   'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale01.0/set/metastable.out' u 1:2 w l ls 11 not,\
   '' u 1:3 w l ls 21 not,\
   '' u 1:4 w l ls 31 not,\
   '' u 1:5 w l ls 41 not,\
   '' u 1:6 w l ls 51 not,\
   'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale05.0/set/metastable.out' u ($1*5):2  w l ls 12 not,\
   '' u ($1*5):3  w l ls 22 not,\
   '' u ($1*5):4  w l ls 32 not,\
   '' u ($1*5):5  w l ls 42 not,\
   '' u ($1*5):6  w l ls 52 not

set out 'fig-meta05-1.eps'
pl \
   'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale01.0/set/metastable.out' u 1:($2+$3) w l ls 11 not,\
   '' u 1:($4+$5) w l ls 31 not,\
   '' u 1:6 w l ls 51 not,\
   'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale05.0/set/metastable.out' u ($1*5):($2+$3)  w l ls 12 not,\
   '' u ($1*5):($4+$5)  w l ls 32 not,\
   '' u ($1*5):6  w l ls 52 not

   
set out 'fig-meta06.eps'
pl \
   'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale01.0/set/metastable.out' u 1:2 w l ls 11 not,\
   '' u 1:3 w l ls 21 not,\
   '' u 1:4 w l ls 31 not,\
   '' u 1:5 w l ls 41 not,\
   '' u 1:6 w l ls 51 not,\
   'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale06.0/set/metastable.out' u ($1*6):2  w l ls 12 not,\
   '' u ($1*6):3  w l ls 22 not,\
   '' u ($1*6):4  w l ls 32 not,\
   '' u ($1*6):5  w l ls 42 not,\
   '' u ($1*6):6  w l ls 52 not


set out 'fig-meta06-1.eps'
pl \
   'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale01.0/set/metastable.out' u 1:($2+$3) w l ls 11 not,\
   '' u 1:($4+$5) w l ls 31 not,\
   '' u 1:6 w l ls 51 not,\
   'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale06.0/set/metastable.out' u ($1*6):($2+$3)  w l ls 12 not,\
   '' u ($1*6):($4+$5)  w l ls 32 not,\
   '' u ($1*6):6  w l ls 52 not


   
set out 'fig-meta08.eps'
pl \
   'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale01.0/set/metastable.out' u 1:2 w l ls 11 not,\
   '' u 1:3 w l ls 21 not,\
   '' u 1:4 w l ls 31 not,\
   '' u 1:5 w l ls 41 not,\
   '' u 1:6 w l ls 51 not,\
   'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale08.0/set/metastable.out' u ($1*8):2  w l ls 12 not,\
   '' u ($1*8):3  w l ls 22 not,\
   '' u ($1*8):4  w l ls 32 not,\
   '' u ($1*8):5  w l ls 42 not,\
   '' u ($1*8):6  w l ls 52 not


set out 'fig-meta08-1.eps'
pl \
   'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale01.0/set/metastable.out' u 1:($2+$3) w l ls 11 not,\
   '' u 1:($4+$5) w l ls 31 not,\
   '' u 1:6 w l ls 51 not,\
   'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale08.0/set/metastable.out' u ($1*8):($2+$3)  w l ls 12 not,\
   '' u ($1*8):($4+$5)  w l ls 32 not,\
   '' u ($1*8):6  w l ls 52 not
   