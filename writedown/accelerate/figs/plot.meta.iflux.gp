set lmargin 4.5
set rmargin 2
unset tmargin
unset bmargin
set term post eps enh color font 14 size 4.2cm,3cm

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

# set style line 11 lc 1 lw 3 pt 2 lt 1
# set style line 12 lc 2 lw 3 pt 2 lt 1
# set style line 13 lc 3 lw 3 pt 2 lt 1
# set style line 14 lc 4 lw 3 pt 2 lt 1
# set style line 15 lc 5 lw 3 pt 2 lt 1

# set style line 21 lc 1 lw 3 pt 2 lt 2
# set style line 22 lc 2 lw 3 pt 2 lt 2
# set style line 23 lc 3 lw 3 pt 2 lt 2
# set style line 24 lc 4 lw 3 pt 2 lt 2
# set style line 25 lc 5 lw 3 pt 2 lt 2

# set style line 31 lc 1 lw 3 pt 2 lt 3
# set style line 32 lc 2 lw 3 pt 2 lt 3
# set style line 33 lc 3 lw 3 pt 2 lt 3
# set style line 34 lc 4 lw 3 pt 2 lt 3
# set style line 35 lc 5 lw 3 pt 2 lt 3

# set style line 41 lc 1 lw 3 pt 2 lt 4
# set style line 42 lc 2 lw 3 pt 2 lt 4
# set style line 43 lc 3 lw 3 pt 2 lt 4
# set style line 44 lc 4 lw 3 pt 2 lt 4
# set style line 45 lc 5 lw 3 pt 2 lt 4

set style line 20 lc 0 lw 3 pt 2 lt 2


unset xrange
set xtics 100
set mxtics 5
# set xlabel 't [ ps ]'
# set yrange [-0.5:0.5]
set mytics 2
set format y "%.1f"
set grid

################################################################################
unset label
set label 'Q_{J,A_1}, T = 40 ps' at 20,0.3
set yrange [-0.5:0.5]
set out 'fig-meta-iflux-02-1.eps'
pl\
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale01.0/set/meta.flux.inte.out' \
u 1:2 ls 11 w l not, \
'' u 1:3 ls 21 w l not, \
'' u 1:4 ls 31 w l not, \
'' u 1:5 ls 41 w l not, \
'' u 1:6 ls 51 w l not, \
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale02.0/set/meta.flux.inte.out' \
u ($1*2):2 ls 12 w l not, \
'' u ($1*2):3  w l ls 22 not,\
'' u ($1*2):4  w l ls 32 not,\
'' u ($1*2):5  w l ls 42 not,\
'' u ($1*2):6  w l ls 52 not

set out 'fig-meta-iflux-01-check-1.eps'
pl\
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale01.0/set/meta.flux.inte.out' \
u 1:2 ls 11 w l not, \
'' u 1:3 ls 21 w l not, \
'' u 1:4 ls 31 w l not, \
'' u 1:5 ls 41 w l not, \
'' u 1:6 ls 51 w l not, \
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale01.0.dtcomp/set/meta.flux.inte.out' \
u ($1*1):2 ls 12 w l not, \
'' u ($1*1):3  w l ls 22 not,\
'' u ($1*1):4  w l ls 32 not,\
'' u ($1*1):5  w l ls 42 not,\
'' u ($1*1):6  w l ls 52 not

set out 'fig-meta-iflux-04-1.eps'
pl\
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale01.0/set/meta.flux.inte.out' \
u 1:2 ls 11 w l not, \
'' u 1:3 ls 21 w l not, \
'' u 1:4 ls 31 w l not, \
'' u 1:5 ls 41 w l not, \
'' u 1:6 ls 51 w l not, \
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale04.0/set/meta.flux.inte.out' \
u ($1*4):2 ls 12 w l not, \
'' u ($1*4):3  w l ls 22 not,\
'' u ($1*4):4  w l ls 32 not,\
'' u ($1*4):5  w l ls 42 not,\
'' u ($1*4):6  w l ls 52 not

set out 'fig-meta-iflux-05-1.eps'
pl\
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale01.0/set/meta.flux.inte.out' \
u 1:2 ls 11 w l not, \
'' u 1:3 ls 21 w l not, \
'' u 1:4 ls 31 w l not, \
'' u 1:5 ls 41 w l not, \
'' u 1:6 ls 51 w l not, \
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale05.0/set/meta.flux.inte.out' \
u ($1*5):2 ls 12 w l not, \
'' u ($1*5):3  w l ls 22 not,\
'' u ($1*5):4  w l ls 32 not,\
'' u ($1*5):5  w l ls 42 not,\
'' u ($1*5):6  w l ls 52 not

set out 'fig-meta-iflux-08-1.eps'
pl\
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale01.0/set/meta.flux.inte.out' \
u 1:2 ls 11 w l not, \
'' u 1:3 ls 21 w l not, \
'' u 1:4 ls 31 w l not, \
'' u 1:5 ls 41 w l not, \
'' u 1:6 ls 51 w l not, \
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale08.0/set/meta.flux.inte.out' \
u ($1*8):2 ls 12 w l not, \
'' u ($1*8):3  w l ls 22 not,\
'' u ($1*8):4  w l ls 32 not,\
'' u ($1*8):5  w l ls 42 not,\
'' u ($1*8):6  w l ls 52 not



################################################################################
unset label
set label 'Q_{J,A_2}, T = 40 ps' at 20,0.3
set yrange [-0.5:0.5]
set out 'fig-meta-iflux-02-2.eps'
pl\
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale01.0/set/meta.flux.inte.out' \
u 1:7 ls 11 w l not, \
'' u 1:8 ls 21 w l not, \
'' u 1:9 ls 31 w l not, \
'' u 1:10 ls 41 w l not, \
'' u 1:11 ls 51 w l not, \
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale02.0/set/meta.flux.inte.out' \
u ($1*2):7 ls 12 w l not, \
'' u ($1*2):8  w l ls 22 not,\
'' u ($1*2):9  w l ls 32 not,\
'' u ($1*2):10  w l ls 42 not,\
'' u ($1*2):11  w l ls 52 not

set out 'fig-meta-iflux-01-check-2.eps'
pl\
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale01.0/set/meta.flux.inte.out' \
u 1:7 ls 11 w l not, \
'' u 1:8 ls 21 w l not, \
'' u 1:9 ls 31 w l not, \
'' u 1:10 ls 41 w l not, \
'' u 1:11 ls 51 w l not, \
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale01.0.dtcomp/set/meta.flux.inte.out' \
u ($1*1):7 ls 12 w l not, \
'' u ($1*1):8  w l ls 22 not,\
'' u ($1*1):9  w l ls 32 not,\
'' u ($1*1):10  w l ls 42 not,\
'' u ($1*1):11  w l ls 52 not

set out 'fig-meta-iflux-04-2.eps'
pl\
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale01.0/set/meta.flux.inte.out' \
u 1:7 ls 11 w l not, \
'' u 1:8 ls 21 w l not, \
'' u 1:9 ls 31 w l not, \
'' u 1:10 ls 41 w l not, \
'' u 1:11 ls 51 w l not, \
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale04.0/set/meta.flux.inte.out' \
u ($1*4):7 ls 12 w l not, \
'' u ($1*4):8  w l ls 22 not,\
'' u ($1*4):9  w l ls 32 not,\
'' u ($1*4):10  w l ls 42 not,\
'' u ($1*4):11  w l ls 52 not

set out 'fig-meta-iflux-05-2.eps'
pl\
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale01.0/set/meta.flux.inte.out' \
u 1:7 ls 11 w l not, \
'' u 1:8 ls 21 w l not, \
'' u 1:9 ls 31 w l not, \
'' u 1:10 ls 41 w l not, \
'' u 1:11 ls 51 w l not, \
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale05.0/set/meta.flux.inte.out' \
u ($1*5):7 ls 12 w l not, \
'' u ($1*5):8  w l ls 22 not,\
'' u ($1*5):9  w l ls 32 not,\
'' u ($1*5):10  w l ls 42 not,\
'' u ($1*5):11  w l ls 52 not

set out 'fig-meta-iflux-08-2.eps'
pl\
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale01.0/set/meta.flux.inte.out' \
u 1:7 ls 11 w l not, \
'' u 1:8 ls 21 w l not, \
'' u 1:9 ls 31 w l not, \
'' u 1:10 ls 41 w l not, \
'' u 1:11 ls 51 w l not, \
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale08.0/set/meta.flux.inte.out' \
u ($1*8):7 ls 12 w l not, \
'' u ($1*8):8  w l ls 22 not,\
'' u ($1*8):9  w l ls 32 not,\
'' u ($1*8):10  w l ls 42 not,\
'' u ($1*8):11  w l ls 52 not


################################################################################
unset label
set label 'Q_{J,B_1}, T = 40 ps' at 20,0.3
set yrange [-0.5:0.5]
set out 'fig-meta-iflux-02-3.eps'
pl\
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale01.0/set/meta.flux.inte.out' \
u 1:12 ls 11 w l not, \
'' u 1:13 ls 21 w l not, \
'' u 1:14 ls 31 w l not, \
'' u 1:15 ls 41 w l not, \
'' u 1:16 ls 51 w l not, \
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale02.0/set/meta.flux.inte.out' \
u ($1*2):12 ls 12 w l not, \
'' u ($1*2):13  w l ls 22 not,\
'' u ($1*2):14  w l ls 32 not,\
'' u ($1*2):15  w l ls 42 not,\
'' u ($1*2):16  w l ls 52 not

set out 'fig-meta-iflux-01-check-3.eps'
pl\
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale01.0/set/meta.flux.inte.out' \
u 1:12 ls 11 w l not, \
'' u 1:13 ls 21 w l not, \
'' u 1:14 ls 31 w l not, \
'' u 1:15 ls 41 w l not, \
'' u 1:16 ls 51 w l not, \
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale01.0.dtcomp/set/meta.flux.inte.out' \
u ($1*1):12 ls 12 w l not, \
'' u ($1*1):13  w l ls 22 not,\
'' u ($1*1):14  w l ls 32 not,\
'' u ($1*1):15  w l ls 42 not,\
'' u ($1*1):16  w l ls 52 not

set out 'fig-meta-iflux-04-3.eps'
pl\
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale01.0/set/meta.flux.inte.out' \
u 1:12 ls 11 w l not, \
'' u 1:13 ls 21 w l not, \
'' u 1:14 ls 31 w l not, \
'' u 1:15 ls 41 w l not, \
'' u 1:16 ls 51 w l not, \
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale04.0/set/meta.flux.inte.out' \
u ($1*4):12 ls 12 w l not, \
'' u ($1*4):13  w l ls 22 not,\
'' u ($1*4):14  w l ls 32 not,\
'' u ($1*4):15  w l ls 42 not,\
'' u ($1*4):16  w l ls 52 not

set out 'fig-meta-iflux-05-3.eps'
pl\
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale01.0/set/meta.flux.inte.out' \
u 1:12 ls 11 w l not, \
'' u 1:13 ls 21 w l not, \
'' u 1:14 ls 31 w l not, \
'' u 1:15 ls 41 w l not, \
'' u 1:16 ls 51 w l not, \
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale05.0/set/meta.flux.inte.out' \
u ($1*5):12 ls 12 w l not, \
'' u ($1*5):13  w l ls 22 not,\
'' u ($1*5):14  w l ls 32 not,\
'' u ($1*5):15  w l ls 42 not,\
'' u ($1*5):16  w l ls 52 not

set out 'fig-meta-iflux-08-3.eps'
pl\
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale01.0/set/meta.flux.inte.out' \
u 1:12 ls 11 w l not, \
'' u 1:13 ls 21 w l not, \
'' u 1:14 ls 31 w l not, \
'' u 1:15 ls 41 w l not, \
'' u 1:16 ls 51 w l not, \
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale08.0/set/meta.flux.inte.out' \
u ($1*8):12 ls 12 w l not, \
'' u ($1*8):13  w l ls 22 not,\
'' u ($1*8):14  w l ls 32 not,\
'' u ($1*8):15  w l ls 42 not,\
'' u ($1*8):16  w l ls 52 not


################################################################################
unset label
set label 'Q_{J,B_2}, T = 40 ps' at 20,0.3
set yrange [-0.5:0.5]
set out 'fig-meta-iflux-02-4.eps'
pl\
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale01.0/set/meta.flux.inte.out' \
u 1:17 ls 11 w l not, \
'' u 1:18 ls 21 w l not, \
'' u 1:19 ls 31 w l not, \
'' u 1:20 ls 41 w l not, \
'' u 1:21 ls 51 w l not, \
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale02.0/set/meta.flux.inte.out' \
u ($1*2):17 ls 12 w l not, \
'' u ($1*2):18  w l ls 22 not,\
'' u ($1*2):19  w l ls 32 not,\
'' u ($1*2):20  w l ls 42 not,\
'' u ($1*2):21  w l ls 52 not

set out 'fig-meta-iflux-01-check-4.eps'
pl\
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale01.0/set/meta.flux.inte.out' \
u 1:17 ls 11 w l not, \
'' u 1:18 ls 21 w l not, \
'' u 1:19 ls 31 w l not, \
'' u 1:20 ls 41 w l not, \
'' u 1:21 ls 51 w l not, \
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale01.0.dtcomp/set/meta.flux.inte.out' \
u ($1*1):17 ls 12 w l not, \
'' u ($1*1):18  w l ls 22 not,\
'' u ($1*1):19  w l ls 32 not,\
'' u ($1*1):20  w l ls 42 not,\
'' u ($1*1):21  w l ls 52 not

set out 'fig-meta-iflux-04-4.eps'
pl\
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale01.0/set/meta.flux.inte.out' \
u 1:17 ls 11 w l not, \
'' u 1:18 ls 21 w l not, \
'' u 1:19 ls 31 w l not, \
'' u 1:20 ls 41 w l not, \
'' u 1:21 ls 51 w l not, \
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale04.0/set/meta.flux.inte.out' \
u ($1*4):17 ls 12 w l not, \
'' u ($1*4):18  w l ls 22 not,\
'' u ($1*4):19  w l ls 32 not,\
'' u ($1*4):20  w l ls 42 not,\
'' u ($1*4):21  w l ls 52 not

set out 'fig-meta-iflux-05-4.eps'
pl\
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale01.0/set/meta.flux.inte.out' \
u 1:17 ls 11 w l not, \
'' u 1:18 ls 21 w l not, \
'' u 1:19 ls 31 w l not, \
'' u 1:20 ls 41 w l not, \
'' u 1:21 ls 51 w l not, \
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale05.0/set/meta.flux.inte.out' \
u ($1*5):17 ls 12 w l not, \
'' u ($1*5):18  w l ls 22 not,\
'' u ($1*5):19  w l ls 32 not,\
'' u ($1*5):20  w l ls 42 not,\
'' u ($1*5):21  w l ls 52 not

set out 'fig-meta-iflux-08-4.eps'
pl\
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale01.0/set/meta.flux.inte.out' \
u 1:17 ls 11 w l not, \
'' u 1:18 ls 21 w l not, \
'' u 1:19 ls 31 w l not, \
'' u 1:20 ls 41 w l not, \
'' u 1:21 ls 51 w l not, \
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale08.0/set/meta.flux.inte.out' \
u ($1*8):17 ls 12 w l not, \
'' u ($1*8):18  w l ls 22 not,\
'' u ($1*8):19  w l ls 32 not,\
'' u ($1*8):20  w l ls 42 not,\
'' u ($1*8):21  w l ls 52 not



################################################################################
unset label
set label 'Q_{J,C}, T = 40 ps' at 20,0.3
set yrange [-0.5:0.5]
set out 'fig-meta-iflux-02-5.eps'
pl\
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale01.0/set/meta.flux.inte.out' \
u 1:22 ls 11 w l not, \
'' u 1:23 ls 21 w l not, \
'' u 1:24 ls 31 w l not, \
'' u 1:25 ls 41 w l not, \
'' u 1:26 ls 51 w l not, \
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale02.0/set/meta.flux.inte.out' \
u ($1*2):22 ls 12 w l not, \
'' u ($1*2):23  w l ls 22 not,\
'' u ($1*2):24  w l ls 32 not,\
'' u ($1*2):25  w l ls 42 not,\
'' u ($1*2):26  w l ls 52 not

set out 'fig-meta-iflux-01-check-5.eps'
pl\
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale01.0/set/meta.flux.inte.out' \
u 1:22 ls 11 w l not, \
'' u 1:23 ls 21 w l not, \
'' u 1:24 ls 31 w l not, \
'' u 1:25 ls 41 w l not, \
'' u 1:26 ls 51 w l not, \
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale01.0.dtcomp/set/meta.flux.inte.out' \
u ($1*1):22 ls 12 w l not, \
'' u ($1*1):23  w l ls 22 not,\
'' u ($1*1):24  w l ls 32 not,\
'' u ($1*1):25  w l ls 42 not,\
'' u ($1*1):26  w l ls 52 not

set out 'fig-meta-iflux-04-5.eps'
pl\
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale01.0/set/meta.flux.inte.out' \
u 1:22 ls 11 w l not, \
'' u 1:23 ls 21 w l not, \
'' u 1:24 ls 31 w l not, \
'' u 1:25 ls 41 w l not, \
'' u 1:26 ls 51 w l not, \
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale04.0/set/meta.flux.inte.out' \
u ($1*4):22 ls 12 w l not, \
'' u ($1*4):23  w l ls 22 not,\
'' u ($1*4):24  w l ls 32 not,\
'' u ($1*4):25  w l ls 42 not,\
'' u ($1*4):26  w l ls 52 not

set out 'fig-meta-iflux-05-5.eps'
pl\
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale01.0/set/meta.flux.inte.out' \
u 1:22 ls 11 w l not, \
'' u 1:23 ls 21 w l not, \
'' u 1:24 ls 31 w l not, \
'' u 1:25 ls 41 w l not, \
'' u 1:26 ls 51 w l not, \
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale05.0/set/meta.flux.inte.out' \
u ($1*5):22 ls 12 w l not, \
'' u ($1*5):23  w l ls 22 not,\
'' u ($1*5):24  w l ls 32 not,\
'' u ($1*5):25  w l ls 42 not,\
'' u ($1*5):26  w l ls 52 not

set out 'fig-meta-iflux-08-5.eps'
pl\
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale01.0/set/meta.flux.inte.out' \
u 1:22 ls 11 w l not, \
'' u 1:23 ls 21 w l not, \
'' u 1:24 ls 31 w l not, \
'' u 1:25 ls 41 w l not, \
'' u 1:26 ls 51 w l not, \
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale08.0/set/meta.flux.inte.out' \
u ($1*8):22 ls 12 w l not, \
'' u ($1*8):23  w l ls 22 not,\
'' u ($1*8):24  w l ls 32 not,\
'' u ($1*8):25  w l ls 42 not,\
'' u ($1*8):26  w l ls 52 not



################################################################################
unset label
set yrange [-0.3:0.3]
set label 'Q_{J,A}, T = 40 ps' at 20,0.25
set out 'fig-meta-iflux-04-1-1.eps'
pl\
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale01.0/set/meta.flux.inte.out' \
u 1:($2+$3+$7+$8) ls 11 w l not, \
'' u 1:($4+$5+$9+$10) ls 31 w l not, \
'' u 1:($6+$11) ls 51 w l not, \
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale04.0/set/meta.flux.inte.out' \
u ($1*4):($2+$3+$7+$8) ls 12 w l not, \
'' u ($1*4):($4+$5+$9+$10)  w l ls 32 not,\
'' u ($1*4):($6+$11)  w l ls 52 not


unset label
set label 'Q_{J,B}, T = 40 ps' at 20,0.25
set out 'fig-meta-iflux-04-1-2.eps'
pl\
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale01.0/set/meta.flux.inte.out' \
u 1:($12+$13+$17+$18) ls 11 w l not, \
'' u 1:($14+$15+$19+$20) ls 31 w l not, \
'' u 1:($16+$21) ls 51 w l not, \
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale04.0/set/meta.flux.inte.out' \
u ($1*4):($12+$13+$17+$18) ls 12 w l not, \
'' u ($1*4):($14+$15+$19+$20)  w l ls 32 not,\
'' u ($1*4):($16+$21)  w l ls 52 not


unset label
set label 'Q_{J,C}, T = 40 ps' at 20,0.25
set out 'fig-meta-iflux-04-1-3.eps'
pl\
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale01.0/set/meta.flux.inte.out' \
u 1:($22+$23) ls 11 w l not, \
'' u 1:($24+$25) ls 31 w l not, \
'' u 1:($26) ls 51 w l not, \
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale04.0/set/meta.flux.inte.out' \
u ($1*4):($22+$23) ls 12 w l not, \
'' u ($1*4):($24+$25)  w l ls 32 not,\
'' u ($1*4):($26)  w l ls 52 not


unset label
set label 'Q_{J,A}, T = 40 ps' at 20,0.25
set out 'fig-meta-iflux-02-1-1.eps'
pl\
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale01.0/set/meta.flux.inte.out' \
u 1:($2+$3+$7+$8) ls 11 w l not, \
'' u 1:($4+$5+$9+$10) ls 31 w l not, \
'' u 1:($6+$11) ls 51 w l not, \
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale02.0/set/meta.flux.inte.out' \
u ($1*2):($2+$3+$7+$8) ls 12 w l not, \
'' u ($1*2):($4+$5+$9+$10)  w l ls 32 not,\
'' u ($1*2):($6+$11)  w l ls 52 not


unset label
set label 'Q_{J,B}, T = 40 ps' at 20,0.25
set out 'fig-meta-iflux-02-1-2.eps'
pl\
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale01.0/set/meta.flux.inte.out' \
u 1:($12+$13+$17+$18) ls 11 w l not, \
'' u 1:($14+$15+$19+$20) ls 31 w l not, \
'' u 1:($16+$21) ls 51 w l not, \
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale02.0/set/meta.flux.inte.out' \
u ($1*2):($12+$13+$17+$18) ls 12 w l not, \
'' u ($1*2):($14+$15+$19+$20)  w l ls 32 not,\
'' u ($1*2):($16+$21)  w l ls 52 not


unset label
set label 'Q_{J,C}, T = 40 ps' at 20,0.25
set out 'fig-meta-iflux-02-1-3.eps'
pl\
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale01.0/set/meta.flux.inte.out' \
u 1:($22+$23) ls 11 w l not, \
'' u 1:($24+$25) ls 31 w l not, \
'' u 1:($26) ls 51 w l not, \
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale02.0/set/meta.flux.inte.out' \
u ($1*2):($22+$23) ls 12 w l not, \
'' u ($1*2):($24+$25)  w l ls 32 not,\
'' u ($1*2):($26)  w l ls 52 not


unset label
set label 'Q_{J,A}, T = 40 ps' at 20,0.25
set out 'fig-meta-iflux-05-1-1.eps'
pl\
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale01.0/set/meta.flux.inte.out' \
u 1:($2+$3+$7+$8) ls 11 w l not, \
'' u 1:($4+$5+$9+$10) ls 31 w l not, \
'' u 1:($6+$11) ls 51 w l not, \
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale05.0/set/meta.flux.inte.out' \
u ($1*5):($2+$3+$7+$8) ls 12 w l not, \
'' u ($1*5):($4+$5+$9+$10)  w l ls 32 not,\
'' u ($1*5):($6+$11)  w l ls 52 not


unset label
set label 'Q_{J,B}, T = 40 ps' at 20,0.25
set out 'fig-meta-iflux-05-1-2.eps'
pl\
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale01.0/set/meta.flux.inte.out' \
u 1:($12+$13+$17+$18) ls 11 w l not, \
'' u 1:($14+$15+$19+$20) ls 31 w l not, \
'' u 1:($16+$21) ls 51 w l not, \
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale05.0/set/meta.flux.inte.out' \
u ($1*5):($12+$13+$17+$18) ls 12 w l not, \
'' u ($1*5):($14+$15+$19+$20)  w l ls 32 not,\
'' u ($1*5):($16+$21)  w l ls 52 not


unset label
set label 'Q_{J,C}, T = 40 ps' at 20,0.25
set out 'fig-meta-iflux-05-1-3.eps'
pl\
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale01.0/set/meta.flux.inte.out' \
u 1:($22+$23) ls 11 w l not, \
'' u 1:($24+$25) ls 31 w l not, \
'' u 1:($26) ls 51 w l not, \
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale05.0/set/meta.flux.inte.out' \
u ($1*5):($22+$23) ls 12 w l not, \
'' u ($1*5):($24+$25)  w l ls 32 not,\
'' u ($1*5):($26)  w l ls 52 not



unset label
set label 'Q_{J,A}, T = 40 ps' at 20,0.25
set out 'fig-meta-iflux-08-1-1.eps'
pl\
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale01.0/set/meta.flux.inte.out' \
u 1:($2+$3+$7+$8) ls 11 w l not, \
'' u 1:($4+$5+$9+$10) ls 31 w l not, \
'' u 1:($6+$11) ls 51 w l not, \
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale08.0/set/meta.flux.inte.out' \
u ($1*8):($2+$3+$7+$8) ls 12 w l not, \
'' u ($1*8):($4+$5+$9+$10)  w l ls 32 not,\
'' u ($1*8):($6+$11)  w l ls 52 not


unset label
set label 'Q_{J,B}, T = 40 ps' at 20,0.25
set out 'fig-meta-iflux-08-1-2.eps'
pl\
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale01.0/set/meta.flux.inte.out' \
u 1:($12+$13+$17+$18) ls 11 w l not, \
'' u 1:($14+$15+$19+$20) ls 31 w l not, \
'' u 1:($16+$21) ls 51 w l not, \
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale08.0/set/meta.flux.inte.out' \
u ($1*8):($12+$13+$17+$18) ls 12 w l not, \
'' u ($1*8):($14+$15+$19+$20)  w l ls 32 not,\
'' u ($1*8):($16+$21)  w l ls 52 not


unset label
set label 'Q_{J,C}, T = 40 ps' at 20,0.25
set out 'fig-meta-iflux-08-1-3.eps'
pl\
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale01.0/set/meta.flux.inte.out' \
u 1:($22+$23) ls 11 w l not, \
'' u 1:($24+$25) ls 31 w l not, \
'' u 1:($26) ls 51 w l not, \
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale08.0/set/meta.flux.inte.out' \
u ($1*8):($22+$23) ls 12 w l not, \
'' u ($1*8):($24+$25)  w l ls 32 not,\
'' u ($1*8):($26)  w l ls 52 not



unset label
set label 'Q_{J,A}, T = 40 ps' at 20,0.25
set out 'fig-meta-iflux-01-check-1-1.eps'
pl\
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale01.0/set/meta.flux.inte.out' \
u 1:($2+$3+$7+$8) ls 11 w l not, \
'' u 1:($4+$5+$9+$10) ls 31 w l not, \
'' u 1:($6+$11) ls 51 w l not, \
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale01.0.dtcomp/set/meta.flux.inte.out' \
u ($1*1):($2+$3+$7+$8) ls 12 w l not, \
'' u ($1*1):($4+$5+$9+$10)  w l ls 32 not,\
'' u ($1*1):($6+$11)  w l ls 52 not


unset label
set label 'Q_{J,B}, T = 40 ps' at 20,0.25
set out 'fig-meta-iflux-01-check-1-2.eps'
pl\
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale01.0/set/meta.flux.inte.out' \
u 1:($12+$13+$17+$18) ls 11 w l not, \
'' u 1:($14+$15+$19+$20) ls 31 w l not, \
'' u 1:($16+$21) ls 51 w l not, \
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale01.0.dtcomp/set/meta.flux.inte.out' \
u ($1*1):($12+$13+$17+$18) ls 12 w l not, \
'' u ($1*1):($14+$15+$19+$20)  w l ls 32 not,\
'' u ($1*1):($16+$21)  w l ls 52 not


unset label
set label 'Q_{J,C}, T = 40 ps' at 20,0.25
set out 'fig-meta-iflux-01-check-1-3.eps'
pl\
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale01.0/set/meta.flux.inte.out' \
u 1:($22+$23) ls 11 w l not, \
'' u 1:($24+$25) ls 31 w l not, \
'' u 1:($26) ls 51 w l not, \
'alanine.charmm.pme.localSD.mode2.040.Ex.01.00.t0200.scale01.0.dtcomp/set/meta.flux.inte.out' \
u ($1*1):($22+$23) ls 12 w l not, \
'' u ($1*1):($24+$25)  w l ls 32 not,\
'' u ($1*1):($26)  w l ls 52 not


