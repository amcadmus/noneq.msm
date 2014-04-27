#set border 4095 front linetype -1 linewidth 1.000
set view map
set samples 50, 50
set isosamples 50, 50
unset surface
set style data pm3d
set style function pm3d
set ticslevel 0
set xlabel "{/Symbol f} [deg.]" 
set ylabel "{/Symbol y} [deg.]"
set xtics 60
set ytics 60
set mxtics 6
set mytics 6

set pm3d implicit at b
# set palette rgbformulae 22,13,-31
set palette gray negative
set size ratio 1.0

set cbrange [0:2e-4]
set cbrange [0:2e-6]
set format cb "%.1e"
set xrange [-180:180]
set yrange [-180:180]

set term post eps enh color solid font 32 size 24cm,16cm 
set out 'fig-equi-ram.eps'
set grid
# set title 't=1000.0 ps'
spl 'equi.dist.out' not
