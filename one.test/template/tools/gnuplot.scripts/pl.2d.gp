set term post eps enh color solid size 17cm,8cm
set out 'fig-2d.eps'

#set border 4095 front linetype -1 linewidth 1.000
set view map
set samples 50, 50
set isosamples 50, 50
unset surface
set style data pm3d
set style function pm3d
set ticslevel 0
set xlabel "q [nm]" 
set ylabel "p [nm/ps]"
set xtics 1.0
set ytics 2.0
set pm3d implicit at b
#set size ratio 0.88
set palette rgbformulae 22,13,-31
set size ratio 1.0

# set yrange [-8:8]
set cbrange [-0.1:0.02]
set format cb "%.2f"

#set size 1.4,1.0
set origin 0, 0

set multiplot 

# set ratio 1
set size 0.24,0.6

set origin 0.0,0.45
unset xlabel
set title 't = 0 ps'
# set title 'linear resp., t=00000 ps' font "tt,16"
spl 'indicator.vx.00000.00.out' not

set origin 0.0, 0.0
unset title
set xlabel "q [nm]" 
#set title 'direct, t=00000 ps' font "tt,16"
spl 'indicator.linear.vx.00000.00.out' not

unset ylabel


set origin 0.25, 0.45
unset xlabel
set title 't = 30 ps'
spl 'indicator.vx.00030.00.out' not

set origin 0.5, 0.45
set title 't = 60 ps'
spl 'indicator.vx.00060.00.out' not

set origin 0.75, 0.45
set title 't = 100 ps'
spl 'indicator.vx.00100.00.out' not


unset title
set xlabel "q [nm]"

set origin 0.25, 0.0
spl 'indicator.linear.vx.00030.00.out' not

set origin 0.5, 0.0
spl 'indicator.linear.vx.00060.00.out' not

set origin 0.75, 0.0
spl 'indicator.linear.vx.00100.00.out' not



unset multiplot
