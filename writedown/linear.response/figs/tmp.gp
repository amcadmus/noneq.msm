set term gif size 960,480
set out 'fig-2d-tog-00030.00.gif'

#set border 4095 front linetype -1 linewidth 1.000
set view map
set samples 50, 50
set isosamples 50, 50
unset surface
set style data pm3d
set style function pm3d
set ticslevel 0
set xlabel "x" 
set ylabel "v" 
set pm3d implicit at b
#set size ratio 0.88
set palette rgbformulae 22,13,-31
set size ratio 1.0

# set yrange [-8:8]
set cbrange [-0.04:0.1]
set format cb "%.2f"

#set size 1.4,1.0
set origin 0, 0

set multiplot 

set size 0.5,1.0

set origin 0.0, 0.0
set title 'direct, t=00030.00 ps' font "tt,16"
spl 'indicator.vx.00030.00.out' u 1:2:(-$3) not

set origin 0.5,0.0
set title 'linear resp., t=00030.00 ps' font "tt,16"
spl 'indicator.linear.corr.vx.00030.00.out' u 1:2:(-$3) not

unset multiplot
