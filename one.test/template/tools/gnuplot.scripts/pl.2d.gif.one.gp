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
set size ratio 0.90

set cbrange [-0.1:0.02]
set format cb "%.2f"

set term gif
set out 'fig-2d-linear-150.gif'
spl 'indicator.linear.vx.00150.00.out' not
