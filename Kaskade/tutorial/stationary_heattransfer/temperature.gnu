set terminal x11
set dgrid3d 50,50 splines
set title "temperature"
set style line 1 lw 0
set pm3d implicit at s
splot "temperature.data" with dots
pause 15
set style fill solid
set pm3d map
splot "temperature.data"
pause 15
