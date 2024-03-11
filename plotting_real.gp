set view map
set pm3d interpolate 4, 4
set dgrid3d 50, 50 gauss 0.6, 0.6

splot "datas_real.dat" using 1:2:3 with pm3d

set terminal qt size 480,480 persist
