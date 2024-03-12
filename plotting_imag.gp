set view map
set pm3d interpolate 7, 7
set dgrid3d 100, 100

splot "datas_imag.dat" using 1:2:3 with pm3d

set terminal qt size 480,480 persist
