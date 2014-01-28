set view map
set contour
unset surface
unset ztics
set cntrparam levels discrete 1e-6, 10, 20, 30, 41, 50, 60, 67, 80, 90, 100
set size ratio -1
set grid
set xlabel "In plane distance from Earth (Mkm)"
set ylabel "Out of plane distance from Earth (Mkm)"
splot 'sun-earth-sail-sigma.dat' with lines title "Sail loading (g/m^2)"
