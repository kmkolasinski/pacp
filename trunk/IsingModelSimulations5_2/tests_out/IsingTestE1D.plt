set terminal png size 800,600 enhanced font "Helvetica,14"

set output  "IsingTestE1D.png"
set ylabel "Energy"
set xlabel "Temperature"
set title "Energy function versus temperature"

plot "IsingTestE1D.txt" u 1:2 w l lw 2 notitle
