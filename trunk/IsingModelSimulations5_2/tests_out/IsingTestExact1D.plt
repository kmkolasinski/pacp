set terminal png size 800,600 enhanced font "Helvetica,14"


set o  "IsingTestExact1D.png"
set yl "Cv(T)/N"
set xl "T"
set yr [0:0.53]
set title "Cv(T) per lattice point for 1D Ising chain. Dots show the analytical solution for given N."
plot "IsingTestExact1DN=2.txt" u 1:2 w l lw 2 t "N=2" , "" u 1:3 every 50 pt 7 lt -1 t "", "IsingTestExact1DN=7.txt"  u 1:2  w l  lw 2 t "N=7" , "" u 1:3 every 50 pt 7 lt -1 t "" ,  "IsingTestExact1DN=22.txt"  u 1:2  w l  lw 2 t "N=22" , "" u 1:3 every 50 pt 7 lt -1 t ""

