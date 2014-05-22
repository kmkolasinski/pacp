set terminal png size 800,600 enhanced font "Helvetica,14"


set o  "IsingTestError1D.png"
set yl "Error(prod. time)"
set xl "prod. time"
set title "Error(prod. time) of Chi for 1D Ising Lattice"
plot "IsingTestError1D.txt" u 1:2 w lp  pt 6 lw 3 t ""
