set terminal png size 800,600 enhanced font "Helvetica,14"


set o  "IsingTestRenormGroup1D.png"
set yl "T (N)"
set xl "N"
set yrange [0:50]
set title "Temperature of length of chain for renormalisation for 1D Ising Lattice"
plot "IsingTestRenormGroup1D.txt" u 1:3 w p  pt 6 lw 3 t ""


