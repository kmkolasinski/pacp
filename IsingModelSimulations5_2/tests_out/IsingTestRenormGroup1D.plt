set terminal png size 800,600 enhanced font "Helvetica,14"


set o  "IsingTestRenormGroup1D.png"
set yl "T (N)"
set xl "N"
set yrange [0:40]
set title "Temperature of length of chain for renormalisation for 1D Ising Lattice"
plot "IsingTestRenormGroup1D.txt" u 1:3 w lp  pt 6 lw 3 t ""

reset
set o  "IsingTestRenormGroup1D_b.png"
set yl "1/T'"
set xl "1/T"
set title "Change of interaction constant for 1 step renormalisation of 1D Ising lattice \n in function of initial 1/T"
plot "IsingTestRenormGroup1D_2.txt" u (1/column(1)):(1/(column(5))) w lp  pt 6 lw 3 t "1/T'", x w l ls 0 t "y=x"
