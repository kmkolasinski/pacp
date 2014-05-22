set terminal png size 800,600 enhanced font "Helvetica,14"


set o  "IsingTestError2D.png"
set yl "Error(prod. time)"
set xl "prod. time"
set title "Error(prod. time) of Chi for 2D Ising Lattice"
plot "IsingTestError2D.txt" u 1:2 w lp  pt 6 lw 3 t ""

#set o  "IsingTestError2D_Chi.png"
#set yl "Chi(prod. time)"
#set xl "prod. time"
#set title "Chi(prod. time) with errorbars for 2D Ising Lattice "
#plot "IsingTestError2D.txt" u 1:3:($3-$2):($3+$2) w lp with yerrorbars pt 6 lw 3 t ""
