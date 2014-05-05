set terminal png size 800,600 enhanced font "Helvetica,14"


set o  "IsingTestM2D.png"
set yl "M(T)/N"
set xl "T"
set title "M(T) per lattice point for 2D Ising NxN Lattice"
plot "IsingTestM2D.txt" u 1:2 w lp  pt 6 lw 3 t "2x2" , ""  u 1:3  w lp  pt 6 lw 3 t "4x4" , ""  u 1:4  w lp  pt 6 lw 3 t "8x8" , ""  u 1:5  w lp  pt 6 lw 3 t "16x16"
