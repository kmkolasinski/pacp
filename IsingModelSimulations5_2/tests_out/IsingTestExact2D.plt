set terminal png size 800,600 enhanced font "Helvetica,14"


set o  "IsingTestExact2D.png"
set yl "CC(T)/N"
set xl "T"
set title "CC(T) per lattice point for 2D Ising NxN Lattice"
plot "IsingTestExact2D2x2.txt" u 1:2 w lp  pt 6 lw 3 t "2x2" , "IsingTestExact2D3x3.txt"  u 1:2  w lp  pt 6 lw 3 t "3x3" , "IsingTestExact2D4x4.txt"  u 1:2  w lp  pt 6 lw 3 t "4x4" , "IsingTestExact2D5x5.txt"  u 1:2  w lp  pt 6 lw 3 t "5x5" 
