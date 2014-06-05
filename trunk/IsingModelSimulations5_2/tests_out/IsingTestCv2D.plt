set terminal png size 800,600 enhanced font "Helvetica,14"


set o  "IsingTestCv2D.png"
set yl "Cv(T)/N"
set xl "T"
set title "Cv(T) per lattice point for 2D Ising NxN Lattice"
plot "IsingTestCv2D.txt" u 1:2 w lp  pt 6 lw 3 t "2x2" , ""  u 1:4  w lp  pt 6 lw 3 t "4x4" , ""  u 1:6  w lp  pt 6 lw 3 t "8x8" , ""  u 1:8  w lp  pt 6 lw 3 t "16x16"
