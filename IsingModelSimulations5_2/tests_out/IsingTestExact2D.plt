set terminal png size 800,600 enhanced font "Helvetica,14"


set o  "IsingTestExact2D.png"
set yl "CC(T)/N"
set xl "T"
set title "Cv(T) per lattice point for 2D Ising NxN Lattice"
set arrow from 2.27,0 to 2.27,1 nohead lt -1 lw 2 lc  rgb 'black' 
set yr [0:1.0]
plot "IsingTestCC2DForExact.txt" u 1:2 w p  pt 6 lw 3 lt -1 t "" , ""  u 1:4  w p  pt 6 lw 3 lt -1 t "" , ""  u 1:6  w p  pt 6 lw 3 lt -1 t "" , ""  u 1:8  w p  pt 6 lw 3 lt -1 t "" , "IsingTestExact2D2x2.txt" u 1:2 w l  lw 3 t "2x2" , "IsingTestExact2D3x3.txt"  u 1:2  w l  lw 3 t "3x3" , "IsingTestExact2D4x4.txt"  u 1:2  w l   lw 3 t "4x4" , "IsingTestExact2D5x5.txt"  u 1:2  w l   lw 3 t "5x5" 