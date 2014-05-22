set terminal png size 800,600 enhanced font "Helvetica,14"


set o  "IsingTestChi2D.png"
set yl "Chi(T)"
set xl "T"
set log y
set title "Numerical and analitycal Chi(T) for 2D Ising chain"
#plot "IsingTestChi2D.txt" u 1:2 w l lw 3 t "Chi(analytical)" , ""  u 1:3 w p ls 6 lc -1 lw 3 t "Chi(numerical)"
plot "IsingTestChi2D.txt"  u 1:2 w p ls 6 lc -1 lw 3 t "Chi(numerical)"


set o  "IsingTestChi2D_N.png"
set yl "Chi(N)"
set xl "N"
set log y
set title "Numerical and analitycal Chi(N) for 2D Ising chain in T=1 K"
plot "IsingTestChi2D_2.txt"  u 1:2 w p ls 6 lc -1 lw 3 t "Chi(numerical)"
