set terminal png size 800,600 enhanced font "Helvetica,14"


set o  "IsingTestChi1D.png"
set yl "Chi(T)"
set xl "T"
set log y
set title "Numerical and analitycal Chi(T) for 1D Ising chain"
plot "IsingTestChi1D.txt" u 1:2 w l lw 3 t "Chi(analytical)" , ""  u 1:3:5:6 w yerrorbars ls 6 lc -1 lw 3 t "Chi(numerical)"
