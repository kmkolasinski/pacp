set terminal png size 800,600 enhanced font "Helvetica,14"

a=0.00001
Tc=2.27
gamma=1.75


#set yrange [0:.1]

set o  "IsingTestChi2D.png"
set yl "Chi(T)"
set xl "T"
#set log y
set title "Numerical Chi(T) for 2D Ising chain"
plot "IsingTestChi2D.txt"  u 1:2  w p pt 7 lc -1 t "Chi(Metropolis)", \
     "IsingTestChi2D_w.txt"  u 1:2  w l lc 2 lw 2 t "Chi(Wolff)" #, a*(abs(x-Tc)/Tc)**(-gamma) t "0.001|t|^{-7/4}"

reset

set o  "IsingTestChi2D_N.png"
set yl "Chi(N)"
set xl "N"
set log y
set title "Numerical Chi(N) for 2D Ising chain in T=1 K"
plot "IsingTestChi2D_2.txt"  u 1:2 w p ls 6 lc -1 lw 3 t "Chi(numerical)"
