set terminal png size 800,600 enhanced font "Helvetica,14"

a=0.001
Tc=2.27
gamma=1.75

#f1(x) = a1*(-(x-Tc1)/Tc1)**(-g1)
#f2(x) = a2*((x-Tc2)/Tc2)**(-g2)

#fit f1(x) "IsingTestChi2D_left.txt" u 1:2 via a1,Tc1,g1
#fit f2(x) "IsingTestChi2D_right.txt" u 1:2 via a2,Tc2,g2


set yrange [0:.6]

set o  "IsingTestChi2D.png"
set yl "Chi(T)"
set xl "T"
#set log y
set title "Numerical Chi(T) for 2D Ising chain"
plot "IsingTestChi2D.txt"  u 1:2:3  w errorbars ls 6 lc -1 lw 3 t "Chi(numerical)", a*(abs(x-Tc)/Tc)**(-gamma) t "0.001|t|^{-7/4}"

reset

set o  "IsingTestChi2D_N.png"
set yl "Chi(N)"
set xl "N"
set log y
set title "Numerical Chi(N) for 2D Ising chain in T=1 K"
plot "IsingTestChi2D_2.txt"  u 1:2 w p ls 6 lc -1 lw 3 t "Chi(numerical)"
