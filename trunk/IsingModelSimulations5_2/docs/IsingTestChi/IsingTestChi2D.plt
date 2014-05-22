set terminal png size 800,600 enhanced font "Helvetica,14"

a=0.000007
Tc=2.27
gamma=1.75

#f1(x) = a1*(-(x-Tc1)/Tc1)**(-g1)
#f2(x) = a2*((x-Tc2)/Tc2)**(-g2)

#fit f1(x) "IsingTestChi2D_left.txt" u 1:2 via a1,Tc1,g1
#fit f2(x) "IsingTestChi2D_right.txt" u 1:2 via a2,Tc2,g2


set yrange [0:.025]

set o  "IsingTestChi2D.png"
set yl "Chi(T)"
set xl "T"
#set log y
set title "Numerical Chi(T) for 2D Ising chain"
plot "300x300.txt"  u 3:5  w p ls 6 lc -1 lw 3 t "Chi(numerical)" , a*((x-Tc)/Tc)**(-gamma) t "0.000007|t|^{-7/4}"

