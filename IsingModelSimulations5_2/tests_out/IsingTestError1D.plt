set terminal png size 800,600 enhanced font "Helvetica,14"

f(x)=a/sqrt(x)
a=100
fit f(x) 'IsingTestError1D.txt' using 1:2 via a
set o  "IsingTestError1D.png"
set yl "Error(prod. time)"
set xl "prod. time"
set title "Error(prod. time) of Chi for 1D Ising Lattice"
plot "IsingTestError1D.txt" u 1:2 w lp  pt 6 lw 3 t "error", a/sqrt(x) lw 2 ti "fit a/sqrt(prod. time)"


set o  "IsingTestError1D_T.png"
set yl "Error(T)"
set xl "T"
set xrange [1.0:3]
set title "Error of Chi and difference between numerical \n and analytical Chi in function of T for 1D Ising Lattice"
plot "IsingTestError1D_T.txt" u 1:2 w lp  pt 6 lw 3 t "error", "" u 1:(abs(column(3)-column(4))) w lp pt 6 lw 3 t "|Chi_{an} - Chi_{num}|"
