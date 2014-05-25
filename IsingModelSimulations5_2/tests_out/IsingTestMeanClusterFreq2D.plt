set terminal png size 800,600 enhanced font "Helvetica,14"

set o  "IsingTestMeanClusterFreq2D.png"
set yl "Cluster size"
set xl "Temperature"
set title "Mean cluster size versus temperature for N=100 NxN 2D lattice"
plot "IsingTestMeanClusterFreq2D.txt" u 1:2 w lp  pt 6 lw 3 notitle
