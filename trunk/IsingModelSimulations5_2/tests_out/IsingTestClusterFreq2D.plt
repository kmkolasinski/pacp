set terminal png size 800,600 enhanced font "Helvetica,14"

set o  "IsingTestClusterFreq2D.png"
set yl "Y"
set xl "X"
set title "Cluster size versus temperature for N=100 NxN 2D lattice"
set xrange [0:10];
set yrange [0:10];
set palette gray;
plot "IsingTestClusterFreq2D.txt" matrix with image notitle
