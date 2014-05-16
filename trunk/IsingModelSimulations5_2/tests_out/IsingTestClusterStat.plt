set terminal png size 800,600 enhanced font "Helvetica,14"

set out  "IsingTestClusterStat.png"
set yl "Cluster size"
set xl "Number of clusters"
set title "Cluster distibution comparition for Wolf and Metopolis algorithmes in fixed temperature"
plot "IsingTestClusterStat.txt" u 1:2 w lp  pt 6 lw 3 t "wolf" , ""  u 1:3  w lp  pt 6 lw 3 t "metropolis" 
