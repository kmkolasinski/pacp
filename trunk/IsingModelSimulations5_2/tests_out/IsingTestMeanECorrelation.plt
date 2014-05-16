set terminal png size 800,600 enhanced font "Helvetica,14"


set o  "IsingTestMeanECorrelation.png"
set yl "Mean energy correlation"
set xl "Correlation distance"

set title "Mean Energy Correlation versus distance in fixed temperature"
plot "IsingTestMeanECorrelation.txt" u 1:2 w lp  pt 6 lw 3 notitle
