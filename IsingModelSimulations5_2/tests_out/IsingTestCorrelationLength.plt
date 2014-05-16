#gnuplot

fit a1*exp(-x/cl1) "IsingTestCorrelationLength.txt" using 1:2 via a1,cl1;
fit a2*exp(-x/cl2) "IsingTestCorrelationLength.txt" using 1:3 via a2,cl2;
fit a3*exp(-x/cl3) "IsingTestCorrelationLength.txt" using 1:4 via a3,cl3;
fit a4*exp(-x/cl4) "IsingTestCorrelationLength.txt" using 1:5 via a4,cl4;
fit a5*exp(-x/cl5) "IsingTestCorrelationLength.txt" using 1:6 via a5,cl5;
fit a6*exp(-x/cl6) "IsingTestCorrelationLength.txt" using 1:7 via a6,cl6;
fit a7*exp(-x/cl7) "IsingTestCorrelationLength.txt" using 1:8 via a7,cl7;
fit a8*exp(-x/cl8) "IsingTestCorrelationLength.txt" using 1:9 via a8,cl8;
fit a9*exp(-x/cl9) "IsingTestCorrelationLength.txt" using 1:10 via a9,cl9;
fit a10*exp(-x/cl10) "IsingTestCorrelationLength.txt" using 1:11 via a10,cl10;

# set title "Tytuł wykresu";
set xlabel "corr";
set ylabel "distance";

set print "cl.txt";
print cl1;
print cl2;
print cl3;
print cl4;
print cl5;
print cl6;
print cl7;
print cl8;
print cl9;
print cl10;

set terminal png;
