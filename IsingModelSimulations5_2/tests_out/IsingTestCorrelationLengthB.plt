#gnuplot
set terminal png;
set output "IsingTestCLvsMCS.png";
set title "Correlation Length in function of mean cluster size for T[0.5;1.5]. N=1000";
set ylabel "Mean Cluster Size";
set xlabel "Correlation length";
fit a*x+b  "IsingTestCLvsMCS.txt" via a,b;
set label "a= %f",a at 6,45;
set label "b= %f",b at 6,40;
plot "IsingTestCLvsMCS.txt" using 1:2 pt 13 notitle, a*x+b notitle;
