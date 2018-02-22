set ylabel "Node Degree Distribution" 
set xlabel "Node degree" 
set terminal postscript eps color "Helvetica" 22 
set output "degree_distr.ps" 
set key bottom right
set logscale
plot 'out.deg' using 1:2 title 'Data' with linespoints
