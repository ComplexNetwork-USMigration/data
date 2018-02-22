set ylabel "CCDF" 
set xlabel "Average Neighbor degree" 
set terminal postscript eps color "Helvetica" 22 
set output "ccdf_avg_nbr.ps" 
set key bottom right
set logscale
plot 'nbr.ccdf' using 1:2 with linespoints
