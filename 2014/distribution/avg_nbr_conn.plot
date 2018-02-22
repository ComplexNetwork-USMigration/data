set ylabel "Average Neighbor Degree" 
set xlabel "Node degree" 
set terminal postscript eps color "Helvetica" 22 
set output "avg_nbr_conn.ps" 
set key bottom right
set logscale
plot 'out.nbr' using 1:2 title 'Data' with linespoints
