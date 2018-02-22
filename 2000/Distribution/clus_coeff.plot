set ylabel "Cluster Coefficient" 
set xlabel "Node degree" 
set terminal postscript eps color "Helvetica" 22 
set output "clus_coeff.ps" 
set key bottom right
set logscale
plot 'out.clus' using 1:2 title 'Observed Value' with linespoints
