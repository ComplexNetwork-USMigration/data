set ylabel "CCDF" 
set xlabel "Node degree" 
set terminal postscript eps color "Helvetica" 22 
set output "ccdf_deg.ps" 
set key bottom right
set logscale
plot 'deg.ccdf' using 2:1 with linespoints
