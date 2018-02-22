set ylabel "Rich Club Connectivity" 
set xlabel "Normalized rank r/N " 
set terminal postscript eps color "Helvetica" 22 
set output "rich_club.ps" 
set key bottom right
set logscale
plot 'rich_club' using 1:2 title 'Data' with linespoints
