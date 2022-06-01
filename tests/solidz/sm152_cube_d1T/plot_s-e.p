# Gnuplot script file for plotting data f-d curve
# Created by G.Guillamet <gerard.guillamet@bsc.es>

#set terminal aqua dashed enhanced
#set terminal x11 dashed nopersist enhanced font "arial,15"
set terminal wxt dashed nopersist enhanced font "arial,15"

set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set title "Cube - d1T"
set xlabel "Strain (-)"
set ylabel "Stress (MPa)"
set xr [0:0.25]
set yr [0:2500]
set key right top

#
# define line styles using explicit rgbcolor names
#
set for [i=1:3] linetype i dashtype i
set style line 1 lt 1 lc rgb "red" lw 1 pt 1
set style line 2 lt 1 lc rgb "black" lw 1 pt 1
set style line 3 lt 1 lc rgb "green" lw 1 pt 2

plot "exp-XT.txt" using 1:2 title 'Exp. XT' with lines ls 1, \
     "alya-solution-200steps.txt" using ($2):(1*$3) title 'ALYA 200 steps' with linespoints ls 2, \
     "alya-solution.txt" using ($2):(1*$3) title 'ALYA' with linespoints ls 3
