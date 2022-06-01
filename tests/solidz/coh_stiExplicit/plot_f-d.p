# Gnuplot script file for plotting data
# Created by G.Guillamet <gerard.guillamet@bsc.es>
# This file is called  plot_f-d.p
set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set title "Mode II"
set xlabel "Displacement (mm)"
set ylabel "Force (N)"
set xr [0:0.04]
set yr [0:120]
plot "analytical.txt" using 2:(1*$5) title 'Analytical' with lines, \
     "alya-solution.txt" using 2:(1*$5) title 'ALYA' with lines           
