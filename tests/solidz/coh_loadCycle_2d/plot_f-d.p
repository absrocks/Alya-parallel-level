# Gnuplot script file for plotting data f-d curve
# Created by G.Guillamet <gerard.guillamet@bsc.es>
set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set title "Mode I"
set xlabel "Displacement (mm)"
set ylabel "Force (N)"
set xr [-1e-3:0.03]
set yr [-100:100]
set grid
plot "post-analyticalsol.txt" using 4:7 title 'Analytical' with lines, \
     "post-alyasol.txt" using 3:6 title 'ALYA 2D' lw 1.5 ps 1.5
