# Gnuplot script file for plotting data f-d curve
# Created by G.Guillamet <gerard.guillamet@bsc.es>
set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set xlabel "Time (s)"
set ylabel "Pressure (MPa)"
set xr [0:0.005]
set yr [0:1000.0]
set key right bottom
set grid
set tics out
set size square

plot "alya-sets.out" using 1:(abs($3)) title 'PRESS' with lines

