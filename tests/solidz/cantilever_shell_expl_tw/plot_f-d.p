# Gnuplot script file for plotting data f-d curve
# Created by G.Guillamet <gerard.guillamet@bsc.es>
set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set title "Cantilever Beam"
set xlabel "Tip Displacement (m)"
set ylabel "RF(N)"
set xr [0:7.0]
set yr [0:90.0]
set key right bottom
plot "alya-sets.out" using 3:4 title 'Alya' with lines
