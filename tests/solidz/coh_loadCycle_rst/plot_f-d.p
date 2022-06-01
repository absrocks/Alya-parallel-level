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
set xr [0:0.03]
set yr [-10.0:100]
plot "analytical.txt" using 1:2 title 'Analytical' with lines, \
     "alya-solution1.txt" using 2:3 title '1st part' lw 1.5 ps 1.5, \
     "alya-solution2.txt" using 2:3 title '2nd part with Restart' pt 7 lc 7 lw 1.5 ps 1.5, \
