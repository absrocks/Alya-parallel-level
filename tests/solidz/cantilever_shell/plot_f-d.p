# Gnuplot script file for plotting data f-d curve
# Created by G.Guillamet <gerard.guillamet@bsc.es>
set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set title "Cantilever Beam"
set xlabel "Tip Displacement (m)"
set ylabel "q/qmax"
set xr [0:7.0]
set yr [0:1.0]
set key right bottom
plot "abaqus-cantilver_y-direction.txt" using 1:2 title 'Reinoso Abaqus y-direction' with lines, \
     "abaqus-cantilver_z-direction.txt" using 1:2 title 'Reinoso Abaqus z-direction' with lines, \
     "alya-solution.txt" using 2:4 title 'ALYA y-direction' with points, \
     "alya-solution.txt" using (1*$3):4 title 'ALYA z-direction' with points
