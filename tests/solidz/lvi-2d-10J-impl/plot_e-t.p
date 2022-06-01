# Gnuplot script file for plotting data e-t curve
# Created by G.Guillamet <gerard.guillamet@bsc.es>

#set terminal aqua dashed enhanced
#set terminal x11 dashed nopersist enhanced font "arial,15"
set terminal wxt dashed nopersist enhanced font "arial,15"

set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set title "Impact E-t"
set xlabel "Time (ms)"
set ylabel "Energy (Joules)"
#set xr [0:3.0]
#set yr [0:2]
set key right top
set grid

#
# define line styles using explicit rgbcolor names
#
set for [i=1:3] linetype i dashtype i
set style line 1 lt 1 lc rgb "black" lw 2
set style line 2 lt 1 lc rgb "red" lw 2
set style line 3 lt 1 lc rgb "blue" lw 2
set style line 4 lt 1 lc rgb "green" lw 1
set style line 5 lt 1 lc rgb "orange" lw 1
set style line 6 lt 2 lc rgb "grey" lw 1
set style line 7 lt 2 lc rgb "green" lw 1

plot "alya-sets-seq-ref.txt" using ($1*1e+3):(10-0.5*2e-6*$3*$3) title 'ALYA Ref.' with lines ls 1, \
     "alya-sets.out" using ($1*1e+3):(10-0.5*2e-6*$3*$3) title 'ALYA current' with lines ls 2
