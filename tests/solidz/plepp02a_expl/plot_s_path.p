# Gnuplot script file for plotting data f-d curve
# Created by G.Guillamet <gerard.guillamet@bsc.es>

set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set xlabel "Path (m)" font "Times-Roman,9"
set ylabel "Stress (Pa)" font "Times-Roman,9"
set title "2D QSI (depth = 0.125)" font "Times-Roman,9"
set xrange [0:1.6]
set yrange [-1.5e+8:5e+7]
set xtics font "Times-Roman, 9"
set ytics font "Times-Roman, 9"
set grid ytics
set grid xtics
set key right bottom

#
# define line styles using explicit rgbcolor names
#
set for [i=1:4] linetype i dashtype i
set style line 1 lt 1 lc rgb "red" lw 2
set style line 2 lt 1 lc rgb "orange" lw 2
set style line 3 lt 1 lc rgb "grey" lw 2
set style line 4 lt 2 lc rgb "red" lw 2
set style line 5 lt 2 lc rgb "orange" lw 2
set style line 6 lt 2 lc rgb "grey" lw 2

plot "code_aster-solution-path-block.txt" using 1:5 title 'Code Aster SXX (Imp.)' with lines ls 1,\
     "code_aster-solution-path-block.txt" using 1:6 title 'Code Aster SYY (Imp.)' with lines ls 2,\
     "code_aster-solution-path-block.txt" using 1:7 title 'Code Aster SXY (Imp.)' with lines ls 3,\
     "alya-solution-path-block.txt" using 1:5 title 'ALYA SXX (Exp.)' with lines ls 4,\
     "alya-solution-path-block.txt" using 1:6 title 'ALYA SYY (Exp.)' with lines ls 5,\
     "alya-solution-path-block.txt" using 1:7 title 'ALYA SXY (Exp.)' with lines ls 6

