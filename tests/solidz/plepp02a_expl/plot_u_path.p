# Gnuplot script file for plotting data f-d curve
# Created by G.Guillamet <gerard.guillamet@bsc.es>

set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set xlabel "Path (m)" font "Times-Roman,9"
set ylabel "UX (m)" font "Times-Roman,9"
set y2label "UY (m)" font "Times-Roman,9"
set title "2D QSI (depth = 0.125)" font "Times-Roman,9"
set xrange [0:1.6]
set x2range [0:1.6]
set yrange [-0.01:0.01]
set y2range [-0.15:0.15]
set y2tics font "Times-Roman, 9"
set xtics font "Times-Roman, 9"
set ytics font "Times-Roman, 9"
set grid ytics
set grid y2tics
set grid xtics
set key top right

#
# define line styles using explicit rgbcolor names
#
set for [i=1:4] linetype i dashtype i
set style line 1 lt 1 lc rgb "red" lw 2
set style line 2 lt 1 lc rgb "black" lw 2
set style line 3 lt 2 lc rgb "red" lw 2
set style line 4 lt 2 lc rgb "black" lw 2
set style line 5 lt 2 lc rgb "green" lw 2
set style line 6 lt 2 lc rgb "blue" lw 2

plot "code_aster-solution-path-block.txt" using 1:3 axes x1y1 title 'Code Aster UX (Imp.)' with lines ls 1,\
     "code_aster-solution-path-block.txt" using 1:4 axes x1y2 title 'Code Aster UY (Imp.)' with lines ls 2,\
     "alya-solution-path-block_old.txt" using 1:3 axes x1y1 title 'ALYA UX (Exp. old tol=1e-3)' with lines ls 3,\
     "alya-solution-path-block_old.txt" using 1:4 axes x1y2 title 'ALYA UY (Exp. old tol=1e-3)' with lines ls 4,\
     "alya-solution-path-block.txt" using 1:3 axes x1y1 title 'ALYA UX (Exp. new tol=1e-4 )' with lines ls 5,\
     "alya-solution-path-block.txt" using 1:4 axes x1y2 title 'ALYA UY (Exp. new tol=1e-4 )' with lines ls 6
