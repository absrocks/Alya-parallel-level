# Gnuplot script file for plotting data f-d curve
# Created by G.Guillamet <gerard.guillamet@bsc.es>
set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set xlabel "True Path (mm)"
set ylabel "x-displacement (mm)"
set y2label "y-displacement (mm)"
set xrange [-1.6:1.6]
set x2range [-1.6:1.6]
set yrange [-0.1:0.1]
set y2range [-0.1:0.1]
set y2tics
plot "alya-solution-path-block-impl.txt" using 1:3 axes x1y1 title 'ALYA UX (Implicit)' with lines lw 2,\
     "alya-solution-path-block-impl.txt" using 1:4 axes x1y2 title 'ALYA UY (Implicit)' with lines lw 2,\
     "alya-solution-path-block.txt" using 1:3 axes x1y1 title 'ALYA UX (Explicit)' with lines lw 2,\
     "alya-solution-path-block.txt" using 1:4 axes x1y2 title 'ALYA UY (Explicit)' with lines lw 2

