#
#  plot all slices
#
reset
#set xrange [0:0.5];
#set yrange [0:500];
set ylabel 'number of particules';
set xlabel 'Time(s)';
set title "number of particules crossing slices over the time"
#set nokey
#set yzeroaxis lt -1
plot "slice1.dat" u 1:2 with lines ti "slice1",\
     "slice2.dat" u 1:2 with lines ti "slice2",\
     "slice3.dat" u 1:2 with lines ti "slice3",\
     "slice4.dat" u 1:2 with lines ti "slice4",\
     "slice5.dat" u 1:2 with lines ti "slice5",\
     "slice6.dat" u 1:2 with lines ti "slice6",\
     "slice7.dat" u 1:2 with lines ti "slice7",\
     "slice8.dat" u 1:2 with lines ti "slice8"

set term post eps enhanced color solid defaultplex 'Helvetica' 24
set output 'allslice.eps'
replot
set output