############################################################################################
#depo_time
############################################################################################
set term postscript enhanced color "Helvetica" 20
set output "depo_time.eps"
set xlabel "Time (second)"
set ylabel "Number of particules deposited"
set title "Particules deposited above Z = 0.04388 over the time "
plot "depo_time.dat" u 1:2 with lines ti "Z >= 0.04388"

