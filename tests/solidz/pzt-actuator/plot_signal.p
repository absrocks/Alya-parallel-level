#
# Set overall margins for the combined set of plots and size them
# to generate a requested inter-plot spacing
#
if (!exists("MP_LEFT"))   MP_LEFT = .1
if (!exists("MP_RIGHT"))  MP_RIGHT = .95
if (!exists("MP_BOTTOM")) MP_BOTTOM = .1
if (!exists("MP_TOP"))    MP_TOP = .9
if (!exists("MP_GAP"))    MP_GAP = 0.05

set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
set format y "%.1f"
set format x "%.1f"
set key box opaque
set yrange [-1:1]
set xrange [0:0.2]
set grid

# Retrieve statistical properties
max(a,b)=a>b?a:b

set xlabel 'Time (ms)'
set ylabel 'Normalized signal'
a=-1e+10
c=-1e+10
plot 'input-signal.txt' using ($1*1e3):(a=max(a,$2),0/0) notitle, '' using ($1*1e3):($2/a) title 'Input Signal' with lines, \
     'alya-sets-1.dat' using ($1*1e3):(c=max(c,$3),0/0) notitle, '' using ($1*1e3):($3/c) title 'Post strains' with lines, \

