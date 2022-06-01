set grid 
set macros 
set format y '%12.3f' 

#set xlabel '' 
#set xtics format '' 
#set x2tics  
#set arrow from X1,graph 0  to X2,graph 1  nohead 

#set term png 
#set output 'BNDRYLAYER01_2016July08.png' 
#set terminal postscript eps enhanced color font 'Helvetica,10' 
#set output 'BNDRYLAYER01_2016July08.eps' 

T0 = 0.0
#T0 = 0.0653876
#T0 = 0.0308702

 cols =  1
 rows =  2 
 set multiplot layout rows,cols rowsfirst
#set tmargin 0
#set bmargin 0

set xrange[:1.1] 

set ylabel '($2)' 


plot [][895:] \
 'vortex2D_witness2.dat' i 0:0 ev 1 u ( ($1>=0.0)?($1):(1/0) ):($2) t '1' w lp ps -1 pt 1  ,\
"../TOUSE/reference_data01.dat" index 0 u ( $1+T0 ) :2 w lp pt 5 t "ref"


plot [][870:] \
 'vortex2D_witness1.dat' i 0:0 ev 1 u ( ($1>=0.0)?($1):(1/0) ):($2) t '1' w lp ps -1 pt 1  ,\
"../TOUSE/reference_data01.dat" index 1 u ( $1+T0 ) :2 w lp pt 5 t "ref"


 unset multiplot  
 pause -1  
#reread 

#
# 2016July08
# /.statelite/tmpfs/gpfs/projects/bsc21/bsc21704/RUNNER_2016/RUNNER/BENCHMARKS02/BNDRYLAYER01
# python /home/bsc21/bsc21704/z2016/REPOSITORY/TOOLs/Plot_multignuplot.py -F 'CPLNG05_02_0*/vortex2D_witness1.dat' -X  '( ($1>=0.0)?($1):(1/0) )' -Y  '($2)' -Rx ':' -Ry ':' 
#
