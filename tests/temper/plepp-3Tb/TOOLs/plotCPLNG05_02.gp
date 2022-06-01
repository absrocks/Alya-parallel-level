 # created by J. MIGUEL ZAVALA-AKE
  set grid 
  set xlabel "t"
  set ylabel "TEMPE "
  #fit a*x+b 'xxx.dat' u 1:2 via a, b 
  set macro 
  Pref = -1
 
  COL = 2
  set multiplot layout 2,1 rowsfirst
  #set multiplot layout 1,1 rowsfirst

  set grid 
 #set yrange [870:900]
  set xrange [0:1.2]

T0 = 0.0 
#T0 = 0.0653876
#T0 = 0.0308702
 
  plot [:][895:]\
"CPLNG05_02/vortex2D_witness2.dat" u ( ($1>=0.0)?($1):(1/0) ):COL w lp t "02" ,\
"CPLNG05_02_03/vortex2D_witness2.dat" u ( ($1>=0.0)?($1):(1/0) ):COL w lp t "03" ,\
"CPLNG05_02_04/vortex2D_witness2.dat" u ( ($1>=0.0)?($1):(1/0) ):COL w lp t "04" ,\
"CPLNG05_02_05/vortex2D_witness2.dat" u ( ($1>=0.0)?($1):(1/0) ):COL w lp t "05" ,\
"CPLNG05_02_06/vortex2D_witness2.dat" u ( ($1>=0.0)?($1):(1/0) ):COL w lp t "06" ,\
"CPLNG05_02_07/vortex2D_witness2.dat" u ( ($1>=0.0)?($1):(1/0) ):COL w lp t "07" ,\
"reference_data01.dat" index 0 u ( $1+T0 ) :2 w lp pt 5 t "ref"


  plot [:][870:]\
"CPLNG05_02/vortex2D_witness1.dat" u ( ($1>=0.0)?($1):(1/0) ):COL w lp t "02" ,\
"CPLNG05_02_03/vortex2D_witness1.dat" u ( ($1>=0.0)?($1):(1/0) ):COL w lp t "03" ,\
"CPLNG05_02_04/vortex2D_witness1.dat" u ( ($1>=0.0)?($1):(1/0) ):COL w lp t "04" ,\
"CPLNG05_02_05/vortex2D_witness1.dat" u ( ($1>=0.0)?($1):(1/0) ):COL w lp t "05" ,\
"CPLNG05_02_06/vortex2D_witness1.dat" u ( ($1>=0.0)?($1):(1/0) ):COL w lp t "06" ,\
"CPLNG05_02_07/vortex2D_witness1.dat" u ( ($1>=0.0)?($1):(1/0) ):COL w lp t "07" ,\
"reference_data01.dat" index 1 u ( $1+T0 ) :2 w lp pt 5 t "ref"


# 0.0653876 
# 0.0308702 



  unset multiplot
  pause -1 
  
