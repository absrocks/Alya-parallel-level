set style data histogram
set style histogram cluster gap 1
set style fill solid border -1
#set boxwidth 0.9
#set xtic rotate by -45 scale 0
#set bmargin 10 
set multiplot layout 1,1  
set key outside
set grid 

set title "DIRIC"
plot './dirichl_000004_000.dat' \
   using 2:xtic(1) t 'CONBLK1' ,\
'' using 4         t 'CONCOU1' ,\
'' using 5         t 'ENDZON1' ,\
'' using 6         t 'ENDSTE1' ,\
'' using 7         t 'TIMSTE1' ,\
'' using 9         t 'BEGZON1' ,\
'' using 11         t 'BEGSTE1' ,\
'' using 12         t 'DOITER1' ,\
'./neumann_000003_000.dat' \
   using 2:xtic(1) t 'CONBLK2' ,\
'' using 4         t 'CONCOU2' ,\
'' using 5         t 'ENDZON2' ,\
'' using 6         t 'ENDSTE2' ,\
'' using 7         t 'TIMSTE2' ,\
'' using 9         t 'BEGZON2' ,\
'' using 11         t 'BEGSTE2' ,\
'' using 12         t 'DOITER2'


unset multiplot

pause -1 

