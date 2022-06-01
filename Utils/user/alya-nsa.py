GNUPLOT = """
set macro 
FILE_NAME = \"%s\"
set multiplot layout 2,2 rowsfirst
plot FILE_NAME u 1:4 t "Density min." w lines ls 1, "" u 1:5 t "Density max." w l ls 2 
plot FILE_NAME u 1:6 t "Pressure min." w lines ls 1,"" u 1:7 t "Pressure max." w lines ls 2 
plot FILE_NAME u 1:8 t "Temper. min." w lines ls 1,"" u 1:9 t "Temper. max." w lines ls 2
plot FILE_NAME u 1:10 t "Mach min." w lines ls 1, "" u 1:11 t "Mach max." w lines ls 2
unset multiplot

pause -1 
"""

import sys
import os 

GNUPLOT = GNUPLOT % sys.argv[1]
FILE = open("nsa_phys.gp", "w")
FILE.write(GNUPLOT)
FILE.close() 

os.system("gnuplot nsa_phys.gp")

