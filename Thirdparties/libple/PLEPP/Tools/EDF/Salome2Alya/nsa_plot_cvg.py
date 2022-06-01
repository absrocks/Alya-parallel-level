GNUPLOT = """
set macro 
FILE_NAME = \"%s\"
set multiplot layout 2,2 rowsfirst
plot FILE_NAME u 4:5 t "Total Res" w lines ls 1 
plot FILE_NAME u 4:6 t "Mom. Res" w lines ls 1 
plot FILE_NAME u 4:7 t "Cont. Res" w lines ls 1 
plot FILE_NAME u 1:8 t "Ener. Res" w lines ls 1 
unset multiplot

pause -1 
"""

import sys
import os 

GNUPLOT = GNUPLOT % sys.argv[1]
FILE = open("nsa_props.gp", "w")
FILE.write(GNUPLOT)
FILE.close() 

os.system("gnuplot nsa_props.gp")

