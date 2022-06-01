## ERROR: MARENOSTRUM 
##   Catastrophic error: could not set locale "" to allow processing of multibyte characters
## SOLUTION: 
##   export LANG=en_US.utf8
##   export LC_ALL=en_US.utf8 
##
## ERROR: 
##   >> ./code_saturne   
##   >>  python_module is not found  
## SOLUTION:
##   vim .bashrc: export PYTHONPATH=$PYTHONPATH:WHERE_IS 'site-packages/code_saturne'? 
## 
./configure CC=mpicc --enable-static --enable-relocatable --without-cgns --without-scotch --disable-gui --prefix=WHERE_WILL_BE_INSTALLED?  
make 
make install 
