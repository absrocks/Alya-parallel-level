#touch ../../Sources/kernel/coupli/mod_commdom_driver.f90  
#rm ./Alya.x
#
##
#export PATH=$PATH:/home/bsc21/bsc21704/z2016/REPOSITORY/HPCTOOLKIT/HPCTK/hpctoolkit/build_gcc_papi/Execs/bin
#
#rm Alya.x.hpcstruct
#hpcstruct ./Alya.x 


INIT()
{
  cp configure.in/config_ifort.in config.in
  ./configure -x nastin parall solidz temper alefor chemic turbul
  make metis4
}

PLEPP_INIT()
{
  make libple
  make libplepp
}

PLEPP_COMPILE()
{
  printf -v idx "%02d" $1

  if [ -d plepp${idx} ]; then 
    echo plepp${idx} "EXIST!"   
  else 
    cp -r unix plepp${idx}
    echo plepp${idx} "CREATED!"   
  fi 

  cd plepp${idx}
  make -j4 COMMDOM=$1
  cd .. 
} 

## INIT 
#cd unix 
#INIT 
#PLEPP_INIT
#cd ..  
#
## COMMDOM=X 
PLEPP_COMPILE -2 
#PLEPP_COMPILE  2
PLEPP_COMPILE  3   
PLEPP_COMPILE  4 

 
