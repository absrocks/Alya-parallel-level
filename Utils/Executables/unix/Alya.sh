#!/bin/csh -f
# -------------------------------------------------------------------
#
#                          A  L  Y  A
#
# -------------------------------------------------------------------

# -------------------------------------------------------------------
#   Identify arguments
# -------------------------------------------------------------------

set version=release
if ( $1 == "-d" ) goto debug
if ( $1 == "-r" ) goto release
if ( $1 == "-g" ) goto debug
if ( $1 == "-x" ) goto release
goto readarguments

debug:
set version=debug
shift
goto readarguments

release:
set version=release
shift
goto readarguments

readarguments:
set name=$1
set directory=$2
set ProblemType=$3
set Gidexe=$4

if ( $1 == "" )       goto help
if ( $1 == "" )       goto help
if ( $1 == "--help" ) goto help
if ( $1 == "-help" )  goto help
if ( $1 == "-h" )     goto help
if ( $1 == "/?" )     goto help
if ( $2 == "" )       goto emptydir

# -------------------------------------------------------------------
#   Script was called from GiD: creates the directory 
# -------------------------------------------------------------------

set gid=1
set input=$directory/data
set output=$directory/results
if ( != -e $input )  mkdir -p $input
if ( != -e $output ) mkdir -p $output

# -------------------------------------------------------------------
#   Move existing files
# -------------------------------------------------------------------

if ( -e $input/$name.dat      ) rm $input/$name.dat    
if ( -e $input/$name.ker.dat  ) rm $input/$name.ker.dat    
if ( -e $input/$name.dom.dat  ) rm $input/$name.dom.dat
if ( -e $input/$name.set.dat  ) rm $input/$name.set.dat
if ( -e $input/$name.geo.dat  ) rm $input/$name.geo.dat
if ( -e $input/$name.nsi.dat  ) rm $input/$name.nsi.dat
if ( -e $input/$name.tem.dat  ) rm $input/$name.tem.dat
if ( -e $input/$name.cdr.dat  ) rm $input/$name.cdr.dat
if ( -e $input/$name.tur.dat  ) rm $input/$name.tur.dat
if ( -e $input/$name.exm.dat  ) rm $input/$name.exm.dat
if ( -e $input/$name.nsa.dat  ) rm $input/$name.nsa.dat
if ( -e $input/$name.ale.dat  ) rm $input/$name.ale.dat
if ( -e $input/$name.sld.dat  ) rm $input/$name.sld.dat
if ( -e $input/$name.got.dat  ) rm $input/$name.got.dat
if ( -e $input/$name.wav.dat  ) rm $input/$name.wav.dat
if ( -e $input/$name.lev.dat  ) rm $input/$name.lev.dat

# -------------------------------------------------------------------
#   Rename the files of the Master
# -------------------------------------------------------------------

if ( -e $directory/$name.dat    )  mv $directory/$name.dat     $input/$name.dat
if ( -e $directory/$name-1.dat  )  mv $directory/$name-1.dat   $input/$name.ker.dat
if ( -e $directory/$name-2.dat  )  mv $directory/$name-2.dat   $input/$name.dom.dat
if ( -e $directory/$name-3.dat  )  mv $directory/$name-3.dat   $input/$name.geo.dat
if ( -e $directory/$name-4.dat  )  mv $directory/$name-4.dat   $input/$name.set.dat
if ( -e $directory/$name-5.dat  )  mv $directory/$name-5.dat   $input/$name.fix.dat
if ( -e $directory/$name-6.dat  )  mv $directory/$name-6.dat   $input/$name.fie.dat

# -------------------------------------------------------------------
#   Rename the files of the modules
# -------------------------------------------------------------------

if ( -e $directory/$name-7.dat  ) mv $directory/$name-7.dat   $input/$name.nsi.dat
if ( -e $directory/$name-8.dat  ) mv $directory/$name-8.dat   $input/$name.tem.dat
if ( -e $directory/$name-9.dat  ) mv $directory/$name-9.dat   $input/$name.cdr.dat
if ( -e $directory/$name-10.dat ) mv $directory/$name-10.dat  $input/$name.tur.dat
if ( -e $directory/$name-11.dat ) mv $directory/$name-11.dat  $input/$name.exm.dat
if ( -e $directory/$name-12.dat ) mv $directory/$name-12.dat  $input/$name.nsa.dat
if ( -e $directory/$name-13.dat ) mv $directory/$name-13.dat  $input/$name.ale.dat
if ( -e $directory/$name-14.dat ) mv $directory/$name-14.dat  $input/$name.sld.dat
if ( -e $directory/$name-15.dat ) mv $directory/$name-15.dat  $input/$name.wav.dat
if ( -e $directory/$name-16.dat ) mv $directory/$name-16.dat  $input/$name.lev.dat

goto continue

# -------------------------------------------------------------------
#   Script was called from DOS: create result directory
# -------------------------------------------------------------------

emptydir:
set directory=.
set gid=0
set input=$directory/data
set output=$directory/results
if( != -e $output ) mkdir -p $output

continue:
# -------------------------------------------------------------------
#   Define environment variables
# -------------------------------------------------------------------

#
# Master
#
setenv ALYA_NAME $name
setenv FOR011  $input/$name.dat           
setenv FOR012 $output/$name.log       
setenv FOR013 $output/$name.mem
setenv FOR014 $output/$name.cvg
setenv FOR015 $output/$name
setenv FOR016 $output/$name.liv
setenv FOR017  $input/$name.rst    
setenv FOR018 $output/$name-latex.tex    
setenv FOR019 $output/$name-latex.plt 
setenv FOR021  $input/$name.dom.dat 
setenv FOR023 $output/$name-DD.post.msh
setenv FOR024 $output/$name.dom.bin
setenv FOR025 $output/$name-elsest.log 
setenv FOR026 $output/$name-elsest.post.msh
setenv FOR027 $output/$name-elsest.post.res
setenv FOR028 $output/$name-system.log    
setenv FOR029 $output/$name-IB.res
setenv FOR030 $output/$name.ensi.case
setenv FOR034 $output/$name
setenv FOR039 $output/$name
setenv FOR045  $input/$name-IB.rst    
setenv FOR046 $output/$name    
setenv FOR047 $output/$name
setenv FOR048 $output/$name    
setenv FOR049 $output/$name

if ( $gid == "1" ) setenv FOR015 $directory/$name
if ( $gid == "1" ) setenv FOR034 $directory/$name
if ( $gid == "1" ) setenv FOR039 $directory/$name
#
# NASTIN: module 1
#
setenv FOR101  $input/$name.nsi.dat        
setenv FOR102 $output/$name.nsi.log  
setenv FOR103 $output/$name.nsi.cvg 
setenv FOR104 $output/$name.nsi.sol  
setenv FOR105  $input/$name.nsi.rst
setenv FOR106 $output/$name-element.nsi.set
setenv FOR107 $output/$name-boundary.nsi.set
setenv FOR108 $output/$name-node.nsi.set
setenv FOR109 $output/$name-ib.nsi.set
setenv FOR110 $output/$name.nsi.bcs
setenv FOR111 $output/$name-solver-mom.nsi.cvg 
setenv FOR112 $output/$name-sgs.nsi.log
setenv FOR113 $output/$name-sgs.nsi.cvg
setenv FOR114 $output/$name-matrix.nsi.ps
setenv FOR115 $output/$name-blogs.nsi.cvg
setenv FOR116 $output/$name-solver-con.nsi.cvg
setenv FOR117 $output/$name-reference.nsi.res
setenv FOR118 $output/$name-reference.nsi.cvg
setenv FOR119 $output/$name-continuity.nsi.sol  
setenv FOR120 $output/$name-low-mach.nsi.cvg  
setenv FOR122 $output/$name-dynamic-coupling.nsi.in
setenv FOR123 $output/$name-dynamic-coupling.nsi.out
setenv FOR124 $output/$name-dynamic-coupling.nsi.log
setenv FOR125 $output/$name-dynamic-coupling.nsi.res
setenv FOR126 $output/$name-boundary-conditions.nsi.sol 
#
# TEMPER: module 2
#
setenv FOR201  $input/$name.tem.dat         
setenv FOR202 $output/$name.tem.log
setenv FOR203 $output/$name.tem.cvg  
setenv FOR204 $output/$name.tem.sol  
setenv FOR205  $input/$name.tem.rst 
setenv FOR206 $output/$name-element.tem.set
setenv FOR207 $output/$name-boundary.tem.set
setenv FOR208 $output/$name-node.tem.set
setenv FOR210 $output/$name.tem.bcs
setenv FOR211 $output/$name.tem.cso
setenv FOR214 $output/$name-matrix.tem.ps 
setenv FOR215  $input/$name-conductivity.tem.fun
setenv FOR216  $input/$name-specificheat.tem.fun
setenv FOR220 $output/$name.tem.wit
setenv FOR221  $input/$name-bcinterpolation.tem.fix
setenv FOR222 $output/$name-dynamic-coupling.tem.in
setenv FOR223 $output/$name-dynamic-coupling.tem.out
setenv FOR224 $output/$name-dynamic-coupling.tem.log
setenv FOR225 $output/$name-dynamic-coupling.tem.res
setenv FOR232 $output/$name.tem.spl
setenv FOR233 $directory/$name-radiation.post.msh
setenv FOR234 $directory/$name-radiation.post.res 
setenv FOR235 $output/$name-low-mach.tem.cvg  
#
# CODIRE: module 3
#
setenv FOR301  $input/$name.cdr.dat
setenv FOR302 $output/$name.cdr.log
setenv FOR303 $output/$name.cdr.cvg
setenv FOR304 $output/$name.cdr.sol
setenv FOR305  $input/$name.cdr.rst
setenv FOR308 $output/$name.cdr.thp
setenv FOR309 $output/$name.cdr.lse
setenv FOR311 $output/$name.cdr.sex
setenv FOR312 $output/$name.cdr.sey
setenv FOR313 $output/$name.cdr.sez
setenv FOR321 $output/$name.cdr.tp1
setenv FOR322 $output/$name.cdr.tp2
setenv FOR323 $output/$name.cdr.tp3
setenv FOR324 $output/$name.cdr.tp4
setenv FOR325 $output/$name.cdr.tp5
setenv FOR330 $output/$name.cdr.flu
setenv FOR331 $output/$name.cdr.frc
setenv FOR332 $output/$name.cdr.spl
#
# TURBUL: module 4
#
setenv FOR401  $input/$name.tur.dat  
setenv FOR402 $output/$name.tur.log
setenv FOR403 $output/$name.tur.cvg 
setenv FOR404 $output/$name-variable1.tur.sol
setenv FOR405  $input/$name.tur.rst
setenv FOR406 $output/$name-variable2.tur.sol
setenv FOR407 $output/$name-variable3.tur.sol
setenv FOR408 $output/$name-variable4.tur.sol
setenv FOR409 $output/$name-wall-distance.tur.sol
setenv FOR410 $output/$name-variable1.tur.cso
setenv FOR411 $output/$name-variable2.tur.cso
setenv FOR412 $output/$name-variable3.tur.cso
setenv FOR413 $output/$name-variable4.tur.cso
setenv FOR414 $output/$name-wall-distance.tur.cso
#
# NASTAL: module 6
#
setenv FOR601  $input/$name.nsa.dat  
setenv FOR602 $output/$name.nsa.log
setenv FOR603 $output/$name.nsa.cvg 
setenv FOR604 $output/$name.nsa.sol
setenv FOR605 $output/$name.nsa.rst
setenv FOR606  $input/$name.nsa.frc
setenv FOR607  $input/$name.chk.in
setenv FOR608  $input/$name.chk.out
setenv FOR609  $input/$name.nsa.p2d
setenv FOR610 $output/$name.nsa.mxm
setenv FOR611 $output/$name-element.nsa.set
setenv FOR612 $output/$name-boundary.nsa.set
setenv FOR613 $output/$name-node.nsa.set
#
# ALEFOR: module 7
#
setenv FOR701  $input/$name.ale.dat  
setenv FOR702 $output/$name.ale.log
setenv FOR703 $output/$name.ale.cvg 
setenv FOR704 $output/$name.ale.sol
setenv FOR704 $input/$name.ale.rst
#
# SOLIDZ: module 10
#
setenv FOR1001  $input/$name.sld.dat         
setenv FOR1002 $output/$name.sld.log
setenv FOR1003 $output/$name.sld.cvg  
setenv FOR1004 $output/$name.sld.sol  
setenv FOR1005  $input/$name.sld.rst 
#
# GOTITA: module 11
#
setenv FOR1101  $input/$name.got.dat        
setenv FOR1102 $output/$name.got.log  
setenv FOR1103 $output/$name.got.cvg 
setenv FOR1104 $output/$name.got.sol  
setenv FOR1105  $input/$name.got.rst
setenv FOR1106 $output/$name-element.got.set
setenv FOR1107 $output/$name-boundary.got.set
setenv FOR1108 $output/$name-node.got.set
setenv FOR1110 $output/$name.got.bcs
setenv FOR1111 $output/$name.got.cso
setenv FOR1112 $output/$name.got.sgs
setenv FOR1113 $output/$name.got.csg
#
# WAVEQU: module 12
#
setenv FOR1201  $input/$name.wav.dat         
setenv FOR1202 $output/$name.wav.log
setenv FOR1203 $output/$name.wav.cvg  
setenv FOR1204 $output/$name.wav.sol  
setenv FOR1205  $input/$name.wav.rst 
setenv FOR1206 $output/$name-element.wav.set
setenv FOR1207 $output/$name-boundary.wav.set
setenv FOR1208 $output/$name-node.wav.set
#
# LEVELS: module 14
#
setenv FOR1401  $input/$name.lev.dat         
setenv FOR1402 $output/$name.lev.log
setenv FOR1403 $output/$name.lev.cvg  
setenv FOR1404 $output/$name.lev.sol  
setenv FOR1405  $input/$name.lev.rst 
setenv FOR1406 $output/$name-element.lev.set
setenv FOR1407 $output/$name-boundary.lev.set
setenv FOR1408 $output/$name-node.lev.set
setenv FOR1409 $output/$name-volume.lev.res
setenv FOR1410 $output/$name-gauge.lev.res
setenv FOR1411 $output/$name-velocity.lev.res
#
# QUANTY: module 15
#
setenv FOR1501  $input/$name.qua.dat         
setenv FOR1502 $output/$name.qua.log
setenv FOR1503 $output/$name.qua.cvg  
setenv FOR1504 $output/$name.qua.sol  
setenv FOR1505  $input/$name.qua.rst 
#
# PARTIS: module 17
#
setenv FOR1701  $input\$name.pts.dat         
setenv FOR1702 $output\$name.pts.log
setenv FOR1703 $output\$name.pts.cvg  
setenv FOR1704 $output\$name.pts.sol  
setenv FOR1705  $input\$name.pts.rst 
setenv FOR1706 $output\$name-element.pts.set
setenv FOR1707 $output\$name-boundary.pts.set
setenv FOR1708 $output\$name-node.pts.set
setenv FOR1709 $output\$name-ib.pts.set
setenv FOR1711 $output\$name.cso
setenv FOR1712  $input\$name-meteo.pts.cvg
setenv FOR1713  $input\$name-source.pts.cvg
setenv FOR1720 $output\$name.wit
#
# CHEMIC: module 19
#
setenv FOR1901  $input\$name.chm.dat         
setenv FOR1902 $output\$name.chm.log
setenv FOR1903 $output\$name.chm.cvg  
setenv FOR1904 $output\$name.chm.sol  
setenv FOR1905  $input\$name.chm.rst 
setenv FOR1906 $output\$name-element.chm.set
setenv FOR1907 $output\$name-boundary.chm.set
setenv FOR1908 $output\$name-node.chm.set
setenv FOR1909 $output\$name-ib.chm.set
setenv FOR1911 $output\$name.chm.cso
setenv FOR1920 $output\$name.chm.wit
#
# PARALL: service 55
#
setenv FOR5502 $output/$name.par.log
setenv FOR5503 $output/$name.par.post.msh
setenv FOR5505 $output/$name.par.mem
setenv FOR5506 $output/$name.par.cvg
setenv FOR5507 $output/$name-partition.par.post.msh
setenv FOR5508 $output/$name-partition.par.post.res
setenv FOR5515 $output/$name.par.post.res 
setenv FOR5516 $output/$name.par.rst
#
# DODEME: service 57
#
setenv FOR5701  $input/$name.dod.int
setenv FOR5702 $output/$name.dod.log

# -------------------------------------------------------------------
#   Execute version
# -------------------------------------------------------------------

if ( $version == "debug" )   goto rundebug
if ( $version == "release" ) goto runrelease
goto help

rundebug:
echo 'ALYA * * * DEBUG VERSION   * * *' $name
db -g %ALYA_DIR%/Executables/dos/Alya.g
goto end

runrelease:
set rexe=$ALYA_DIR/Executables/unix/Alya.x
if ( ! -e $rexe ) goto stoprexe
echo Running $rexe
$rexe
goto end

help:
echo '--|'
echo '--|   Usage:'
echo '--| '
echo '--|   Alya.bat [-g,-x,-help] problem [option1] [option2]'
echo '--| '
echo '--|   -x ........ Release version (default)'
echo '--|   -g ........ Debug version '
echo '--|   -help ..... Displays this help '
echo '--| '
echo '--|   option1 ... If it is present, the script understands that it was called by GiD.'
echo '--|               It is the problem directory.'
echo '--|   option2 ... If option1 is present, it is the GiD problem type.'
echo '--|  '
echo '--|   Runs Alya. '
echo '--|   The problem data must be located in directory: ./data/'
echo '--|   The problem results are dumped in directory:   ./results/'
echo '--|  '
stop

stoprexe:
echo 'Could not find executable' $rexe

end:
