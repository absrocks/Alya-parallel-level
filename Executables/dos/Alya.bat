@echo off
rem -------------------------------------------------------------------
rem
rem                          A  L  Y  A
rem
rem -------------------------------------------------------------------

rem -------------------------------------------------------------------
rem   Identify arguments
rem -------------------------------------------------------------------

set version=release
if "%1%" == "-d" goto debug
if "%1%" == "-r" goto release
if "%1%" == "-g" goto debug
if "%1%" == "-x" goto release
goto readarguments

:debug
set version=debug
shift
goto readarguments

:release
set version=release
shift
goto readarguments

:readarguments
set name=%1
set directory=%2
set ProblemType=%3

if "%1%" == ""       goto help
if "%1%" == "--help" goto help
if "%1%" == "-help"  goto help
if "%1%" == "-h"     goto help
if "%1%" == "/?"     goto help
if "%2%" == ""       goto emptydir

rem -------------------------------------------------------------------
rem   Script was called from GiD: creates the directory 
rem -------------------------------------------------------------------

set gid=1
set input=%directory%\data
set output=%directory%\results
rem OutputFile: %output%\%name%.liv 
if not exist %input%  mkdir %input%
if not exist %output% mkdir %output%
cd %input%

rem -------------------------------------------------------------------
rem   Remove existing files
rem -------------------------------------------------------------------

if   exist %input%\%name%.dat       del /Q /S /F %input%\%name%.dat    
if   exist %input%\%name%.dom.dat   del /Q /S /F %input%\%name%.dom.dat
if   exist %input%\%name%.set.dat   del /Q /S /F %input%\%name%.set.dat
if   exist %input%\%name%.geo.dat   del /Q /S /F %input%\%name%.geo.dat
if   exist %input%\%name%.nsi.dat   del /Q /S /F %input%\%name%.nsi.dat
if   exist %input%\%name%.nsi.fix   del /Q /S /F %input%\%name%.nsi.fix
if   exist %input%\%name%.tem.dat   del /Q /S /F %input%\%name%.tem.dat
if   exist %input%\%name%.tem.fix   del /Q /S /F %input%\%name%.tem.fix
if   exist %input%\%name%.cdr.dat   del /Q /S /F %input%\%name%.cdr.dat
if   exist %input%\%name%.cdr.fix   del /Q /S /F %input%\%name%.cdr.fix
if   exist %input%\%name%.tur.dat   del /Q /S /F %input%\%name%.tur.dat
if   exist %input%\%name%.tur.fix   del /Q /S /F %input%\%name%.tur.fix
if   exist %input%\%name%.exm.dat   del /Q /S /F %input%\%name%.exm.dat
if   exist %input%\%name%.exm.fix   del /Q /S /F %input%\%name%.exm.fix
if   exist %input%\%name%.nsa.dat   del /Q /S /F %input%\%name%.nsa.dat
if   exist %input%\%name%.nsa.fix   del /Q /S /F %input%\%name%.nsa.fix
if   exist %input%\%name%.ale.dat   del /Q /S /F %input%\%name%.ale.dat
if   exist %input%\%name%.ale.fix   del /Q /S /F %input%\%name%.ale.fix
if   exist %input%\%name%.sld.dat   del /Q /S /F %input%\%name%.sld.dat
if   exist %input%\%name%.sld.fix   del /Q /S /F %input%\%name%.sld.fix
if   exist %input%\%name%.got.dat   del /Q /S /F %input%\%name%.got.dat
if   exist %input%\%name%.got.fix   del /Q /S /F %input%\%name%.got.fix
if   exist %input%\%name%.dod.int   del /Q /S /F %input%\%name%.dod.int
if   exist %input%\%name%.dod.fix   del /Q /S /F %input%\%name%.dod.fix
if   exist %input%\%name%.wav.dat   del /Q /S /F %input%\%name%.wav.dat
if   exist %input%\%name%.wav.fix   del /Q /S /F %input%\%name%.wav.fix
if   exist %input%\%name%.lev.dat   del /Q /S /F %input%\%name%.lev.dat
if   exist %input%\%name%.lev.fix   del /Q /S /F %input%\%name%.lev.fix

rem -------------------------------------------------------------------
rem   Rename the files of the Master
rem -------------------------------------------------------------------

move %directory%\%name%.dat     %input%\%name%.dat
move %directory%\%name%-1.dat   %input%\%name%.dom.dat
move %directory%\%name%-2.dat   %input%\%name%.geo.dat
move %directory%\%name%-3.dat   %input%\%name%.set.dat

rem -------------------------------------------------------------------
rem   Rename the files of the modules
rem -------------------------------------------------------------------

move %directory%\%name%-4.dat   %input%\%name%.nsi.dat
move %directory%\%name%-5.dat   %input%\%name%.nsi.fix
move %directory%\%name%-6.dat   %input%\%name%.tem.dat
move %directory%\%name%-7.dat   %input%\%name%.tem.fix
move %directory%\%name%-8.dat   %input%\%name%.cdr.dat
move %directory%\%name%-9.dat   %input%\%name%.cdr.fix
move %directory%\%name%-10.dat  %input%\%name%.tur.dat
move %directory%\%name%-11.dat  %input%\%name%.tur.fix
move %directory%\%name%-12.dat  %input%\%name%.exm.dat
move %directory%\%name%-13.dat  %input%\%name%.exm.fix
move %directory%\%name%-14.dat  %input%\%name%.nsa.dat
move %directory%\%name%-15.dat  %input%\%name%.nsa.fix
move %directory%\%name%-16.dat  %input%\%name%.ale.dat
move %directory%\%name%-17.dat  %input%\%name%.ale.fix
move %directory%\%name%-18.dat  %input%\%name%.sld.dat
move %directory%\%name%-19.dat  %input%\%name%.sld.fix
move %directory%\%name%-20.dat  %input%\%name%.got.dat
move %directory%\%name%-21.dat  %input%\%name%.got.fix
move %directory%\%name%-22.dat  %input%\%name%.dod.int
move %directory%\%name%-23.dat  %input%\%name%.dod.fix
move %directory%\%name%-24.dat  %input%\%name%.wav.dat
move %directory%\%name%-25.dat  %input%\%name%.wav.fix
move %directory%\%name%-26.dat  %input%\%name%.lev.dat
move %directory%\%name%-27.dat  %input%\%name%.lev.fix

goto continue

rem -------------------------------------------------------------------
rem   Script was called from DOS: create result directory
rem -------------------------------------------------------------------
:emptydir
set directory=%cd%
set gid=0
set input=%directory%\data
set output=%directory%\results
if not exist %output% mkdir %output%

:continue
rem -------------------------------------------------------------------
rem   Define environment variables
rem -------------------------------------------------------------------

rem
rem Master
rem
set ALYA_NAME=%name%
set FOR011= %input%\%name%.dat           
set FOR012=%output%\%name%.log       
set FOR013=%output%\%name%.mem
set FOR014=%output%\%name%.cvg
set FOR015=%output%\%name%
set FOR016=%output%\%name%.liv
set FOR017= %input%\%name%.rst    
set FOR018=%output%\%name%-latex.tex    
set FOR019=%output%\%name%-latex.plt 
set FOR021= %input%\%name%.dom.dat 
set FOR023=%output%\%name%-DD.post.msh
set FOR024=%output%\%name%.dom.bin
set FOR025=%output%\%name%-elsest.log 
set FOR026=%output%\%name%-elsest.post.msh
set FOR027=%output%\%name%-elsest.post.res
set FOR028=%output%\%name%-system.log
set FOR029=%output%\%name%-IB.res
set FOR030=%output%\%name%.ensi.case
set FOR034=%output%\%name%
set FOR039=%output%\%name%
set FOR045= %input%\%name%-IB.rst    
set FOR046=%output%\%name%
set FOR047=%output%\%name%
set FOR048=%output%\%name%
set FOR049=%output%\%name%
if "%gid%" == "1" set FOR015=%directory%\%name%
if "%gid%" == "1" set FOR034=%directory%\%name%
if "%gid%" == "1" set FOR039=%directory%\%name%
rem
rem NASTIN: module 1
rem
set FOR101= %input%\%name%.nsi.dat        
set FOR102=%output%\%name%.nsi.log  
set FOR103=%output%\%name%.nsi.cvg 
set FOR104=%output%\%name%.nsi.sol  
set FOR105= %input%\%name%.nsi.rst
set FOR106=%output%\%name%-element.nsi.set
set FOR107=%output%\%name%-boundary.nsi.set
set FOR108=%output%\%name%-node.nsi.set
set FOR109=%output%\%name%-ib.nsi.set
set FOR110=%output%\%name%.nsi.bcs
set FOR111=%output%\%name%-solver-mom.nsi.cvg
set FOR112=%output%\%name%-sgs.nsi.log
set FOR113=%output%\%name%-sgs.nsi.cvg
set FOR114=%output%\%name%-matrix.nsi.ps
set FOR115=%output%\%name%-blogs.nsi.cvg
set FOR116=%output%\%name%-solver-con.nsi.cvg
set FOR117=%output%\%name%-reference.nsi.res
set FOR118=%output%\%name%-reference.nsi.cvg
set FOR119=%output%\%name%-continuity.nsi.sol
set FOR120=%output%\%name%-low-mach.nsi.cvg
set FOR122=%output%\%name%-dynamic-coupling.nsi.in
set FOR123=%output%\%name%-dynamic-coupling.nsi.out
set FOR124=%output%\%name%-dynamic-coupling.nsi.log
set FOR125=%output%\%name%-dynamic-coupling.nsi.res
set FOR126=%output%\%name%-boundary-conditions.nsi.sol
rem
rem TEMPER: module 2
rem
set FOR201= %input%\%name%.tem.dat         
set FOR202=%output%\%name%.tem.log
set FOR203=%output%\%name%.tem.cvg  
set FOR204=%output%\%name%.tem.sol  
set FOR205= %input%\%name%.tem.rst 
set FOR206=%output%\%name%-element.tem.set
set FOR207=%output%\%name%-boundary.tem.set
set FOR208=%output%\%name%-node.tem.set
set FOR210=%output%\%name%.tem.bcs
set FOR211=%output%\%name%.tem.cso
set FOR214=%output%\%name%-matrix.tem.ps 
set FOR215= %input%\%name%-conductivity.tem.fun
set FOR216= %input%\%name%-specificheat.tem.fun
set FOR220=%output%\%name%.tem.wit
set FOR221= %input%\%name%-bcinterpolation.tem.fix
set FOR222=%output%\%name%-dynamic-coupling.tem.in
set FOR223=%output%\%name%-dynamic-coupling.tem.out
set FOR224=%output%\%name%-dynamic-coupling.tem.log
set FOR225=%output%\%name%-dynamic-coupling.tem.res
set FOR232=%output%\%name%.tem.spl
set FOR233=%directory%\%name%-radiation.post.msh
set FOR234=%directory%\%name%-radiation.post.res 
rem
rem CODIRE: module 3
rem
set FOR301= %input%\%name%.cdr.dat
set FOR302=%output%\%name%.cdr.log
set FOR303=%output%\%name%.cdr.cvg
set FOR304=%output%\%name%.cdr.sol
set FOR305= %input%\%name%.cdr.rst
set FOR308=%output%\%name%.cdr.thp
set FOR309=%output%\%name%.cdr.lse
set FOR311=%output%\%name%.cdr.sex
set FOR312=%output%\%name%.cdr.sey
set FOR313=%output%\%name%.cdr.sez
set FOR321=%output%\%name%.cdr.tp1
set FOR322=%output%\%name%.cdr.tp2
set FOR323=%output%\%name%.cdr.tp3
set FOR324=%output%\%name%.cdr.tp4
set FOR325=%output%\%name%.cdr.tp5
set FOR330=%output%\%name%.cdr.flu
set FOR331=%output%\%name%.cdr.frc
set FOR332=%output%\%name%.cdr.spl
rem
rem TURBUL: module 4
rem
set FOR401= %input%\%name%.tur.dat  
set FOR402=%output%\%name%.tur.log
set FOR403=%output%\%name%.tur.cvg 
set FOR404=%output%\%name%-variable1.tur.sol
set FOR405= %input%\%name%.tur.rst
set FOR406=%output%\%name%-variable2.tur.sol
set FOR407=%output%\%name%-variable3.tur.sol
set FOR408=%output%\%name%-variable4.tur.sol
set FOR409=%output%\%name%-wall-distance.tur.sol
set FOR410=%output%\%name%-variable1.tur.cso
set FOR411=%output%\%name%-variable2.tur.cso
set FOR412=%output%\%name%-variable3.tur.cso
set FOR413=%output%\%name%-variable4.tur.cso
set FOR414=%output%\%name%-wall-distance.tur.cso
rem
rem NASTAL: module 6
rem
set FOR601= %input%\%name%.nsa.dat  
set FOR602=%output%\%name%.nsa.log
set FOR603=%output%\%name%.nsa.cvg 
set FOR604=%output%\%name%.nsa.sol
set FOR605=%output%\%name%.nsa.rst
set FOR606= %input%\%name%.nsa.frc
set FOR607= %input%\%name%.chk.in
set FOR608= %input%\%name%.chk.out
set FOR609= %input%\%name%.nsa.p2d
set FOR610=%output%\%name%.nsa.mxm
set FOR611=%output%\%name%-element.nsa.set
set FOR612=%output%\%name%-boundary.nsa.set
set FOR613=%output%\%name%-node.nsa.set
rem
rem ALEFOR: module 7
rem
set FOR701= %input%\%name%.ale.dat  
set FOR702=%output%\%name%.ale.log
set FOR703=%output%\%name%.ale.cvg 
set FOR704=%output%\%name%.ale.sol
set FOR704= %input%\%name%.ale.rst
rem
rem SOLIDZ: module 10
rem
set FOR1001= %input%\%name%.sld.dat         
set FOR1002=%output%\%name%.sld.log
set FOR1003=%output%\%name%.sld.cvg  
set FOR1004=%output%\%name%.sld.sol  
set FOR1005= %input%\%name%.sld.rst 
rem
rem GOTITA: module 11
rem
set FOR1101= %input%\%name%.got.dat        
set FOR1102=%output%\%name%.got.log  
set FOR1103=%output%\%name%.got.cvg 
set FOR1104=%output%\%name%.got.sol  
set FOR1105= %input%\%name%.got.rst
set FOR1106=%output%\%name%-element.got.set
set FOR1107=%output%\%name%-boundary.got.set
set FOR1108=%output%\%name%-node.got.set
set FOR1110=%output%\%name%.got.bcs
set FOR1111=%output%\%name%.got.cso
set FOR1112=%output%\%name%.got.sgs
set FOR1113=%output%\%name%.got.csg
rem
rem WAVEQU: module 12
rem
set FOR1201= %input%\%name%.wav.dat         
set FOR1202=%output%\%name%.wav.log
set FOR1203=%output%\%name%.wav.cvg  
set FOR1204=%output%\%name%.wav.sol  
set FOR1205= %input%\%name%.wav.rst 
set FOR1206=%output%\%name%-element.wav.set
set FOR1207=%output%\%name%-boundary.wav.set
set FOR1208=%output%\%name%-node.wav.set
rem
rem LEVELS: module 14
rem
set FOR1401= %input%\%name%.lev.dat         
set FOR1402=%output%\%name%.lev.log
set FOR1403=%output%\%name%.lev.cvg  
set FOR1404=%output%\%name%.lev.sol  
set FOR1405= %input%\%name%.lev.rst 
set FOR1406=%output%\%name%-element.lev.set
set FOR1407=%output%\%name%-boundary.lev.set
set FOR1408=%output%\%name%-node.lev.set
set FOR1409=%output%\%name%-volume.lev.res
set FOR1410=%output%\%name%-gauge.lev.res
set FOR1411=%output%\%name%-velocity.lev.res
rem
rem QUANTY: module 15
rem
set FOR1501= %input%\%name%.qua.dat         
set FOR1502=%output%\%name%.qua.log
set FOR1503=%output%\%name%.qua.cvg  
set FOR1504=%output%\%name%.qua.sol  
set FOR1505= %input%\%name%.qua.rst 
rem
rem PARTIS: module 17
rem
set FOR1701= %input%\%name%.pts.dat         
set FOR1702=%output%\%name%.pts.log
set FOR1703=%output%\%name%.pts.cvg  
set FOR1704=%output%\%name%.pts.sol  
set FOR1705= %input%\%name%.pts.rst 
set FOR1706=%output%\%name%-element.pts.set
set FOR1707=%output%\%name%-boundary.pts.set
set FOR1708=%output%\%name%-node.pts.set
set FOR1709=%output%\%name%-ib.pts.set
set FOR1711=%output%\%name%.pts.cso
set FOR1712= %input%\%name%-meteo.pts.cvg
set FOR1713= %input%\%name%-source.pts.cvg
set FOR1720=%output%\%name%.pts.wit
rem
rem CHEMIC: module 19
rem
set FOR1901= %input%\%name%.chm.dat         
set FOR1902=%output%\%name%.chm.log
set FOR1903=%output%\%name%.chm.cvg  
set FOR1904=%output%\%name%.chm.sol  
set FOR1905= %input%\%name%.chm.rst 
set FOR1906=%output%\%name%-element.chm.set
set FOR1907=%output%\%name%-boundary.chm.set
set FOR1908=%output%\%name%-node.chm.set
set FOR1909=%output%\%name%-ib.chm.set
set FOR1911=%output%\%name%.cso
set FOR1920=%output%\%name%.wit
rem
rem PARALL: service 55
rem
set FOR5502=%output%\%name%.par.log
set FOR5503=%output%\%name%.par.post.msh
set FOR5505=%output%\%name%.par.mem
set FOR5506=%output%\%name%.par.cvg
set FOR5507=%output%\%name%-partition.par.post.msh
set FOR5508=%output%\%name%-partition.par.post.res
set FOR5515=%output%\%name%.par.post.res 
set FOR5516=%output%\%name%.par.rst
rem
rem DODEME: service 57
rem
set FOR5701= %input%\%name%.dod.int
set FOR5702=%output%\%name%.dod.log

rem -------------------------------------------------------------------
rem   Execute version
rem -------------------------------------------------------------------

if "%version%" == "debug"   goto rundebug
if "%version%" == "release" goto runrelease
goto help

:rundebug
echo ALYA * * * DEBUG VERSION   * * * %name%
dfdev.exe %ALYA_DIR%\Executables\dos\Alya.dsw
goto end

:runrelease
echo ALYA * * * RELEASE VERSION * * * %name%
%ALYA_DIR%\Executables\dos\release\Alya.exe
goto end

:help
echo $---------------------------------------------------------------------------------
echo $ 
echo $  Usage:
echo $
echo $  Alya.bat [-g,-x,-help] problem [option1] [option2]
echo $
echo $  -x ........ Release version (default)
echo $  -g ........ Debug version 
echo $  -help ..... Displays this help 
echo $
echo $  option1 ... If it is present, the script understands that it was called by GiD.
echo $              It is the problem directory.
echo $  option2 ... If option1 is present, it is the GiD problem type.
echo $ 
echo $  Runs Alya. 
echo $  The problem data must be located in directory: .\data\
echo $  The problem results are dumped in directory:   .\results\
echo $ 
echo $---------------------------------------------------------------------------------

:end
