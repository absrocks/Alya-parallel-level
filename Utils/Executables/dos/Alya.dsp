# Microsoft Developer Studio Project File - Name="Alya" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=Alya - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "Alya.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "Alya.mak" CFG="Alya - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "Alya - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "Alya - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 1
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "Alya - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE F90 /compile_only /nologo /warn:nofileopt
# ADD F90 /alignment:commons /assume:noaccuracy_sensitive /browser /compile_only /define:MPI_OFF=1 /define:DOS=1 /fpp /include:"..\..\..\Mumps\libseq" /include:"..\..\..\Mumps\include" /nologo /real_size:64 /transform_loops /warn:nofileopt
# SUBTRACT F90 /assume:buffered_io /check:arg_temp_created /check:bounds /check:format /check:output_conversion /check:overflow /check:underflow /traceback /fast
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /I "services\solmum" /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /FR /YX /FD /c
# ADD BASE RSC /l 0xc0a /d "NDEBUG"
# ADD RSC /l 0xc0a /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib libgoto_p4_256-r0.94.lib gidpost.lib libmumps.lib libmetis.a /nologo /stack:0x3b9aca00 /subsystem:console /profile /map /machine:I386 /libpath:".\Libraries\\" /libpath:"..\..\..\Mumps\Release" /libpath:"..\..\..\Blas\Goto" /libpath:"..\..\..\Metis"
# SUBTRACT LINK32 /verbose

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE F90 /check:bounds /compile_only /dbglibs /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD F90 /browser /check:arg_temp_created /check:bounds /check:format /check:power /check:output_conversion /check:overflow /check:underflow /compile_only /dbglibs /debug:full /define:MPI_OFF=1 /define:DOS=1 /fpe:0 /fpp /include:".\services\solmum" /include:"..\..\..\Mumps\libseq" /include:"..\..\..\Mumps\include" /include:"services\solmum" /nologo /real_size:64 /stand:f90 /traceback /warn:argument_checking /warn:declarations /warn:nofileopt /warn:unused
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /FR /YX /FD /GZ /c
# ADD BASE RSC /l 0xc0a /d "_DEBUG"
# ADD RSC /l 0xc0a /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib gidpost.lib libgoto_p4_256-r0.94.lib libmumps.lib libmetis.a PLS.lib /nologo /stack:0x3b9aca00 /subsystem:console /incremental:no /debug /machine:I386 /pdbtype:sept /libpath:".\Libraries\\" /libpath:"..\..\..\Mumps\Debug" /libpath:"..\..\..\Blas\Goto" /libpath:"..\..\..\Metis" /libpath:".\..\..\..\..\MsDev\PLS\Debug"
# SUBTRACT LINK32 /profile

!ENDIF 

# Begin Target

# Name "Alya - Win32 Release"
# Name "Alya - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat;f90;for;f;fpp"
# Begin Group "kernel"

# PROP Default_Filter "f90"
# Begin Group "wtools"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\Sources\kernel\wtools\codfix.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CODFI=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CODFI=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\wtools\cputim.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CPUTI=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CPUTI=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\wtools\ecoute.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ECOUT=\
	".\Release\def_inpout.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ECOUT=\
	".\Debug\def_inpout.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\wtools\iexcha.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_IEXCH=\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_IEXCH=\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\wtools\matrea.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_MATRE=\
	".\Release\def_inpout.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_MATRE=\
	".\Debug\def_inpout.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\wtools\mod_htable.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_MOD_H=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_MOD_H=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\wtools\rexcha.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_REXCH=\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_REXCH=\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\wtools\vecrea.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_VECRE=\
	".\Release\def_inpout.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_VECRE=\
	".\Debug\def_inpout.mod"\
	

!ENDIF 

# End Source File
# End Group
# Begin Group "soldir"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\Sources\kernel\soldir\fvecdo.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_FVECD=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_FVECD=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\soldir\rencon.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_RENCO=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_RENCO=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\soldir\renum0.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_RENUM=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_RENUM=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\soldir\renum1.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_RENUM1=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_RENUM1=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\soldir\renum2.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_RENUM2=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_RENUM2=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\soldir\renum3.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_RENUM3=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_RENUM3=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\soldir\renum4.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_RENUM4=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_RENUM4=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\soldir\renumb.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_RENUMB=\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_RENUMB=\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\soldir\renumn.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_RENUMN=\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_RENUMN=\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\soldir\skyase.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SKYAS=\
	".\Release\def_kintyp.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SKYAS=\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\soldir\skybak.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SKYBA=\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SKYBA=\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\soldir\skycek.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SKYCE=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SKYCE=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\soldir\skydia.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SKYDI=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SKYDI=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\soldir\skyini.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SKYIN=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SKYIN=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\soldir\skylin.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SKYLI=\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SKYLI=\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\soldir\skylpo.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SKYLP=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SKYLP=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\soldir\skypak.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SKYPA=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SKYPA=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\soldir\skyplu.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SKYPL=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SKYPL=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\soldir\skyren.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SKYRE=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SKYRE=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\soldir\skytri.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SKYTR=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SKYTR=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\soldir\soldir.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SOLDI=\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SOLDI=\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# End Group
# Begin Group "outrut"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\Sources\kernel\outrut\chanum.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CHANU=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CHANU=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\outrut\engold_initia.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ENGOL=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_postpr.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ENGOL=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_postpr.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\outrut\engold_openfi.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ENGOLD=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_postpr.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ENGOLD=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_postpr.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\outrut\engold_wrcase.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ENGOLD_=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_postpr.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ENGOLD_=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_postpr.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\outrut\geobin.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GEOBI=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_postpr.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GEOBI=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_postpr.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\outrut\geoens.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GEOEN=\
	".\Release\def_domain.mod"\
	".\Release\def_elmtyp.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_postpr.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GEOEN=\
	".\Debug\def_domain.mod"\
	".\Debug\def_elmtyp.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_postpr.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\outrut\geofem.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GEOFE=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GEOFE=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\outrut\geogid.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GEOGI=\
	".\Release\def_domain.mod"\
	".\Release\def_elmtyp.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GEOGI=\
	".\Debug\def_domain.mod"\
	".\Debug\def_elmtyp.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\outrut\geomvu.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GEOMV=\
	".\Release\def_domain.mod"\
	".\Release\def_elmtyp.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_postpr.mod"\
	".\Release\Mod_iofile.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GEOMV=\
	".\Debug\def_domain.mod"\
	".\Debug\def_elmtyp.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_postpr.mod"\
	".\Debug\Mod_iofile.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\outrut\livinf.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_LIVIN=\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_LIVIN=\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\outrut\mod_output.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_MOD_O=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_MOD_O=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\outrut\mod_postpr.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_MOD_P=\
	".\Release\def_domain.mod"\
	".\Release\def_elmtyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_postpr.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_MOD_P=\
	".\Debug\def_domain.mod"\
	".\Debug\def_elmtyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_postpr.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\outrut\outcpu.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_OUTCP=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_OUTCP=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\outrut\outerr.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_OUTER=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_OUTER=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\outrut\outfor.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_OUTFO=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_OUTFO=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\outrut\outhel.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_OUTHE=\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_OUTHE=\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\outrut\outinf.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_OUTIN=\
	".\Release\def_domain.mod"\
	".\Release\def_elmtyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_OUTIN=\
	".\Debug\def_domain.mod"\
	".\Debug\def_elmtyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\outrut\outlat.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_OUTLA=\
	".\Release\def_domain.mod"\
	".\Release\def_elmtyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_OUTLA=\
	".\Debug\def_domain.mod"\
	".\Debug\def_elmtyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\outrut\outmem.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_OUTME=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_meshin.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_OUTME=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_meshin.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\outrut\outres.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_OUTRE=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_postpr.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_OUTRE=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_postpr.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\outrut\outset.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_OUTSE=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_OUTSE=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\outrut\outvar.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_OUTVA=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_postpr.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_OUTVA=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_postpr.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\outrut\posr3p.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_POSR3=\
	".\Release\def_domain.mod"\
	".\Release\def_elmtyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_POSR3=\
	".\Debug\def_domain.mod"\
	".\Debug\def_elmtyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\outrut\pspltm.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PSPLT=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PSPLT=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\outrut\realen.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_REALE=\
	".\Release\def_inpout.mod"\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_REALE=\
	".\Debug\def_inpout.mod"\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\outrut\reasol.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_REASO=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_REASO=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\outrut\reatau.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_REATA=\
	".\Release\def_inpout.mod"\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_REATA=\
	".\Debug\def_inpout.mod"\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\outrut\strcub.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_STRCU=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_STRCU=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\outrut\strfun.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_STRFU=\
	".\Release\def_domain.mod"\
	".\Release\def_elmtyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_STRFU=\
	".\Debug\def_domain.mod"\
	".\Debug\def_elmtyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\outrut\strloc.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_STRLO=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_STRLO=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\outrut\suplot.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SUPLO=\
	".\Release\def_domain.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SUPLO=\
	".\Debug\def_domain.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# End Group
# Begin Group "memory"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\Sources\kernel\memory\memctr.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_MEMCT=\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_MEMCT=\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\memory\memerr.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_MEMER=\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_MEMER=\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\memory\memgen.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_MEMGE=\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_MEMGE=\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\memory\mod_memchk.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_MOD_M=\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_MOD_M=\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# End Group
# Begin Group "mathru"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\Sources\kernel\mathru\adrmat.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ADRMA=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ADRMA=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\mathru\adrres.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ADRRE=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ADRRE=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\mathru\adrtes.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ADRTE=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ADRTE=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\mathru\aitken.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_AITKE=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_AITKE=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\mathru\assmat.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ASSMA=\
	".\Release\def_kintyp.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ASSMA=\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\mathru\assrhs.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ASSRH=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ASSRH=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\mathru\bcsrsc.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_BCSRS=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_BCSRS=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\mathru\bcsrsm.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_BCSRSM=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_BCSRSM=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\mathru\btdbma.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_BTDBM=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_BTDBM=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\mathru\elmadr.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ELMAD=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ELMAD=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\mathru\elmce2.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ELMCE=\
	".\Release\def_domain.mod"\
	".\Release\def_elmtyp.mod"\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ELMCE=\
	".\Debug\def_domain.mod"\
	".\Debug\def_elmtyp.mod"\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\mathru\elmcel.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ELMCEL=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ELMCEL=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\mathru\elmdir.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ELMDI=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ELMDI=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\mathru\elmshc.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ELMSH=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ELMSH=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\mathru\elmtss.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ELMTS=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ELMTS=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\mathru\frivel.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_FRIVE=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_FRIVE=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\mathru\funbcs.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_FUNBC=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_FUNBC=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\mathru\funcrd.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_FUNCR=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_FUNCR=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\mathru\funcre.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_FUNCRE=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_FUNCRE=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\mathru\gather.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GATHE=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GATHE=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\mathru\heapsorti1.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_HEAPS=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_HEAPS=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\mathru\invert.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_INVER=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_INVER=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\mathru\invmtx.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_INVMT=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_INVMT=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\mathru\linint.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_LININ=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_LININ=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\mathru\matove.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_MATOV=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_MATOV=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\mathru\mbmab0.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_MBMAB=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_MBMAB=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\mathru\mbmabt.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_MBMABT=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_MBMABT=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\mathru\mbmatb.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_MBMAT=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_MBMAT=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\mathru\mbvab0.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_MBVAB=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_MBVAB=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\mathru\mbvab1.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_MBVAB1=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_MBVAB1=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\mathru\mbvab2.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_MBVAB2=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_MBVAB2=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\mathru\mbvatb.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_MBVAT=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_MBVAT=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\mathru\minmax.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_MINMA=\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_MINMA=\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\mathru\mixing.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_MIXIN=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_MIXIN=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\mathru\mixmax.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_MIXMA=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_MIXMA=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\mathru\mod_gradie.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_MOD_G=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_MOD_G=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\mathru\ordena.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ORDEN=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ORDEN=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\mathru\parbdf.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PARBD=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PARBD=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\mathru\residu.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_RESID=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_RESID=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\mathru\scapro.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SCAPR=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SCAPR=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\mathru\sortin.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SORTI=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SORTI=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\mathru\tauadr.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TAUAD=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TAUAD=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\mathru\typarr.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TYPAR=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TYPAR=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\mathru\vecasi.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_VECAS=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_VECAS=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\mathru\vecnor.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_VECNO=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_VECNO=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\mathru\vecpro.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_VECPR=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_VECPR=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\mathru\vecres.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_VECRES=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_VECRES=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\mathru\vecuni.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_VECUN=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_VECUN=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\mathru\vetoma.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_VETOM=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_VETOM=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\mathru\vortic.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_VORTI=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_VORTI=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# End Group
# Begin Group "solite"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\Sources\kernel\solite\bcgpls.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_BCGPL=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_solver.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_BCGPL=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_solver.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\solite\bcsrax.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_BCSRA=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_BCSRA=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\solite\bsymax.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_BSYMA=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_BSYMA=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\solite\cgrpls.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CGRPL=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_solver.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CGRPL=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_solver.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\solite\csrase.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CSRAS=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CSRAS=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\solite\csrper.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CSRPE=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CSRPE=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\solite\deflcg.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_DEFLC=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_solver.mod"\
	".\Release\Mod_memchk.mod"\
	".\Release\Mod_postpr.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_DEFLC=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_solver.mod"\
	".\Debug\Mod_memchk.mod"\
	".\Debug\Mod_postpr.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\soldir\dpstrf.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_DPSTR=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_DPSTR=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\soldir\dpstrs.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_DPSTRS=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_DPSTRS=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\solite\gmres.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GMRES=\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GMRES=\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\solite\gmrort.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GMROR=\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GMROR=\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\solite\gmrpls.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GMRPL=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_solver.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GMRPL=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_solver.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\solite\linlet.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_LINLE=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	".\Release\Mod_memchk.mod"\
	".\Release\Mod_postpr.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_LINLE=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	".\Debug\Mod_memchk.mod"\
	".\Debug\Mod_postpr.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\solite\nor22x.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NOR22=\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NOR22=\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\solite\norm1x.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NORM1=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NORM1=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\solite\norm2x.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NORM2=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NORM2=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\solite\prodts.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PRODT=\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PRODT=\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\solite\prodxy.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PRODX=\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PRODX=\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\solite\resvec.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_RESVE=\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_RESVE=\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\solite\rhsper.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_RHSPE=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_RHSPE=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\solite\soldia.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SOLDIA=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_solver.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SOLDIA=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_solver.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\solite\solite.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SOLIT=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SOLIT=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\solite\solric.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SOLRI=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_solver.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SOLRI=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_solver.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# End Group
# Begin Group "domain"

# PROP Default_Filter "f90"
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\arrinb.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ARRIN=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ARRIN=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\arrind.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ARRIND=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ARRIND=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\boucar.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_BOUCA=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_BOUCA=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\bouder.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_BOUDE=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_BOUDE=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\bouele.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_BOUEL=\
	".\Release\def_elmtyp.mod"\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_BOUEL=\
	".\Debug\def_elmtyp.mod"\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\bounor.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_BOUNO=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_BOUNO=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\cartbo.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CARTB=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CARTB=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\cderda.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CDERD=\
	".\Release\def_domain.mod"\
	".\Release\def_elmtyp.mod"\
	".\Release\def_master.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CDERD=\
	".\Debug\def_domain.mod"\
	".\Debug\def_elmtyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\chaord.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CHAOR=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CHAOR=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\chenor.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CHENO=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CHENO=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\conmas.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CONMA=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CONMA=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\connbo.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CONNB=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CONNB=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\connpo.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CONNP=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CONNP=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\corner.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CORNE=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CORNE=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\cregro.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CREGR=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	".\Release\Mod_memchk.mod"\
	".\Release\Mod_postpr.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CREGR=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	".\Debug\Mod_memchk.mod"\
	".\Debug\Mod_postpr.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\cresla.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CRESL=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CRESL=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\cshder.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CSHDE=\
	".\Release\def_domain.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CSHDE=\
	".\Debug\def_domain.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\d2sdx2.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_D2SDX=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_D2SDX=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\defmsh.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_DEFMS=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_DEFMS=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\domain.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_DOMAI=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_DOMAI=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\dombou.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_DOMBO=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_DOMBO=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\domfa2.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_DOMFA=\
	".\Release\def_elmtyp.mod"\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_DOMFA=\
	".\Debug\def_elmtyp.mod"\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\domfac.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_DOMFAC=\
	".\Release\def_domain.mod"\
	".\Release\def_elmtyp.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_DOMFAC=\
	".\Debug\def_domain.mod"\
	".\Debug\def_elmtyp.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\domgra.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_DOMGR=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_DOMGR=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\dompla.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_DOMPL=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	".\Release\Mod_postpr.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_DOMPL=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	".\Debug\Mod_postpr.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\domset.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_DOMSE=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	".\Release\Mod_postpr.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_DOMSE=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	".\Debug\Mod_postpr.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\domvar.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_DOMVA=\
	".\Release\def_domain.mod"\
	".\Release\def_elmtyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_DOMVA=\
	".\Debug\def_domain.mod"\
	".\Debug\def_elmtyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\elmcar.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ELMCA=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ELMCA=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\elmchl.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ELMCH=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ELMCH=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\elmdel.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ELMDE=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ELMDE=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\elmder.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ELMDER=\
	".\Release\def_kintyp.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ELMDER=\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\elmgra.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ELMGR=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\mod_htable.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ELMGR=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\mod_htable.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\elmhes.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ELMHE=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ELMHE=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\elmle1.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ELMLE=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ELMLE=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\elmlen.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ELMLEN=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ELMLEN=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\elmmas.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ELMMA=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ELMMA=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\elmstr.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ELMST=\
	".\Release\def_elmtyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ELMST=\
	".\Debug\def_elmtyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\elmtyp.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ELMTY=\
	".\Release\def_domain.mod"\
	".\Release\def_elmtyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ELMTY=\
	".\Debug\def_domain.mod"\
	".\Debug\def_elmtyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\elsini.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ELSIN=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ELSIN=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\exampl.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_EXAMP=\
	".\Release\def_domain.mod"\
	".\Release\def_elmtyp.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_EXAMP=\
	".\Debug\def_domain.mod"\
	".\Debug\def_elmtyp.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\extbcs.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_EXTBC=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_EXTBC=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\extcog.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_EXTCO=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_EXTCO=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\extnor.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_EXTNO=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_EXTNO=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\finbou.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_FINBO=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_FINBO=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\fintyp.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_FINTY=\
	".\Release\def_elmtyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_FINTY=\
	".\Debug\def_elmtyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\geonor.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GEONO=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GEONO=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\gloloc.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GLOLO=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GLOLO=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\gpcabo.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GPCAB=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GPCAB=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\hangin.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_HANGI=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_HANGI=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\immers.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_IMMER=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_IMMER=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\inidom.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_INIDO=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_meshin.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_INIDO=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_meshin.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\massma.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_MASSM=\
	".\Release\def_domain.mod"\
	".\Release\def_elmtyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	".\Release\Mod_postpr.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_MASSM=\
	".\Debug\def_domain.mod"\
	".\Debug\def_elmtyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	".\Debug\Mod_postpr.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\materi.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_MATER=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_MATER=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\membcs.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_MEMBC=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_iofile.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_MEMBC=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_iofile.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\memgeo.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_MEMGEO=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_iofile.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_MEMGEO=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_iofile.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\memose.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_MEMOS=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_iofile.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_MEMOS=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_iofile.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\mergli.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_MERGL=\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_MERGL=\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\mescek.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_MESCE=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_MESCE=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\miniel.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_MINIE=\
	".\Release\def_domain.mod"\
	".\Release\def_elmtyp.mod"\
	".\Release\def_master.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_MINIE=\
	".\Debug\def_domain.mod"\
	".\Debug\def_elmtyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\nortri.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NORTR=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NORTR=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\nzecob.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NZECO=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NZECO=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\nzecof.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NZECOF=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NZECOF=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\opfdom.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_OPFDO=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_OPFDO=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\outdom.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_OUTDO=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_postpr.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_OUTDO=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_postpr.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\poscog.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_POSCO=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_POSCO=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\reabcs.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_REABC=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_REABC=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\reacod.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_REACO=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_REACO=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\readim.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_READI=\
	".\Release\def_domain.mod"\
	".\Release\def_elmtyp.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_READI=\
	".\Debug\def_domain.mod"\
	".\Debug\def_elmtyp.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\reageo.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_REAGE=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_iofile.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_REAGE=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_iofile.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\reaset.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_REASE=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_REASE=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\reasla.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_REASL=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_REASL=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\reastl.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_REAST=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_iofile.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_REAST=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_iofile.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\reastr.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_REASTR=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_iofile.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_REASTR=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_iofile.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\reatyp.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_REATY=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_iofile.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_REATY=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_iofile.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\rhsmod.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_RHSMO=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_RHSMO=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\rubclo.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_RUBCL=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_RUBCL=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\rubope.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_RUBOP=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_RUBOP=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\rulepw.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_RULEP=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_RULEP=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\rupclo.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_RUPCL=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_RUPCL=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\rupoin.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_RUPOI=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_RUPOI=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\rupope.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_RUPOP=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_RUPOP=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\ruqclo.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_RUQCL=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_RUQCL=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\ruqope.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_RUQOP=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_RUQOP=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\rutclo.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_RUTCL=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_RUTCL=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\rutope.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_RUTOP=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_RUTOP=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\ruyclo.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_RUYCL=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_RUYCL=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\ruyope.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_RUYOP=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_RUYOP=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\setext.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SETEX=\
	".\Release\def_domain.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SETEX=\
	".\Debug\def_domain.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\setlbe.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SETLB=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SETLB=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\shafga.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SHAFG=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SHAFG=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\shafun.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SHAFU=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SHAFU=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\shaga1.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SHAGA=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SHAGA=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\shaga2.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SHAGA2=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SHAGA2=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\shaga3.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SHAGA3=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SHAGA3=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\shape0.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SHAPE=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SHAPE=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\shape1.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SHAPE1=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SHAPE1=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\shape2.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SHAPE2=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SHAPE2=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\shape3.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SHAPE3=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SHAPE3=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\sslcon.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SSLCO=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SSLCO=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\sslcr2.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SSLCR=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SSLCR=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\sslcro.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SSLCRO=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SSLCRO=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\sslgct.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SSLGC=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SSLGC=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\sslmsh.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SSLMS=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SSLMS=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\sslord.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SSLOR=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SSLOR=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\sslpry.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SSLPR=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SSLPR=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\symgra.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SYMGR=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SYMGR=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\tanvec.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TANVE=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TANVE=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\updmsh.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_UPDMS=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_UPDMS=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\velchl.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_VELCH=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_VELCH=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\domain\witnes.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_WITNE=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_WITNE=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# End Group
# Begin Group "defmod"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\Sources\kernel\defmod\def_domain.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_DEF_D=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_DEF_D=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\defmod\def_elmope.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_DEF_E=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_DEF_E=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\defmod\def_elmtyp.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_DEF_EL=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_DEF_EL=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\defmod\def_elsest.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_DEF_ELS=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_DEF_ELS=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\defmod\def_inpout.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_DEF_I=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_DEF_I=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\defmod\def_kintyp.f90
# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\defmod\def_master.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_DEF_M=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_DEF_M=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\defmod\def_meshin.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_DEF_ME=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_DEF_ME=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\defmod\def_parame.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_DEF_P=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_DEF_P=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\defmod\def_postpr.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_DEF_PO=\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_DEF_PO=\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\defmod\def_solver.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_DEF_S=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_DEF_S=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# End Group
# Begin Group "master"

# PROP Default_Filter "f90"
# Begin Source File

SOURCE=..\..\Sources\kernel\master\algres.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ALGRE=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_solver.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ALGRE=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_solver.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\master\Alya.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ALYA_=\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ALYA_=\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\master\Begste.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_BEGST=\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_BEGST=\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\master\commun.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_COMMU=\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_COMMU=\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\master\Conblk.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CONBL=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CONBL=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\master\Concou.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CONCO=\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CONCO=\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\master\Doiter.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_DOITE=\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_DOITE=\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\master\Endste.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ENDST=\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ENDST=\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\master\iniomp.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_INIOM=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_INIOM=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\master\inirun.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_INIRU=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_meshin.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_INIRU=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_meshin.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\master\inisol.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_INISO=\
	".\Release\def_master.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_INISO=\
	".\Debug\def_master.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\master\iniste.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_INIST=\
	".\Release\def_master.mod"\
	".\Release\def_postpr.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_INIST=\
	".\Debug\def_master.mod"\
	".\Debug\def_postpr.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\master\inivar.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_INIVA=\
	".\Release\def_domain.mod"\
	".\Release\def_elmtyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_postpr.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_INIVA=\
	".\Debug\def_domain.mod"\
	".\Debug\def_elmtyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_postpr.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\master\mediso.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_MEDIS=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_MEDIS=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\master\memunk.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_MEMUN=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_MEMUN=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\master\mod_iofile.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_MOD_I=\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_MOD_I=\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\master\modser.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_MODSE=\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_MODSE=\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\master\moduls.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_MODUL=\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_MODUL=\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\master\Newmsh.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NEWMS=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NEWMS=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\master\openfi.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_OPENF=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_postpr.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_OPENF=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_postpr.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\master\readat.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_READA=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_READA=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\master\reamod.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_REAMO=\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_REAMO=\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\master\Reapro.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_REAPR=\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_REAPR=\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\master\restar.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_RESTA=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_RESTA=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\master\rrudat.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_RRUDA=\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_RRUDA=\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\master\runend.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_RUNEN=\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_RUNEN=\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\master\setgts.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SETGT=\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SETGT=\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\master\soldef.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SOLDE=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SOLDE=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\master\solver.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SOLVE=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_solver.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SOLVE=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_solver.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\master\timfun.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TIMFU=\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TIMFU=\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\master\Timste.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TIMST=\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TIMST=\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\master\Turnof.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TURNO=\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TURNO=\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\master\Turnon.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TURNON=\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TURNON=\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\master\vocabu.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_VOCAB=\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_VOCAB=\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# End Group
# Begin Group "elsest"

# PROP Default_Filter "f90"
# Begin Source File

SOURCE=..\..\Sources\kernel\elsest\elsest.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ELSES=\
	".\Release\def_elsest.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ELSES=\
	".\Debug\def_elsest.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\elsest\elsest_alloca.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ELSEST=\
	".\Release\def_elsest.mod"\
	".\Release\mod_elsest.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ELSEST=\
	".\Debug\def_elsest.mod"\
	".\Debug\mod_elsest.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\elsest\elsest_bindea.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ELSEST_=\
	".\Release\def_elsest.mod"\
	".\Release\mod_elsest.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ELSEST_=\
	".\Debug\def_elsest.mod"\
	".\Debug\mod_elsest.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\elsest\elsest_binlis.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ELSEST_B=\
	".\Release\def_elsest.mod"\
	".\Release\mod_elsest.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ELSEST_B=\
	".\Debug\def_elsest.mod"\
	".\Debug\mod_elsest.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\elsest\elsest_binpoi.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ELSEST_BI=\
	".\Release\def_elsest.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ELSEST_BI=\
	".\Debug\def_elsest.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\elsest\elsest_binpos.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ELSEST_BIN=\
	".\Release\def_elsest.mod"\
	".\Release\mod_elsest.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ELSEST_BIN=\
	".\Debug\def_elsest.mod"\
	".\Debug\mod_elsest.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\elsest\elsest_binpre.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ELSEST_BINP=\
	".\Release\def_elsest.mod"\
	".\Release\mod_elsest.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ELSEST_BINP=\
	".\Debug\def_elsest.mod"\
	".\Debug\mod_elsest.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\elsest\elsest_binpro.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ELSEST_BINPR=\
	".\Release\def_elsest.mod"\
	".\Release\mod_elsest.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ELSEST_BINPR=\
	".\Debug\def_elsest.mod"\
	".\Debug\mod_elsest.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\elsest\elsest_bintyp.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ELSEST_BINT=\
	".\Release\def_elsest.mod"\
	".\Release\mod_elsest.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ELSEST_BINT=\
	".\Debug\def_elsest.mod"\
	".\Debug\mod_elsest.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\elsest\elsest_boubox.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ELSEST_BO=\
	".\Release\def_elsest.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ELSEST_BO=\
	".\Debug\def_elsest.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\elsest\elsest_boxcoo.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ELSEST_BOX=\
	".\Release\def_elsest.mod"\
	".\Release\mod_elsest.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ELSEST_BOX=\
	".\Debug\def_elsest.mod"\
	".\Debug\mod_elsest.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\elsest\elsest_boxijk.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ELSEST_BOXI=\
	".\Release\def_elsest.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ELSEST_BOXI=\
	".\Debug\def_elsest.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\elsest\elsest_boxnum.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ELSEST_BOXN=\
	".\Release\def_elsest.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ELSEST_BOXN=\
	".\Debug\def_elsest.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\elsest\elsest_chkelm.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ELSEST_C=\
	".\Release\def_elsest.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ELSEST_C=\
	".\Debug\def_elsest.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\elsest\elsest_cputim.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ELSEST_CP=\
	".\Release\def_elsest.mod"\
	
NODEP_F90_ELSEST_CP=\
	".\Release\omp_lib.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ELSEST_CP=\
	".\Debug\def_elsest.mod"\
	
NODEP_F90_ELSEST_CP=\
	".\Debug\omp_lib.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\elsest\elsest_deallo.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ELSEST_D=\
	".\Release\def_elsest.mod"\
	".\Release\mod_elsest.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ELSEST_D=\
	".\Debug\def_elsest.mod"\
	".\Debug\mod_elsest.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\elsest\elsest_elq1p1.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ELSEST_E=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ELSEST_E=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\elsest\elsest_elq2p2.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ELSEST_EL=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ELSEST_EL=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\elsest\elsest_geogid.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ELSEST_G=\
	".\Release\def_elsest.mod"\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ELSEST_G=\
	".\Debug\def_elsest.mod"\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\elsest\elsest_inbyte.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ELSEST_I=\
	".\Release\def_elsest.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ELSEST_I=\
	".\Debug\def_elsest.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\elsest\elsest_invert.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ELSEST_IN=\
	".\Release\def_elsest.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ELSEST_IN=\
	".\Debug\def_elsest.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\elsest\elsest_invmtx.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ELSEST_INV=\
	".\Release\def_elsest.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ELSEST_INV=\
	".\Debug\def_elsest.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\elsest\elsest_memerr.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ELSEST_M=\
	".\Release\def_elsest.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ELSEST_M=\
	".\Debug\def_elsest.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\elsest\elsest_neighb.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ELSEST_N=\
	".\Release\def_elsest.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ELSEST_N=\
	".\Debug\def_elsest.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\elsest\elsest_newrap.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ELSEST_NE=\
	".\Release\def_elsest.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ELSEST_NE=\
	".\Debug\def_elsest.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\elsest\elsest_nextbo.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ELSEST_NEX=\
	".\Release\def_elsest.mod"\
	".\Release\mod_elsest.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ELSEST_NEX=\
	".\Debug\def_elsest.mod"\
	".\Debug\mod_elsest.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\elsest\elsest_octdea.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ELSEST_O=\
	".\Release\def_elsest.mod"\
	".\Release\mod_elsest.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ELSEST_O=\
	".\Debug\def_elsest.mod"\
	".\Debug\mod_elsest.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\elsest\elsest_octdep.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ELSEST_OC=\
	".\Release\def_elsest.mod"\
	".\Release\mod_elsest.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ELSEST_OC=\
	".\Debug\def_elsest.mod"\
	".\Debug\mod_elsest.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\elsest\elsest_octdes.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ELSEST_OCT=\
	".\Release\def_elsest.mod"\
	".\Release\mod_elsest.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ELSEST_OCT=\
	".\Debug\def_elsest.mod"\
	".\Debug\mod_elsest.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\elsest\elsest_octfin.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ELSEST_OCTF=\
	".\Release\def_elsest.mod"\
	".\Release\mod_elsest.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ELSEST_OCTF=\
	".\Debug\def_elsest.mod"\
	".\Debug\mod_elsest.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\elsest\elsest_octpoi.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ELSEST_OCTP=\
	".\Release\def_elsest.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ELSEST_OCTP=\
	".\Debug\def_elsest.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\elsest\elsest_octpos.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ELSEST_OCTPO=\
	".\Release\def_elsest.mod"\
	".\Release\mod_elsest.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ELSEST_OCTPO=\
	".\Debug\def_elsest.mod"\
	".\Debug\mod_elsest.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\elsest\elsest_octpre.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ELSEST_OCTPR=\
	".\Release\def_elsest.mod"\
	".\Release\mod_elsest.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ELSEST_OCTPR=\
	".\Debug\def_elsest.mod"\
	".\Debug\mod_elsest.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\elsest\elsest_octpro.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ELSEST_OCTPRO=\
	".\Release\def_elsest.mod"\
	".\Release\mod_elsest.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ELSEST_OCTPRO=\
	".\Debug\def_elsest.mod"\
	".\Debug\mod_elsest.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\elsest\elsest_octsub.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ELSEST_OCTS=\
	".\Release\def_elsest.mod"\
	".\Release\mod_elsest.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ELSEST_OCTS=\
	".\Debug\def_elsest.mod"\
	".\Debug\mod_elsest.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\elsest\elsest_recomm.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ELSEST_R=\
	".\Release\def_elsest.mod"\
	".\Release\mod_elsest.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ELSEST_R=\
	".\Debug\def_elsest.mod"\
	".\Debug\mod_elsest.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\elsest\elsest_resmat.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ELSEST_RE=\
	".\Release\def_elsest.mod"\
	".\Release\mod_elsest.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ELSEST_RE=\
	".\Debug\def_elsest.mod"\
	".\Debug\mod_elsest.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\elsest\elsest_resvec.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ELSEST_RES=\
	".\Release\def_elsest.mod"\
	".\Release\mod_elsest.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ELSEST_RES=\
	".\Debug\def_elsest.mod"\
	".\Debug\mod_elsest.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\elsest\elsest_runend.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ELSEST_RU=\
	".\Release\def_elsest.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ELSEST_RU=\
	".\Debug\def_elsest.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\elsest\elsest_segpro.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ELSEST_S=\
	".\Release\def_elsest.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ELSEST_S=\
	".\Debug\def_elsest.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\elsest\elsest_statis.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ELSEST_ST=\
	".\Release\def_elsest.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ELSEST_ST=\
	".\Debug\def_elsest.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\elsest\mod_elsest.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_MOD_E=\
	".\Release\def_elsest.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_MOD_E=\
	".\Debug\def_elsest.mod"\
	

!ENDIF 

# End Source File
# End Group
# Begin Group "meshin"

# PROP Default_Filter "f90"
# Begin Source File

SOURCE=..\..\Sources\kernel\meshin\cargeo.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CARGE=\
	".\Release\def_domain.mod"\
	".\Release\def_elmtyp.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_meshin.mod"\
	".\Release\mod_cart.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CARGE=\
	".\Debug\def_domain.mod"\
	".\Debug\def_elmtyp.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_meshin.mod"\
	".\Debug\mod_cart.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\meshin\fakept.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_FAKEP=\
	".\Release\def_inpout.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_meshin.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_FAKEP=\
	".\Debug\def_inpout.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_meshin.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\meshin\memmsh.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_MEMMS=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_meshin.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_iofile.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_MEMMS=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_meshin.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_iofile.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\meshin\meshin.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_MESHI=\
	".\Release\def_domain.mod"\
	".\Release\def_elmtyp.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_meshin.mod"\
	".\Release\def_parame.mod"\
	".\Release\mod_cart.mod"\
	".\Release\mod_extrpr.mod"\
	".\Release\Mod_memchk.mod"\
	".\Release\mod_source.mod"\
	".\Release\mod_surf.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_MESHI=\
	".\Debug\def_domain.mod"\
	".\Debug\def_elmtyp.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_meshin.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\mod_cart.mod"\
	".\Debug\mod_extrpr.mod"\
	".\Debug\Mod_memchk.mod"\
	".\Debug\mod_source.mod"\
	".\Debug\mod_surf.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\meshin\mod_cart.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_MOD_C=\
	".\Release\def_elmtyp.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_meshin.mod"\
	".\Release\Mod_memchk.mod"\
	".\Release\mod_mshtol.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_MOD_C=\
	".\Debug\def_elmtyp.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_meshin.mod"\
	".\Debug\Mod_memchk.mod"\
	".\Debug\mod_mshtol.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\meshin\mod_extrpr.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_MOD_EX=\
	".\Release\def_kintyp.mod"\
	".\Release\def_meshin.mod"\
	".\Release\Mod_memchk.mod"\
	".\Release\mod_mshtol.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_MOD_EX=\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_meshin.mod"\
	".\Debug\Mod_memchk.mod"\
	".\Debug\mod_mshtol.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\meshin\mod_interp.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_MOD_IN=\
	".\Release\def_kintyp.mod"\
	".\Release\def_meshin.mod"\
	".\Release\mod_cart.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_MOD_IN=\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_meshin.mod"\
	".\Debug\mod_cart.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\meshin\mod_mshtol.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_MOD_MS=\
	".\Release\def_kintyp.mod"\
	".\Release\def_meshin.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_MOD_MS=\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_meshin.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\meshin\mod_sort.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_MOD_S=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_MOD_S=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\meshin\mod_source.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_MOD_SO=\
	".\Release\def_kintyp.mod"\
	".\Release\def_meshin.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_MOD_SO=\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_meshin.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\meshin\mod_surf.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_MOD_SU=\
	".\Release\def_kintyp.mod"\
	".\Release\def_meshin.mod"\
	".\Release\mod_cart.mod"\
	".\Release\Mod_memchk.mod"\
	".\Release\mod_sort.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_MOD_SU=\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_meshin.mod"\
	".\Debug\mod_cart.mod"\
	".\Debug\Mod_memchk.mod"\
	".\Debug\mod_sort.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\meshin\read_geogrid.c
# End Source File
# Begin Source File

SOURCE=..\..\Sources\kernel\meshin\reamsh.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_REAMS=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_meshin.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_REAMS=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_meshin.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# End Group
# End Group
# Begin Group "modules"

# PROP Default_Filter "f90"
# Begin Group "temper"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\Sources\modules\temper\def_temper.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_DEF_T=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_DEF_T=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_assexp.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_A=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_A=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_averag.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_AV=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_temper.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_AV=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_temper.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_bcntoe.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_B=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_temper.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_B=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_temper.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_begite.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_BE=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_temper.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_BE=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_temper.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_begste.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_BEG=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_temper.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_BEG=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_temper.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_boumat.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_BO=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_BO=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_bounib.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_BOU=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_temper.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_BOU=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_temper.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_bouope.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_BOUO=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_temper.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_BOUO=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_temper.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_bouset.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_BOUS=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_temper.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_BOUS=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_temper.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_bouwal.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_BOUW=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_BOUW=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_conblk.f90
# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_concou.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_C=\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_temper.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_C=\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_temper.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_cvgunk.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_CV=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_temper.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_CV=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_temper.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_diffus.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_D=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_temper.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_D=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_temper.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_dodeme.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_DO=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_temper.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_DO=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_temper.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_doiter.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_DOI=\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	".\Release\def_temper.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_DOI=\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	".\Debug\def_temper.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_dync01.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_DY=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_temper.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_DY=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_temper.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_dyncou.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_DYN=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_temper.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_DYN=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_temper.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_elmadr.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_E=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_temper.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_E=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_temper.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_elmbub.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_EL=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_temper.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_EL=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_temper.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_elmchl.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_ELM=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_temper.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_ELM=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_temper.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_elmdir.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_ELMD=\
	".\Release\def_kintyp.mod"\
	".\Release\def_temper.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_ELMD=\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_temper.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_elmexa.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_ELME=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_temper.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_ELME=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_temper.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_elmgat.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_ELMG=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_temper.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_ELMG=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_temper.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_elmmat.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_ELMM=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_temper.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_ELMM=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_temper.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_elmop2.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_ELMO=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_temper.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_ELMO=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_temper.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_elmope.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_ELMOP=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_temper.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_ELMOP=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_temper.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_elmpre.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_ELMP=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_temper.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_ELMP=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_temper.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_elmpro.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_ELMPR=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_temper.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_ELMPR=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_temper.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_elmres.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_ELMR=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_temper.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_ELMR=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_temper.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_elmrhs.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_ELMRH=\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_temper.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_ELMRH=\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_temper.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_elmset.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_ELMS=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_temper.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_ELMS=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_temper.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_elmsgs.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_ELMSG=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_temper.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_ELMSG=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_temper.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_elmshc.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_ELMSH=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_temper.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_ELMSH=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_temper.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_elmtes.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_ELMT=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_temper.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_ELMT=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_temper.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_elmtss.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_ELMTS=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_temper.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_ELMTS=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_temper.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_endite.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_EN=\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_temper.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_EN=\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_temper.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_endste.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_END=\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_temper.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_END=\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_temper.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_exabcs.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_EX=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_temper.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_EX=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_temper.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_exacso.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_EXA=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_temper.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_EXA=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_temper.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_exaerr.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_EXAE=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_temper.mod"\
	".\Release\Mod_memchk.mod"\
	".\Release\Mod_postpr.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_EXAE=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_temper.mod"\
	".\Debug\Mod_memchk.mod"\
	".\Debug\Mod_postpr.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_exasol.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_EXAS=\
	".\Release\def_kintyp.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_EXAS=\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_heatfl.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_H=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_temper.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_H=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_temper.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_inicnd.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_I=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_temper.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_I=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_temper.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_inisol.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_IN=\
	".\Release\def_domain.mod"\
	".\Release\def_solver.mod"\
	".\Release\def_temper.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_IN=\
	".\Debug\def_domain.mod"\
	".\Debug\def_solver.mod"\
	".\Debug\def_temper.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_iniunk.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_INI=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_temper.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_INI=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_temper.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_inivar.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_INIV=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	".\Release\def_temper.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_INIV=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	".\Debug\def_temper.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_intbcs.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_INT=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_temper.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_INT=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_temper.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_intphy.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_INTP=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_temper.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_INTP=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_temper.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_matrix.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_M=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_temper.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_M=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_temper.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_memall.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_ME=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	".\Release\def_temper.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_ME=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	".\Debug\def_temper.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_membcs.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_MEM=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_temper.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_MEM=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_temper.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_membub.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_MEMB=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_temper.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_MEMB=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_temper.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_memose.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_MEMO=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_temper.mod"\
	".\Release\Mod_iofile.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_MEMO=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_temper.mod"\
	".\Debug\Mod_iofile.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_memphy.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_MEMP=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_temper.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_MEMP=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_temper.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_newmsh.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_N=\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_N=\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_openfi.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_O=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_temper.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_O=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_temper.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_outbcs.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_OU=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_temper.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_OU=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_temper.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_outerr.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_OUT=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_temper.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_OUT=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_temper.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_outhfl.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_OUTH=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_temper.mod"\
	".\Release\mod_gradie.mod"\
	".\Release\Mod_memchk.mod"\
	".\Release\Mod_postpr.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_OUTH=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_temper.mod"\
	".\Debug\mod_gradie.mod"\
	".\Debug\Mod_memchk.mod"\
	".\Debug\Mod_postpr.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_outinf.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_OUTI=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_temper.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_OUTI=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_temper.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_outlat.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_OUTL=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_temper.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_OUTL=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_temper.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_output.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_OUTP=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_temper.mod"\
	".\Release\Mod_iofile.mod"\
	".\Release\Mod_postpr.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_OUTP=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_temper.mod"\
	".\Debug\Mod_iofile.mod"\
	".\Debug\Mod_postpr.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_outset.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_OUTS=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_temper.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_OUTS=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_temper.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_outvar.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_OUTV=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_temper.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_OUTV=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_temper.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_outwit.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_OUTW=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_temper.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_OUTW=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_temper.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_parall.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_P=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_temper.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_P=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_temper.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_radpos.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_R=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_temper.mod"\
	".\Release\Mod_iofile.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_R=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_temper.mod"\
	".\Debug\Mod_iofile.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_radvuf.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_RA=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_temper.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_RA=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_temper.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_reabcs.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_RE=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_temper.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_RE=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_temper.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_reanut.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_REA=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	".\Release\def_temper.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_REA=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	".\Debug\def_temper.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_reaous.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_REAO=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_temper.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_REAO=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_temper.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_reaphy.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_REAP=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_temper.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_REAP=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_temper.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_reapro.f90
# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_restar.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_RES=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_temper.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_RES=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_temper.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_sendat.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_S=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	".\Release\def_temper.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_S=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	".\Debug\def_temper.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_sgsope.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_SG=\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_temper.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_SG=\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_temper.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_solite.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_SO=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_temper.mod"\
	".\Release\mod_gradie.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_SO=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_temper.mod"\
	".\Debug\mod_gradie.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_timste.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_T=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_temper.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_T=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_temper.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_tistep.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_TI=\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_temper.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_TI=\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_temper.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_turnof.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_TU=\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	".\Release\def_temper.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_TU=\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	".\Debug\def_temper.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_turnon.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_TUR=\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_temper.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_TUR=\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_temper.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_updbcs.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_U=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_temper.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_U=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_temper.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_updhfl.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_UP=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_temper.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_UP=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_temper.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_updtss.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_UPD=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_temper.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_UPD=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_temper.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_updunk.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_UPDU=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_temper.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_UPDU=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_temper.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\tem_velfun.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEM_V=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_temper.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEM_V=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_temper.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\temper\Temper.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TEMPE=\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TEMPE=\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# End Group
# Begin Group "nastin"

# PROP Default_Filter "f90"
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\def_nastin.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_DEF_N=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_DEF_N=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\mod_nsi_solsch.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_MOD_N=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_MOD_N=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\Nastin.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NASTI=\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NASTI=\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_addcst.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_A=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_nastin.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_A=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_nastin.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_assdia.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_AS=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_AS=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_assma3.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_ASS=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_nastin.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_ASS=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_nastin.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_assmat.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_ASSM=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_nastin.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_ASSM=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_nastin.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_assrhs.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_ASSR=\
	".\Release\def_kintyp.mod"\
	".\Release\def_nastin.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_ASSR=\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_nastin.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_autbcs.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_AU=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_AU=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_begite.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_B=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_B=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_begste.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_BE=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_BE=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_bibset.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_BI=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_BI=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_bouass.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_BO=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_BO=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_bouave.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_BOU=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_BOU=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_bougat.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_BOUG=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_BOUG=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_bougau.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_BOUGA=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_BOUGA=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_bounib.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_BOUN=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_BOUN=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_bounod.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_BOUNO=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_BOUNO=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_bouopb.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_BOUO=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_nastin.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_BOUO=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_nastin.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_bouope.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_BOUOP=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_BOUOP=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_bouopp.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_BOUOPP=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_BOUOPP=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_boupre.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_BOUP=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_BOUP=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_bouset.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_BOUS=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_BOUS=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_bouwal.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_BOUW=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_nastin.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_BOUW=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_nastin.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_chkbcs.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_C=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_C=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_conblk.f90
# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_concou.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_CO=\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_CO=\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_cvgunk.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_CV=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_parame.mod"\
	".\Release\mod_nsi_solsch.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_CV=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\mod_nsi_solsch.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_denvis.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_D=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_D=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_doiter.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_DO=\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_DO=\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_dommas.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_DOM=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_DOM=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_dync01.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_DY=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_DY=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_dyncou.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_DYN=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_DYN=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_elmco2.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_E=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_nastin.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_E=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_nastin.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_elmcof.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_EL=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_nastin.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_EL=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_nastin.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_elmcon.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_ELM=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_nastin.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_ELM=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_nastin.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_elmcor.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_ELMC=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_parame.mod"\
	".\Release\mod_gradie.mod"\
	".\Release\mod_nsi_solsch.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_ELMC=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\mod_gradie.mod"\
	".\Debug\mod_nsi_solsch.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_elmcp1.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_ELMCP=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_nastin.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_ELMCP=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_nastin.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_elmcst.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_ELMCS=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_ELMCS=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_elmdi2.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_ELMD=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_nastin.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_ELMD=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_nastin.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_elmdi3.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_ELMDI=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_nastin.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_ELMDI=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_nastin.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_elmdir.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_ELMDIR=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_nastin.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_ELMDIR=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_nastin.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_elmexa.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_ELME=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_nastin.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_ELME=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_nastin.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_elmga2.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_ELMG=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_ELMG=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_elmga3.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_ELMGA=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_ELMGA=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_elmgat.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_ELMGAT=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_ELMGAT=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_elmibm.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_ELMI=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_ELMI=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_elmini.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_ELMIN=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_ELMIN=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_elmlap.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_ELML=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_nastin.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_ELML=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_nastin.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_elmlbc.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_ELMLB=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_nastin.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_ELMLB=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_nastin.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_elmma2.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_ELMM=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_nastin.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_ELMM=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_nastin.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_elmma3.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_ELMMA=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_nastin.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_ELMMA=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_nastin.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_elmmat.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_ELMMAT=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_ELMMAT=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_elmmo2.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_ELMMO=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_nastin.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_ELMMO=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_nastin.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_elmmoc.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_ELMMOC=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_nastin.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_ELMMOC=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_nastin.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_elmmof.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_ELMMOF=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_nastin.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_ELMMOF=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_nastin.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_elmmom.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_ELMMOM=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_nastin.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_ELMMOM=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_nastin.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_elmmp1.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_ELMMP=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_nastin.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_ELMMP=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_nastin.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_elmop2.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_ELMO=\
	".\Release\def_domain.mod"\
	".\Release\def_elmtyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_ELMO=\
	".\Debug\def_domain.mod"\
	".\Debug\def_elmtyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_elmop3.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_ELMOP=\
	".\Release\def_domain.mod"\
	".\Release\def_elmtyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_ELMOP=\
	".\Debug\def_domain.mod"\
	".\Debug\def_elmtyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_elmope.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_ELMOPE=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_ELMOPE=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_elmp13.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_ELMP=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_nastin.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_ELMP=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_nastin.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_elmpr2.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_ELMPR=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_ELMPR=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_elmprc.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_ELMPRC=\
	".\Release\def_kintyp.mod"\
	".\Release\def_nastin.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_ELMPRC=\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_nastin.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_elmpre.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_ELMPRE=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_nastin.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_ELMPRE=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_nastin.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_elmprl.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_ELMPRL=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_nastin.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_ELMPRL=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_nastin.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_elmpro.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_ELMPRO=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_ELMPRO=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_elmre2.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_ELMR=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_nastin.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_ELMR=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_nastin.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_elmre3.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_ELMRE=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_nastin.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_ELMRE=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_nastin.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_elmres.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_ELMRES=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_nastin.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_ELMRES=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_nastin.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_elmrh2.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_ELMRH=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_nastin.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_ELMRH=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_nastin.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_elmrhc.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_ELMRHC=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_ELMRHC=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_elmrhs.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_ELMRHS=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	
# ADD F90 /fpp
# SUBTRACT F90 /nodefine

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_ELMRHS=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_elmrp1.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_ELMRP=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_nastin.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_ELMRP=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_nastin.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_elmsch.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_ELMS=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_nastin.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_ELMS=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_nastin.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_elmset.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_ELMSE=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_ELMSE=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_elmsg2.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_ELMSG=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_nastin.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_ELMSG=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_nastin.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_elmsgs.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_ELMSGS=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_nastin.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_ELMSGS=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_nastin.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_elmshc.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_ELMSH=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_nastin.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_ELMSH=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_nastin.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_elmte2.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_ELMT=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_nastin.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_ELMT=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_nastin.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_elmtes.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_ELMTE=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_nastin.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_ELMTE=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_nastin.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_elmtss.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_ELMTS=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_ELMTS=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_endite.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_EN=\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_EN=\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_endste.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_END=\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_END=\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_esspre.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_ES=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_ES=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_exacso.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_EX=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_EX=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_exaerr.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_EXA=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_EXA=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_ifconf.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_I=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_I=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_inisol.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_IN=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_IN=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_iniunk.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_INI=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_INI=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_inivar.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_INIV=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_INIV=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_intbcs.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_INT=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_INT=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_intrst.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_INTR=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_iofile.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_INTR=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_iofile.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_levels.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_L=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_L=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_matrix.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_M=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_M=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_memall.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_ME=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_ME=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_membcs.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_MEM=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_MEM=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_memose.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_MEMO=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_iofile.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_MEMO=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_iofile.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_memphy.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_MEMP=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_MEMP=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_minmax.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_MI=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_MI=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_newmsh.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_N=\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_N=\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_oldmat.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_O=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_O=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_openda.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_OP=\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_OP=\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_openfi.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_OPE=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_postpr.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_OPE=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_postpr.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_opesch.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_OPES=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_OPES=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_outbcs.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_OU=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_OU=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_outcpu.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_OUT=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_OUT=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_outerr.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_OUTE=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_OUTE=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_outinf.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_OUTI=\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_OUTI=\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_outlat.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_OUTL=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_OUTL=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_outpmv.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_OUTP=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_OUTP=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_output.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_OUTPU=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_postpr.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_OUTPU=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_postpr.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_outres.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_OUTR=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_OUTR=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_outset.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_OUTS=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_OUTS=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_outtan.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_OUTT=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_parame.mod"\
	".\Release\mod_gradie.mod"\
	".\Release\Mod_memchk.mod"\
	".\Release\Mod_postpr.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_OUTT=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\mod_gradie.mod"\
	".\Debug\Mod_memchk.mod"\
	".\Debug\Mod_postpr.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_outtau.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_OUTTA=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_OUTTA=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_outvar.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_OUTV=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_OUTV=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_parall.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_P=\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_P=\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_pmvppd.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_PM=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_PM=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_potent.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_PO=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_solver.mod"\
	".\Release\mod_gradie.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_PO=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_solver.mod"\
	".\Debug\mod_gradie.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_precon.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_PR=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_solver.mod"\
	".\Release\Mod_postpr.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_PR=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_solver.mod"\
	".\Debug\Mod_postpr.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_prelev.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_PRE=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_PRE=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_promsh.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_PRO=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_PRO=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_proper.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_PROP=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_PROP=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_reabcp.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_R=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_R=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_reabcs.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_RE=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_RE=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_reanut.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_REA=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_REA=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_reaous.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_REAO=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_REAO=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_reaphy.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_REAP=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_REAP=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_reapro.f90
# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_refere.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_REF=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_REF=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_rescon.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_RES=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_RES=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_residu.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_RESI=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_RESI=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_restar.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_REST=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_REST=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_rotma3.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_RO=\
	".\Release\def_domain.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_RO=\
	".\Debug\def_domain.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_rotmat.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_ROT=\
	".\Release\def_domain.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_ROT=\
	".\Debug\def_domain.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_rotunk.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_ROTU=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_ROTU=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_schpre.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_S=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_solver.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_S=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_solver.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_sendat.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_SE=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_SE=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_solinp.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_SO=\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_SO=\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_solite.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_SOL=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_parame.mod"\
	".\Release\mod_nsi_solsch.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_SOL=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\mod_nsi_solsch.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_solmon.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_SOLM=\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_SOLM=\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_timste.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_T=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_T=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_tistep.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_TI=\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_TI=\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_turnof.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_TU=\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_TU=\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_turnon.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_TUR=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_TUR=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_updbcs.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_U=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_U=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_updfor.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_UP=\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_UP=\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_updmsh.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_UPD=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_UPD=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_updrel.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_UPDR=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_UPDR=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_updthe.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_UPDT=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_UPDT=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_updtss.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_UPDTS=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_UPDTS=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_updunk.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_UPDU=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_UPDU=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastin\nsi_usrbcs.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSI_US=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSI_US=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# End Group
# Begin Group "codire"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\Sources\modules\codire\cdr_assmat.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CDR_A=\
	".\Release\def_kintyp.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CDR_A=\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\cdr_assrhs.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CDR_AS=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CDR_AS=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\cdr_begite.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CDR_B=\
	".\Release\def_codire.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CDR_B=\
	".\Debug\def_codire.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\cdr_begste.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CDR_BE=\
	".\Release\def_codire.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CDR_BE=\
	".\Debug\def_codire.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\cdr_bouope.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CDR_BO=\
	".\Release\def_cdrloc.mod"\
	".\Release\def_codire.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CDR_BO=\
	".\Debug\def_cdrloc.mod"\
	".\Debug\def_codire.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\cdr_conblk.f90
# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\cdr_concou.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CDR_C=\
	".\Release\def_codire.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CDR_C=\
	".\Debug\def_codire.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\cdr_cvgsca.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CDR_CV=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CDR_CV=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\cdr_cvgunk.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CDR_CVG=\
	".\Release\def_codire.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CDR_CVG=\
	".\Debug\def_codire.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\cdr_cvgvec.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CDR_CVGV=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CDR_CVGV=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\cdr_defbou.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CDR_D=\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CDR_D=\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\cdr_defdif.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CDR_DE=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CDR_DE=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\cdr_defgen.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CDR_DEF=\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CDR_DEF=\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\cdr_deflom.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CDR_DEFL=\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CDR_DEFL=\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\cdr_defmat.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CDR_DEFM=\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CDR_DEFM=\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\cdr_defnst.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CDR_DEFN=\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CDR_DEFN=\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\cdr_defpla.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CDR_DEFP=\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CDR_DEFP=\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\cdr_denvel.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CDR_DEN=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CDR_DEN=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\cdr_doiter.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CDR_DO=\
	".\Release\def_codire.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CDR_DO=\
	".\Debug\def_codire.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\cdr_elmdir.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CDR_E=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CDR_E=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\cdr_elmgal.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CDR_EL=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CDR_EL=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\cdr_elmgat.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CDR_ELM=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CDR_ELM=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\cdr_elmmul.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CDR_ELMM=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CDR_ELMM=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\cdr_elmope.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CDR_ELMO=\
	".\Release\def_cdrloc.mod"\
	".\Release\def_codire.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CDR_ELMO=\
	".\Debug\def_cdrloc.mod"\
	".\Debug\def_codire.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\cdr_elmres.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CDR_ELMR=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CDR_ELMR=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\cdr_elmtes.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CDR_ELMT=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CDR_ELMT=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\cdr_endite.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CDR_EN=\
	".\Release\def_codire.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CDR_EN=\
	".\Debug\def_codire.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\cdr_endste.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CDR_END=\
	".\Release\def_codire.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CDR_END=\
	".\Debug\def_codire.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\cdr_forces.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CDR_F=\
	".\Release\def_cdrloc.mod"\
	".\Release\def_codire.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_iofile.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CDR_F=\
	".\Debug\def_cdrloc.mod"\
	".\Debug\def_codire.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_iofile.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\cdr_iniofl.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CDR_I=\
	".\Release\def_codire.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CDR_I=\
	".\Debug\def_codire.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\cdr_inisol.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CDR_IN=\
	".\Release\def_codire.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CDR_IN=\
	".\Debug\def_codire.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\cdr_initpr.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CDR_INI=\
	".\Release\def_codire.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CDR_INI=\
	".\Debug\def_codire.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\cdr_iniunk.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CDR_INIU=\
	".\Release\def_codire.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CDR_INIU=\
	".\Debug\def_codire.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\cdr_linsea.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CDR_L=\
	".\Release\def_codire.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CDR_L=\
	".\Debug\def_codire.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\cdr_lsarmj.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CDR_LS=\
	".\Release\def_codire.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CDR_LS=\
	".\Debug\def_codire.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\cdr_lspics.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CDR_LSP=\
	".\Release\def_codire.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CDR_LSP=\
	".\Debug\def_codire.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\cdr_lstarg.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CDR_LST=\
	".\Release\def_cdrloc.mod"\
	".\Release\def_codire.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CDR_LST=\
	".\Debug\def_cdrloc.mod"\
	".\Debug\def_codire.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\cdr_matrix.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CDR_M=\
	".\Release\def_cdrloc.mod"\
	".\Release\def_codire.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CDR_M=\
	".\Debug\def_cdrloc.mod"\
	".\Debug\def_codire.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\cdr_memall.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CDR_ME=\
	".\Release\def_codire.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	".\Release\Mod_iofile.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CDR_ME=\
	".\Debug\def_codire.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	".\Debug\Mod_iofile.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\cdr_memass.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CDR_MEM=\
	".\Release\def_cdrloc.mod"\
	".\Release\def_codire.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CDR_MEM=\
	".\Debug\def_cdrloc.mod"\
	".\Debug\def_codire.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\cdr_membou.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CDR_MEMB=\
	".\Release\def_cdrloc.mod"\
	".\Release\def_codire.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CDR_MEMB=\
	".\Debug\def_cdrloc.mod"\
	".\Debug\def_codire.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\cdr_memelm.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CDR_MEME=\
	".\Release\def_cdrloc.mod"\
	".\Release\def_codire.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CDR_MEME=\
	".\Debug\def_cdrloc.mod"\
	".\Debug\def_codire.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\cdr_newmsh.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CDR_N=\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CDR_N=\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\cdr_openfi.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CDR_O=\
	".\Release\def_codire.mod"\
	".\Release\def_master.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CDR_O=\
	".\Debug\def_codire.mod"\
	".\Debug\def_master.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\cdr_outerr.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CDR_OU=\
	".\Release\def_codire.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CDR_OU=\
	".\Debug\def_codire.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\cdr_outinf.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CDR_OUT=\
	".\Release\def_codire.mod"\
	".\Release\def_master.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CDR_OUT=\
	".\Debug\def_codire.mod"\
	".\Debug\def_master.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\cdr_output.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CDR_OUTP=\
	".\Release\def_codire.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_iofile.mod"\
	".\Release\Mod_postpr.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CDR_OUTP=\
	".\Debug\def_codire.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_iofile.mod"\
	".\Debug\Mod_postpr.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\cdr_outsec.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CDR_OUTS=\
	".\Release\def_codire.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CDR_OUTS=\
	".\Debug\def_codire.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\cdr_outtpo.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CDR_OUTT=\
	".\Release\def_codire.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CDR_OUTT=\
	".\Debug\def_codire.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\cdr_reabcs.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CDR_R=\
	".\Release\def_codire.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CDR_R=\
	".\Debug\def_codire.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\cdr_reanut.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CDR_RE=\
	".\Release\def_codire.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CDR_RE=\
	".\Debug\def_codire.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\cdr_reaous.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CDR_REA=\
	".\Release\def_codire.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CDR_REA=\
	".\Debug\def_codire.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\cdr_reaphy.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CDR_REAP=\
	".\Release\def_codire.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_iofile.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CDR_REAP=\
	".\Debug\def_codire.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_iofile.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\cdr_reapro.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CDR_REAPR=\
	".\Release\def_codire.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CDR_REAPR=\
	".\Debug\def_codire.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\cdr_restar.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CDR_RES=\
	".\Release\def_codire.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CDR_RES=\
	".\Debug\def_codire.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\cdr_shflow.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CDR_S=\
	".\Release\def_cdrloc.mod"\
	".\Release\def_codire.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CDR_S=\
	".\Debug\def_cdrloc.mod"\
	".\Debug\def_codire.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\cdr_smflow.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CDR_SM=\
	".\Release\def_cdrloc.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CDR_SM=\
	".\Debug\def_cdrloc.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\cdr_solite.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CDR_SO=\
	".\Release\def_codire.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CDR_SO=\
	".\Debug\def_codire.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\cdr_timste.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CDR_T=\
	".\Release\def_codire.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CDR_T=\
	".\Debug\def_codire.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\cdr_tistep.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CDR_TI=\
	".\Release\def_codire.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CDR_TI=\
	".\Debug\def_codire.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\cdr_turnof.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CDR_TU=\
	".\Release\def_codire.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CDR_TU=\
	".\Debug\def_codire.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\cdr_turnon.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CDR_TUR=\
	".\Release\def_codire.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CDR_TUR=\
	".\Debug\def_codire.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\cdr_updbcs.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CDR_U=\
	".\Release\def_codire.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CDR_U=\
	".\Debug\def_codire.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\cdr_updtpr.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CDR_UP=\
	".\Release\def_codire.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CDR_UP=\
	".\Debug\def_codire.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\cdr_updunk.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CDR_UPD=\
	".\Release\def_codire.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CDR_UPD=\
	".\Debug\def_codire.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\cdr_vinvte.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CDR_V=\
	".\Release\def_cdrloc.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CDR_V=\
	".\Debug\def_cdrloc.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\cdr_vislaw.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CDR_VI=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CDR_VI=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\Codire.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CODIR=\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CODIR=\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\def_cdrloc.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_DEF_C=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_DEF_C=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\codire\def_codire.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_DEF_CO=\
	"..\..\..\Mumps\include\dmumps_root.h"\
	"..\..\..\Mumps\include\dmumps_struc.h"\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_DEF_CO=\
	"..\..\..\Mumps\include\dmumps_root.h"\
	"..\..\..\Mumps\include\dmumps_struc.h"\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# End Group
# Begin Group "alefor"

# PROP Default_Filter "f90"
# Begin Source File

SOURCE=..\..\Sources\modules\alefor\ale_begite.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ALE_B=\
	".\Release\def_alefor.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ALE_B=\
	".\Debug\def_alefor.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\alefor\ale_begste.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ALE_BE=\
	".\Release\def_alefor.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ALE_BE=\
	".\Debug\def_alefor.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\alefor\ale_concou.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ALE_C=\
	".\Release\def_alefor.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ALE_C=\
	".\Debug\def_alefor.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\alefor\ale_doiter.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ALE_D=\
	".\Release\def_alefor.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ALE_D=\
	".\Debug\def_alefor.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\alefor\ale_endste.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ALE_E=\
	".\Release\def_alefor.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ALE_E=\
	".\Debug\def_alefor.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\alefor\ale_inirot.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ALE_I=\
	".\Release\def_alefor.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ALE_I=\
	".\Debug\def_alefor.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\alefor\ale_iniunk.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ALE_IN=\
	".\Release\def_alefor.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ALE_IN=\
	".\Debug\def_alefor.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\alefor\ale_memall.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ALE_M=\
	".\Release\def_alefor.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ALE_M=\
	".\Debug\def_alefor.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\alefor\ale_openfi.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ALE_O=\
	".\Release\def_alefor.mod"\
	".\Release\def_master.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ALE_O=\
	".\Debug\def_alefor.mod"\
	".\Debug\def_master.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\alefor\ale_outinf.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ALE_OU=\
	".\Release\def_alefor.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ALE_OU=\
	".\Debug\def_alefor.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\alefor\ale_output.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ALE_OUT=\
	".\Release\def_alefor.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_iofile.mod"\
	".\Release\Mod_postpr.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ALE_OUT=\
	".\Debug\def_alefor.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_iofile.mod"\
	".\Debug\Mod_postpr.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\alefor\ale_reabcs.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ALE_R=\
	".\Release\def_alefor.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ALE_R=\
	".\Debug\def_alefor.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\alefor\ale_reaous.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ALE_RE=\
	".\Release\def_alefor.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ALE_RE=\
	".\Debug\def_alefor.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\alefor\ale_reaphy.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ALE_REA=\
	".\Release\def_alefor.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ALE_REA=\
	".\Debug\def_alefor.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\alefor\ale_reapro.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ALE_REAP=\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ALE_REAP=\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\alefor\ale_rotate.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ALE_RO=\
	".\Release\def_alefor.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastin.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ALE_RO=\
	".\Debug\def_alefor.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastin.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\alefor\ale_turnof.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ALE_T=\
	".\Release\def_alefor.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ALE_T=\
	".\Debug\def_alefor.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\alefor\ale_turnon.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ALE_TU=\
	".\Release\def_alefor.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ALE_TU=\
	".\Debug\def_alefor.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\alefor\Alefor.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_ALEFO=\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_ALEFO=\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\alefor\def_alefor.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_DEF_A=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_DEF_A=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# End Group
# Begin Group "turbul"

# PROP Default_Filter "f90"
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\def_turbul.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_DEF_TU=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_DEF_TU=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_adapti.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_A=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_turbul.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_A=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_turbul.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_addarr.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_AD=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_turbul.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_AD=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_turbul.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_assexp.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_AS=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_AS=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_assmat.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_ASS=\
	".\Release\def_kintyp.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_ASS=\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_assrhs.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_ASSR=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_ASSR=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_begite.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_B=\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_turbul.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_B=\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_turbul.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_begste.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_BE=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_turbul.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_BE=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_turbul.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_bouave.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_BO=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_turbul.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_BO=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_turbul.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_brepen.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_BR=\
	".\Release\def_kintyp.mod"\
	".\Release\def_turbul.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_BR=\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_turbul.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_clippi.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_C=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_turbul.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_C=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_turbul.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_conblk.f90
# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_concou.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_CO=\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_turbul.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_CO=\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_turbul.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_cvgunk.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_CV=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_turbul.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_CV=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_turbul.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_doiter.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_D=\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	".\Release\def_turbul.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_D=\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	".\Debug\def_turbul.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_elmadj.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_E=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_turbul.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_E=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_turbul.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_elmal3.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_EL=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_turbul.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_EL=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_turbul.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_elmco2.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_ELM=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_turbul.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_ELM=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_turbul.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_elmcoe.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_ELMC=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_turbul.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_ELMC=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_turbul.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_elmdif.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_ELMD=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_turbul.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_ELMD=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_turbul.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_elmdir.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_ELMDI=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_turbul.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_ELMDI=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_turbul.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_elmgat.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_ELMG=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_turbul.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_ELMG=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_turbul.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_elmla2.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_ELML=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_turbul.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_ELML=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_turbul.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_elmlap.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_ELMLA=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_turbul.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_ELMLA=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_turbul.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_elmma2.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_ELMM=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_turbul.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_ELMM=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_turbul.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_elmmat.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_ELMMA=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_turbul.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_ELMMA=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_turbul.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_elmop2.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_ELMO=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_turbul.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_ELMO=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_turbul.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_elmope.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_ELMOP=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_turbul.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_ELMOP=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_turbul.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_elmpro.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_ELMP=\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_turbul.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_ELMP=\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_turbul.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_elmshc.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_ELMS=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_turbul.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_ELMS=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_turbul.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_elmsta.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_ELMST=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_turbul.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_ELMST=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_turbul.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_elmtss.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_ELMT=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_turbul.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_ELMT=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_turbul.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_elmust.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_ELMU=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_turbul.mod"\
	".\Release\mod_gradie.mod"\
	".\Release\Mod_memchk.mod"\
	".\Release\Mod_postpr.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_ELMU=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_turbul.mod"\
	".\Debug\mod_gradie.mod"\
	".\Debug\Mod_memchk.mod"\
	".\Debug\Mod_postpr.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_endite.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_EN=\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_turbul.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_EN=\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_turbul.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_endste.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_END=\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_turbul.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_END=\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_turbul.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_frivel.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_F=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_turbul.mod"\
	".\Release\mod_gradie.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_F=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_turbul.mod"\
	".\Debug\mod_gradie.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_grave2.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_G=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_turbul.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_G=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_turbul.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_grsqki.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_GR=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_turbul.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_GR=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_turbul.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_inisol.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_I=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_solver.mod"\
	".\Release\def_turbul.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_I=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_solver.mod"\
	".\Debug\def_turbul.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_iniunk.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_IN=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_turbul.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_IN=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_turbul.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_inivar.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_INI=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	".\Release\def_turbul.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_INI=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	".\Debug\def_turbul.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_kchien.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_K=\
	".\Release\def_kintyp.mod"\
	".\Release\def_turbul.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_K=\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_turbul.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_kepsil.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_KE=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_turbul.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_KE=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_turbul.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_kepspf.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_KEP=\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_turbul.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_KEP=\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_turbul.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_kepsv2.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_KEPS=\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_turbul.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_KEPS=\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_turbul.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_kjawhw.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_KJ=\
	".\Release\def_kintyp.mod"\
	".\Release\def_turbul.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_KJ=\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_turbul.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_komega.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_KO=\
	".\Release\def_kintyp.mod"\
	".\Release\def_turbul.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_KO=\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_turbul.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_lambre.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_L=\
	".\Release\def_kintyp.mod"\
	".\Release\def_turbul.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_L=\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_turbul.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_lausha.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_LA=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_turbul.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_LA=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_turbul.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_matrix.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_M=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_turbul.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_M=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_turbul.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_memall.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_ME=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	".\Release\def_turbul.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_ME=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	".\Debug\def_turbul.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_memarr.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_MEM=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_turbul.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_MEM=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_turbul.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_membcs.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_MEMB=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_turbul.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_MEMB=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_turbul.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_nagano.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_N=\
	".\Release\def_kintyp.mod"\
	".\Release\def_turbul.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_N=\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_turbul.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_newmsh.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_NE=\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_NE=\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_nut2nd.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_NU=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_turbul.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_NU=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_turbul.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_openfi.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_O=\
	".\Release\def_master.mod"\
	".\Release\def_turbul.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_O=\
	".\Debug\def_master.mod"\
	".\Debug\def_turbul.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_outerr.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_OU=\
	".\Release\def_master.mod"\
	".\Release\def_turbul.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_OU=\
	".\Debug\def_master.mod"\
	".\Debug\def_turbul.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_outinf.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_OUT=\
	".\Release\def_master.mod"\
	".\Release\def_turbul.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_OUT=\
	".\Debug\def_master.mod"\
	".\Debug\def_turbul.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_output.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_OUTP=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_turbul.mod"\
	".\Release\Mod_iofile.mod"\
	".\Release\Mod_postpr.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_OUTP=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_turbul.mod"\
	".\Debug\Mod_iofile.mod"\
	".\Debug\Mod_postpr.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_outvar.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_OUTV=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_turbul.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_OUTV=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_turbul.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_parall.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_P=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_P=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_projec.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_PR=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_turbul.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_PR=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_turbul.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_reabcs.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_R=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_turbul.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_R=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_turbul.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_reanut.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_RE=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	".\Release\def_turbul.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_RE=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	".\Debug\def_turbul.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_reaous.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_REA=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_turbul.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_REA=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_turbul.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_reaphy.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_REAP=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_turbul.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_REAP=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_turbul.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_reapro.f90
# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_restar.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_RES=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_turbul.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_RES=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_turbul.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_secvar.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_S=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_turbul.mod"\
	".\Release\Mod_iofile.mod"\
	".\Release\Mod_memchk.mod"\
	".\Release\Mod_postpr.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_S=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_turbul.mod"\
	".\Debug\Mod_iofile.mod"\
	".\Debug\Mod_memchk.mod"\
	".\Debug\Mod_postpr.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_sendat.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_SE=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	".\Release\def_turbul.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_SE=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	".\Debug\def_turbul.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_solexp.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_SO=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_turbul.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_SO=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_turbul.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_solite.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_SOL=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_turbul.mod"\
	".\Release\mod_gradie.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_SOL=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_turbul.mod"\
	".\Debug\mod_gradie.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_spaalm.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_SP=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_SP=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_spanut.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_SPA=\
	".\Release\def_kintyp.mod"\
	".\Release\def_turbul.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_SPA=\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_turbul.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_stdkep.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_ST=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_turbul.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_ST=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_turbul.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_timste.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_T=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_turbul.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_T=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_turbul.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_tistep.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_TI=\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_turbul.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_TI=\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_turbul.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_turnof.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_TU=\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	".\Release\def_turbul.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_TU=\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	".\Debug\def_turbul.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_turnon.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_TUR=\
	".\Release\def_master.mod"\
	".\Release\def_turbul.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_TUR=\
	".\Debug\def_master.mod"\
	".\Debug\def_turbul.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_twonat.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_TW=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_TW=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_tworod.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_TWO=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_TWO=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_updbcs.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_U=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_turbul.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_U=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_turbul.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_updedd.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_UP=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_turbul.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_UP=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_turbul.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_updrel.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_UPD=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_turbul.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_UPD=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_turbul.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_updtss.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_UPDT=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_turbul.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_UPDT=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_turbul.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_updunk.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_UPDU=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_turbul.mod"\
	".\Release\Mod_postpr.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_UPDU=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_turbul.mod"\
	".\Debug\Mod_postpr.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_usrbcs.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_US=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_turbul.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_US=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_turbul.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_walgen.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_W=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	".\Release\def_turbul.mod"\
	".\Release\mod_gradie.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_W=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	".\Debug\def_turbul.mod"\
	".\Debug\mod_gradie.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\tur_xuchen.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TUR_X=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TUR_X=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\turbul\Turbul.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_TURBU=\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_TURBU=\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# End Group
# Begin Group "exmedi"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\Sources\modules\exmedi\def_exmedi.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_DEF_EX=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_DEF_EX=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\exmedi\exm_assmat.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_EXM_A=\
	".\Release\def_kintyp.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_EXM_A=\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\exmedi\exm_begite.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_EXM_B=\
	".\Release\def_domain.mod"\
	".\Release\def_exmedi.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_EXM_B=\
	".\Debug\def_domain.mod"\
	".\Debug\def_exmedi.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\exmedi\exm_begste.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_EXM_BE=\
	".\Release\def_domain.mod"\
	".\Release\def_exmedi.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_EXM_BE=\
	".\Debug\def_domain.mod"\
	".\Debug\def_exmedi.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\exmedi\exm_ceauxi.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_EXM_C=\
	".\Release\def_domain.mod"\
	".\Release\def_exmedi.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_EXM_C=\
	".\Debug\def_domain.mod"\
	".\Debug\def_exmedi.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\exmedi\exm_ceconc.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_EXM_CE=\
	".\Release\def_domain.mod"\
	".\Release\def_exmedi.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_EXM_CE=\
	".\Debug\def_domain.mod"\
	".\Debug\def_exmedi.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\exmedi\exm_ceicur.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_EXM_CEI=\
	".\Release\def_domain.mod"\
	".\Release\def_exmedi.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_EXM_CEI=\
	".\Debug\def_domain.mod"\
	".\Debug\def_exmedi.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\exmedi\exm_chkpar.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_EXM_CH=\
	".\Release\def_exmedi.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_EXM_CH=\
	".\Debug\def_exmedi.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\exmedi\exm_chkpoi.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_EXM_CHK=\
	".\Release\def_domain.mod"\
	".\Release\def_exmedi.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_iofile.mod"\
	".\Release\Mod_memchk.mod"\
	".\Release\Mod_postpr.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_EXM_CHK=\
	".\Debug\def_domain.mod"\
	".\Debug\def_exmedi.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_iofile.mod"\
	".\Debug\Mod_memchk.mod"\
	".\Debug\Mod_postpr.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\exmedi\exm_comapp.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_EXM_CO=\
	".\Release\def_domain.mod"\
	".\Release\def_exmedi.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_EXM_CO=\
	".\Debug\def_domain.mod"\
	".\Debug\def_exmedi.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\exmedi\exm_comcnd.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_EXM_COM=\
	".\Release\def_domain.mod"\
	".\Release\def_exmedi.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_EXM_COM=\
	".\Debug\def_domain.mod"\
	".\Debug\def_exmedi.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\exmedi\exm_concou.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_EXM_CON=\
	".\Release\def_exmedi.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_EXM_CON=\
	".\Debug\def_exmedi.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\exmedi\exm_cvgunk.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_EXM_CV=\
	".\Release\def_domain.mod"\
	".\Release\def_exmedi.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_EXM_CV=\
	".\Debug\def_domain.mod"\
	".\Debug\def_exmedi.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\exmedi\exm_doiter.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_EXM_D=\
	".\Release\def_domain.mod"\
	".\Release\def_exmedi.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_EXM_D=\
	".\Debug\def_domain.mod"\
	".\Debug\def_exmedi.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\exmedi\exm_eapelm.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_EXM_E=\
	".\Release\def_domain.mod"\
	".\Release\def_exmedi.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_EXM_E=\
	".\Debug\def_domain.mod"\
	".\Debug\def_exmedi.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\exmedi\exm_eapsol.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_EXM_EA=\
	".\Release\def_domain.mod"\
	".\Release\def_exmedi.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_EXM_EA=\
	".\Debug\def_domain.mod"\
	".\Debug\def_exmedi.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\exmedi\exm_eapupd.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_EXM_EAP=\
	".\Release\def_domain.mod"\
	".\Release\def_exmedi.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_EXM_EAP=\
	".\Debug\def_domain.mod"\
	".\Debug\def_exmedi.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\exmedi\exm_endite.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_EXM_EN=\
	".\Release\def_exmedi.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_EXM_EN=\
	".\Debug\def_exmedi.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\exmedi\exm_endste.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_EXM_END=\
	".\Release\def_exmedi.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_EXM_END=\
	".\Debug\def_exmedi.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\exmedi\exm_iapelm.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_EXM_I=\
	".\Release\def_domain.mod"\
	".\Release\def_exmedi.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_EXM_I=\
	".\Debug\def_domain.mod"\
	".\Debug\def_exmedi.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\exmedi\exm_iapode.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_EXM_IA=\
	".\Release\def_domain.mod"\
	".\Release\def_exmedi.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_EXM_IA=\
	".\Debug\def_domain.mod"\
	".\Debug\def_exmedi.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\exmedi\exm_iapsol.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_EXM_IAP=\
	".\Release\def_domain.mod"\
	".\Release\def_exmedi.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_EXM_IAP=\
	".\Debug\def_domain.mod"\
	".\Debug\def_exmedi.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\exmedi\exm_iapupd.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_EXM_IAPU=\
	".\Release\def_domain.mod"\
	".\Release\def_exmedi.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_EXM_IAPU=\
	".\Debug\def_domain.mod"\
	".\Debug\def_exmedi.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\exmedi\exm_iniunk.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_EXM_IN=\
	".\Release\def_domain.mod"\
	".\Release\def_exmedi.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_EXM_IN=\
	".\Debug\def_domain.mod"\
	".\Debug\def_exmedi.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\exmedi\exm_ionicu.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_EXM_IO=\
	".\Release\def_domain.mod"\
	".\Release\def_exmedi.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_EXM_IO=\
	".\Debug\def_domain.mod"\
	".\Debug\def_exmedi.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\exmedi\exm_memall.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_EXM_M=\
	".\Release\def_domain.mod"\
	".\Release\def_exmedi.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_EXM_M=\
	".\Debug\def_domain.mod"\
	".\Debug\def_exmedi.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\exmedi\exm_membcs.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_EXM_ME=\
	".\Release\def_domain.mod"\
	".\Release\def_exmedi.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_EXM_ME=\
	".\Debug\def_domain.mod"\
	".\Debug\def_exmedi.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\exmedi\exm_memose.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_EXM_MEM=\
	".\Release\def_domain.mod"\
	".\Release\def_exmedi.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_iofile.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_EXM_MEM=\
	".\Debug\def_domain.mod"\
	".\Debug\def_exmedi.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_iofile.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\exmedi\exm_memphy.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_EXM_MEMP=\
	".\Release\def_domain.mod"\
	".\Release\def_exmedi.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_EXM_MEMP=\
	".\Debug\def_domain.mod"\
	".\Debug\def_exmedi.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\exmedi\exm_openfi.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_EXM_O=\
	".\Release\def_exmedi.mod"\
	".\Release\def_master.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_EXM_O=\
	".\Debug\def_exmedi.mod"\
	".\Debug\def_master.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\exmedi\exm_outinf.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_EXM_OU=\
	".\Release\def_exmedi.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_EXM_OU=\
	".\Debug\def_exmedi.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\exmedi\exm_output.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_EXM_OUT=\
	".\Release\def_domain.mod"\
	".\Release\def_exmedi.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_iofile.mod"\
	".\Release\mod_output.mod"\
	".\Release\Mod_postpr.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_EXM_OUT=\
	".\Debug\def_domain.mod"\
	".\Debug\def_exmedi.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_iofile.mod"\
	".\Debug\mod_output.mod"\
	".\Debug\Mod_postpr.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\exmedi\exm_outrep.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_EXM_OUTR=\
	".\Release\def_exmedi.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_EXM_OUTR=\
	".\Debug\def_exmedi.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\exmedi\exm_outset.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_EXM_OUTS=\
	".\Release\def_domain.mod"\
	".\Release\def_exmedi.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_EXM_OUTS=\
	".\Debug\def_domain.mod"\
	".\Debug\def_exmedi.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\exmedi\exm_outvar.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_EXM_OUTV=\
	".\Release\def_domain.mod"\
	".\Release\def_exmedi.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_postpr.mod"\
	".\Release\Mod_iofile.mod"\
	".\Release\Mod_memchk.mod"\
	".\Release\Mod_postpr.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_EXM_OUTV=\
	".\Debug\def_domain.mod"\
	".\Debug\def_exmedi.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_postpr.mod"\
	".\Debug\Mod_iofile.mod"\
	".\Debug\Mod_memchk.mod"\
	".\Debug\Mod_postpr.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\exmedi\exm_parall.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_EXM_P=\
	".\Release\def_domain.mod"\
	".\Release\def_exmedi.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_EXM_P=\
	".\Debug\def_domain.mod"\
	".\Debug\def_exmedi.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\exmedi\exm_rcpupd.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_EXM_R=\
	".\Release\def_domain.mod"\
	".\Release\def_exmedi.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_EXM_R=\
	".\Debug\def_domain.mod"\
	".\Debug\def_exmedi.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\exmedi\exm_reabcs.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_EXM_RE=\
	".\Release\def_domain.mod"\
	".\Release\def_exmedi.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_EXM_RE=\
	".\Debug\def_domain.mod"\
	".\Debug\def_exmedi.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\exmedi\exm_reanut.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_EXM_REA=\
	".\Release\def_domain.mod"\
	".\Release\def_exmedi.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_EXM_REA=\
	".\Debug\def_domain.mod"\
	".\Debug\def_exmedi.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\exmedi\exm_reaous.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_EXM_REAO=\
	".\Release\def_domain.mod"\
	".\Release\def_exmedi.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_EXM_REAO=\
	".\Debug\def_domain.mod"\
	".\Debug\def_exmedi.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\exmedi\exm_reaphy.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_EXM_REAP=\
	".\Release\def_domain.mod"\
	".\Release\def_exmedi.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_EXM_REAP=\
	".\Debug\def_domain.mod"\
	".\Debug\def_exmedi.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\exmedi\exm_reapro.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_EXM_REAPR=\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_EXM_REAPR=\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\exmedi\exm_sendat.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_EXM_S=\
	".\Release\def_domain.mod"\
	".\Release\def_exmedi.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_EXM_S=\
	".\Debug\def_domain.mod"\
	".\Debug\def_exmedi.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\exmedi\exm_solite.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_EXM_SO=\
	".\Release\def_domain.mod"\
	".\Release\def_exmedi.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_EXM_SO=\
	".\Debug\def_domain.mod"\
	".\Debug\def_exmedi.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\exmedi\exm_solsys.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_EXM_SOL=\
	".\Release\def_domain.mod"\
	".\Release\def_exmedi.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_EXM_SOL=\
	".\Debug\def_domain.mod"\
	".\Debug\def_exmedi.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\exmedi\exm_timste.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_EXM_T=\
	".\Release\def_exmedi.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_EXM_T=\
	".\Debug\def_exmedi.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\exmedi\exm_turnof.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_EXM_TU=\
	".\Release\def_exmedi.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_EXM_TU=\
	".\Debug\def_exmedi.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\exmedi\exm_turnon.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_EXM_TUR=\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_EXM_TUR=\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\exmedi\exm_updtss.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_EXM_U=\
	".\Release\def_domain.mod"\
	".\Release\def_exmedi.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_EXM_U=\
	".\Debug\def_domain.mod"\
	".\Debug\def_exmedi.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\exmedi\exm_updunk.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_EXM_UP=\
	".\Release\def_domain.mod"\
	".\Release\def_exmedi.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_EXM_UP=\
	".\Debug\def_domain.mod"\
	".\Debug\def_exmedi.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\exmedi\Exmedi.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_EXMED=\
	".\Release\def_domain.mod"\
	".\Release\def_exmedi.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_EXMED=\
	".\Debug\def_domain.mod"\
	".\Debug\def_exmedi.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# End Group
# Begin Group "latbol"

# PROP Default_Filter ""
# End Group
# Begin Group "nastal"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\def_nastal.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_DEF_NA=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_DEF_NA=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\Nastal.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NASTA=\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NASTA=\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_aebody.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_A=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastal.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_A=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastal.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_autbcs.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_AU=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastal.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_AU=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastal.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_begite.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_B=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastal.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_B=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastal.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_begste.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_BE=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastal.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_BE=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastal.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_bounod.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_BO=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastal.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_BO=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastal.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_bouset.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_BOU=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastal.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_BOU=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastal.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_chkpar.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_C=\
	".\Release\def_master.mod"\
	".\Release\def_nastal.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_C=\
	".\Debug\def_master.mod"\
	".\Debug\def_nastal.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_chkpoi.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_CH=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastal.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_iofile.mod"\
	".\Release\Mod_memchk.mod"\
	".\Release\Mod_postpr.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_CH=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastal.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_iofile.mod"\
	".\Debug\Mod_memchk.mod"\
	".\Debug\Mod_postpr.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_conblk.f90
# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_concou.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_CO=\
	".\Release\def_master.mod"\
	".\Release\def_nastal.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_CO=\
	".\Debug\def_master.mod"\
	".\Debug\def_nastal.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_cvgunk.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_CV=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastal.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_CV=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastal.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_doiter.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_D=\
	".\Release\def_master.mod"\
	".\Release\def_nastal.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_D=\
	".\Debug\def_master.mod"\
	".\Debug\def_nastal.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_elchea.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_E=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastal.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_E=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastal.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_elcons.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_EL=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastal.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_EL=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastal.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_elfmom.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_ELF=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastal.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_ELF=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastal.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_elmchl.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_ELM=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_ELM=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_elmset.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_ELMS=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastal.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_ELMS=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastal.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_endite.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_EN=\
	".\Release\def_master.mod"\
	".\Release\def_nastal.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_EN=\
	".\Debug\def_master.mod"\
	".\Debug\def_nastal.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_endste.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_END=\
	".\Release\def_master.mod"\
	".\Release\def_nastal.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_END=\
	".\Debug\def_master.mod"\
	".\Debug\def_nastal.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_exchan.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_EX=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastal.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_EX=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastal.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_funcre.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_F=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_F=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_gacons.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_G=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_G=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_gafmom.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_GA=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_GA=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_gocons.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_GO=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastal.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_GO=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastal.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_gofmom.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_GOF=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastal.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_GOF=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastal.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_inimet.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_I=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastal.mod"\
	".\Release\Mod_postpr.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_I=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastal.mod"\
	".\Debug\Mod_postpr.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_inisol.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_IN=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastal.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_IN=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastal.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_iniunk.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_INI=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastal.mod"\
	".\Release\Mod_postpr.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_INI=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastal.mod"\
	".\Debug\Mod_postpr.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_inivar.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_INIV=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastal.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_INIV=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastal.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_itexco.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_IT=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastal.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_IT=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastal.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_itexin.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_ITE=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastal.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_ITE=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastal.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_jacset.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_J=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastal.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_J=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastal.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_jacvar.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_JA=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastal.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_JA=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastal.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_lawvis.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_L=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastal.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_L=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastal.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_memall.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_M=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastal.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_M=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastal.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_membcs.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_ME=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastal.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_ME=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastal.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_memose.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_MEM=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastal.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_iofile.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_MEM=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastal.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_iofile.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_memphy.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_MEMP=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastal.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_MEMP=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastal.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_modspe.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_MO=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastal.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_MO=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastal.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_openfi.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_O=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastal.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_O=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastal.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_outinf.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_OU=\
	".\Release\def_master.mod"\
	".\Release\def_nastal.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_OU=\
	".\Debug\def_master.mod"\
	".\Debug\def_nastal.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_outpro.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_OUT=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastal.mod"\
	".\Release\Mod_iofile.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_OUT=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastal.mod"\
	".\Debug\Mod_iofile.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_output.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_OUTP=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastal.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_iofile.mod"\
	".\Release\Mod_postpr.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_OUTP=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastal.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_iofile.mod"\
	".\Debug\Mod_postpr.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_outset.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_OUTS=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastal.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_OUTS=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastal.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_outvar.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_OUTV=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastal.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_postpr.mod"\
	".\Release\Mod_iofile.mod"\
	".\Release\Mod_memchk.mod"\
	".\Release\Mod_postpr.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_OUTV=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastal.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_postpr.mod"\
	".\Debug\Mod_iofile.mod"\
	".\Debug\Mod_memchk.mod"\
	".\Debug\Mod_postpr.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_parall.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_P=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastal.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_P=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastal.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_partim.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_PA=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastal.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_PA=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastal.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_reabcp.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_R=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastal.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_R=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastal.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_reabcs.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_RE=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastal.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_RE=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastal.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_reanut.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_REA=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastal.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_REA=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastal.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_reaous.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_REAO=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastal.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_iofile.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_REAO=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastal.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_iofile.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_reaphy.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_REAP=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastal.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_REAP=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastal.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_reapro.f90
# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_relvar.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_REL=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastal.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_REL=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastal.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_rotunk.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_RO=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastal.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_RO=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastal.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_rotvec.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_ROT=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_ROT=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_sendat.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_S=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastal.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_S=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastal.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_setleo.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_SE=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastal.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_SE=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastal.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_setmac.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_SET=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastal.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_SET=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastal.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_setvar.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_SETV=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastal.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_SETV=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastal.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_shocap.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_SH=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_SH=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_solite.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_SO=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastal.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_SO=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastal.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_stalaw.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_ST=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastal.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_ST=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastal.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_timste.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_T=\
	".\Release\def_master.mod"\
	".\Release\def_nastal.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_T=\
	".\Debug\def_master.mod"\
	".\Debug\def_nastal.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_tistep.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_TI=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastal.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_TI=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastal.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_turnof.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_TU=\
	".\Release\def_master.mod"\
	".\Release\def_nastal.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_TU=\
	".\Debug\def_master.mod"\
	".\Debug\def_nastal.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_turnon.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_TUR=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastal.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_TUR=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastal.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_upcons.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_U=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastal.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_U=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastal.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_updbcs.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_UP=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastal.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_UP=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastal.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_updtss.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_UPD=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastal.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_UPD=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastal.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_updunk.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_UPDU=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastal.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_UPDU=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastal.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_upfmom.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_UPF=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastal.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_UPF=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastal.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_uptimi.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_UPT=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_nastal.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_UPT=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_nastal.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nastal\nsa_vortic.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NSA_V=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NSA_V=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# End Group
# Begin Group "solidz"

# PROP Default_Filter "f90"
# Begin Source File

SOURCE=..\..\Sources\modules\solidz\def_solidz.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_DEF_SO=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_DEF_SO=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\solidz\sld_begite.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SLD_B=\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solidz.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SLD_B=\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solidz.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\solidz\sld_begste.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SLD_BE=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solidz.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SLD_BE=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solidz.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\solidz\sld_bouope.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SLD_BO=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_solidz.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SLD_BO=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_solidz.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\solidz\sld_conblk.f90
# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\solidz\sld_concou.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SLD_C=\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solidz.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SLD_C=\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solidz.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\solidz\sld_cvgunk.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SLD_CV=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solidz.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SLD_CV=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solidz.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\solidz\sld_doiter.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SLD_D=\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solidz.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SLD_D=\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solidz.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\solidz\sld_elmcla.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SLD_E=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_solidz.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SLD_E=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_solidz.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\solidz\sld_elmgat.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SLD_EL=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_solidz.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SLD_EL=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_solidz.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\solidz\sld_elmope.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SLD_ELM=\
	".\Release\def_domain.mod"\
	".\Release\def_elmope.mod"\
	".\Release\def_master.mod"\
	".\Release\def_solidz.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SLD_ELM=\
	".\Debug\def_domain.mod"\
	".\Debug\def_elmope.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_solidz.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\solidz\sld_elmpre.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SLD_ELMP=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SLD_ELMP=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\solidz\sld_elmrhs.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SLD_ELMR=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SLD_ELMR=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\solidz\sld_endite.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SLD_EN=\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solidz.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SLD_EN=\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solidz.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\solidz\sld_endste.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SLD_END=\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solidz.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SLD_END=\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solidz.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\solidz\sld_inisol.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SLD_I=\
	".\Release\def_domain.mod"\
	".\Release\def_solidz.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SLD_I=\
	".\Debug\def_domain.mod"\
	".\Debug\def_solidz.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\solidz\sld_iniunk.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SLD_IN=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solidz.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SLD_IN=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solidz.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\solidz\sld_inivar.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SLD_INI=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solidz.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SLD_INI=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solidz.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\solidz\sld_matrix.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SLD_M=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solidz.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SLD_M=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solidz.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\solidz\sld_memall.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SLD_ME=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solidz.mod"\
	".\Release\def_solver.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SLD_ME=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solidz.mod"\
	".\Debug\def_solver.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\solidz\sld_membcs.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SLD_MEM=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solidz.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SLD_MEM=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solidz.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\solidz\sld_memphy.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SLD_MEMP=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_solidz.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SLD_MEMP=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_solidz.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\solidz\sld_openfi.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SLD_O=\
	".\Release\def_master.mod"\
	".\Release\def_solidz.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SLD_O=\
	".\Debug\def_master.mod"\
	".\Debug\def_solidz.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\solidz\sld_outinf.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SLD_OU=\
	".\Release\def_master.mod"\
	".\Release\def_solidz.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SLD_OU=\
	".\Debug\def_master.mod"\
	".\Debug\def_solidz.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\solidz\sld_outlat.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SLD_OUT=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_solidz.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SLD_OUT=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_solidz.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\solidz\sld_output.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SLD_OUTP=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solidz.mod"\
	".\Release\Mod_iofile.mod"\
	".\Release\Mod_postpr.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SLD_OUTP=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solidz.mod"\
	".\Debug\Mod_iofile.mod"\
	".\Debug\Mod_postpr.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\solidz\sld_outvar.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SLD_OUTV=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solidz.mod"\
	".\Release\Mod_postpr.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SLD_OUTV=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solidz.mod"\
	".\Debug\Mod_postpr.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\solidz\sld_reabcp.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SLD_R=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solidz.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SLD_R=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solidz.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\solidz\sld_reabcs.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SLD_RE=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solidz.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SLD_RE=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solidz.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\solidz\sld_reanut.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SLD_REA=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solidz.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SLD_REA=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solidz.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\solidz\sld_reaous.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SLD_REAO=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solidz.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SLD_REAO=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solidz.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\solidz\sld_reaphy.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SLD_REAP=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solidz.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SLD_REAP=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solidz.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\solidz\sld_reapro.f90
# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\solidz\sld_restar.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SLD_RES=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solidz.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SLD_RES=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solidz.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\solidz\sld_solexp.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SLD_S=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solidz.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SLD_S=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solidz.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\solidz\sld_solite.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SLD_SO=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solidz.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SLD_SO=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solidz.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\solidz\sld_timste.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SLD_T=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solidz.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SLD_T=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solidz.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\solidz\sld_tistep.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SLD_TI=\
	".\Release\def_master.mod"\
	".\Release\def_solidz.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SLD_TI=\
	".\Debug\def_master.mod"\
	".\Debug\def_solidz.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\solidz\sld_turnof.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SLD_TU=\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solidz.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SLD_TU=\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solidz.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\solidz\sld_turnon.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SLD_TUR=\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solidz.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SLD_TUR=\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solidz.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\solidz\sld_updtss.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SLD_U=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solidz.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SLD_U=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solidz.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\solidz\sld_updunk.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SLD_UP=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solidz.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SLD_UP=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solidz.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\solidz\Solidz.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SOLID=\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SOLID=\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# End Group
# Begin Group "gotita"

# PROP Default_Filter "f90"
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\def_gotita.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_DEF_G=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_DEF_G=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_begite.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_B=\
	".\Release\def_domain.mod"\
	".\Release\def_gotita.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_B=\
	".\Debug\def_domain.mod"\
	".\Debug\def_gotita.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_begste.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_BE=\
	".\Release\def_domain.mod"\
	".\Release\def_gotita.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_BE=\
	".\Debug\def_domain.mod"\
	".\Debug\def_gotita.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_conblk.f90
# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_concou.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_C=\
	".\Release\def_gotita.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_C=\
	".\Debug\def_gotita.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_cvgunk.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_CV=\
	".\Release\def_domain.mod"\
	".\Release\def_gotita.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_CV=\
	".\Debug\def_domain.mod"\
	".\Debug\def_gotita.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_diffun.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_D=\
	".\Release\def_gotita.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_D=\
	".\Debug\def_gotita.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_diffus.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_DI=\
	".\Release\def_domain.mod"\
	".\Release\def_gotita.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_DI=\
	".\Debug\def_domain.mod"\
	".\Debug\def_gotita.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_doiter.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_DO=\
	".\Release\def_gotita.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_DO=\
	".\Debug\def_gotita.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_dragco.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_DR=\
	".\Release\def_gotita.mod"\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_DR=\
	".\Debug\def_gotita.mod"\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_elmdir.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_E=\
	".\Release\def_gotita.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_E=\
	".\Debug\def_gotita.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_elmexa.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_EL=\
	".\Release\def_domain.mod"\
	".\Release\def_gotita.mod"\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_EL=\
	".\Debug\def_domain.mod"\
	".\Debug\def_gotita.mod"\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_elmgat.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_ELM=\
	".\Release\def_domain.mod"\
	".\Release\def_gotita.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_ELM=\
	".\Debug\def_domain.mod"\
	".\Debug\def_gotita.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_elmma2.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_ELMM=\
	".\Release\def_domain.mod"\
	".\Release\def_gotita.mod"\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_ELMM=\
	".\Debug\def_domain.mod"\
	".\Debug\def_gotita.mod"\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_elmma3.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_ELMMA=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_ELMMA=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_elmmat.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_ELMMAT=\
	".\Release\def_domain.mod"\
	".\Release\def_gotita.mod"\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_ELMMAT=\
	".\Debug\def_domain.mod"\
	".\Debug\def_gotita.mod"\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_elmop2.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_ELMO=\
	".\Release\def_domain.mod"\
	".\Release\def_gotita.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_ELMO=\
	".\Debug\def_domain.mod"\
	".\Debug\def_gotita.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_elmop3.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_ELMOP=\
	".\Release\def_domain.mod"\
	".\Release\def_gotita.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_ELMOP=\
	".\Debug\def_domain.mod"\
	".\Debug\def_gotita.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_elmope.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_ELMOPE=\
	".\Release\def_domain.mod"\
	".\Release\def_gotita.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_ELMOPE=\
	".\Debug\def_domain.mod"\
	".\Debug\def_gotita.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_elmpre.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_ELMP=\
	".\Release\def_domain.mod"\
	".\Release\def_gotita.mod"\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_ELMP=\
	".\Debug\def_domain.mod"\
	".\Debug\def_gotita.mod"\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_elmpro.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_ELMPR=\
	".\Release\def_domain.mod"\
	".\Release\def_gotita.mod"\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_ELMPR=\
	".\Debug\def_domain.mod"\
	".\Debug\def_gotita.mod"\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_elmre2.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_ELMR=\
	".\Release\def_domain.mod"\
	".\Release\def_gotita.mod"\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_ELMR=\
	".\Debug\def_domain.mod"\
	".\Debug\def_gotita.mod"\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_elmre3.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_ELMRE=\
	".\Release\def_domain.mod"\
	".\Release\def_gotita.mod"\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_ELMRE=\
	".\Debug\def_domain.mod"\
	".\Debug\def_gotita.mod"\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_elmres.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_ELMRES=\
	".\Release\def_domain.mod"\
	".\Release\def_gotita.mod"\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_ELMRES=\
	".\Debug\def_domain.mod"\
	".\Debug\def_gotita.mod"\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_elmrhs.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_ELMRH=\
	".\Release\def_domain.mod"\
	".\Release\def_gotita.mod"\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_ELMRH=\
	".\Debug\def_domain.mod"\
	".\Debug\def_gotita.mod"\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_elmsc3.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_ELMS=\
	".\Release\def_domain.mod"\
	".\Release\def_gotita.mod"\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_ELMS=\
	".\Debug\def_domain.mod"\
	".\Debug\def_gotita.mod"\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_elmsgs.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_ELMSG=\
	".\Release\def_domain.mod"\
	".\Release\def_gotita.mod"\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_ELMSG=\
	".\Debug\def_domain.mod"\
	".\Debug\def_gotita.mod"\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_elmshc.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_ELMSH=\
	".\Release\def_domain.mod"\
	".\Release\def_gotita.mod"\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_ELMSH=\
	".\Debug\def_domain.mod"\
	".\Debug\def_gotita.mod"\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_elmshm.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_ELMSHM=\
	".\Release\def_domain.mod"\
	".\Release\def_gotita.mod"\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_ELMSHM=\
	".\Debug\def_domain.mod"\
	".\Debug\def_gotita.mod"\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_elmsm2.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_ELMSM=\
	".\Release\def_domain.mod"\
	".\Release\def_gotita.mod"\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_ELMSM=\
	".\Debug\def_domain.mod"\
	".\Debug\def_gotita.mod"\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_elmte2.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_ELMT=\
	".\Release\def_domain.mod"\
	".\Release\def_gotita.mod"\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_ELMT=\
	".\Debug\def_domain.mod"\
	".\Debug\def_gotita.mod"\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_elmte3.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_ELMTE=\
	".\Release\def_domain.mod"\
	".\Release\def_gotita.mod"\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_ELMTE=\
	".\Debug\def_domain.mod"\
	".\Debug\def_gotita.mod"\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_elmtes.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_ELMTES=\
	".\Release\def_domain.mod"\
	".\Release\def_gotita.mod"\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_ELMTES=\
	".\Debug\def_domain.mod"\
	".\Debug\def_gotita.mod"\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_elmtss.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_ELMTS=\
	".\Release\def_domain.mod"\
	".\Release\def_gotita.mod"\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_ELMTS=\
	".\Debug\def_domain.mod"\
	".\Debug\def_gotita.mod"\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_endite.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_EN=\
	".\Release\def_gotita.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_EN=\
	".\Debug\def_gotita.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_endste.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_END=\
	".\Release\def_gotita.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_END=\
	".\Debug\def_gotita.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_exacso.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_EX=\
	".\Release\def_domain.mod"\
	".\Release\def_gotita.mod"\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_EX=\
	".\Debug\def_domain.mod"\
	".\Debug\def_gotita.mod"\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_exaerr.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_EXA=\
	".\Release\def_domain.mod"\
	".\Release\def_gotita.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	".\Release\Mod_postpr.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_EXA=\
	".\Debug\def_domain.mod"\
	".\Debug\def_gotita.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	".\Debug\Mod_postpr.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_inisol.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_I=\
	".\Release\def_domain.mod"\
	".\Release\def_gotita.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_I=\
	".\Debug\def_domain.mod"\
	".\Debug\def_gotita.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_iniunk.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_IN=\
	".\Release\def_domain.mod"\
	".\Release\def_gotita.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_IN=\
	".\Debug\def_domain.mod"\
	".\Debug\def_gotita.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_inivar.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_INI=\
	".\Release\def_domain.mod"\
	".\Release\def_gotita.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_INI=\
	".\Debug\def_domain.mod"\
	".\Debug\def_gotita.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_intrst.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_INT=\
	".\Release\def_domain.mod"\
	".\Release\def_gotita.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_iofile.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_INT=\
	".\Debug\def_domain.mod"\
	".\Debug\def_gotita.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_iofile.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_matrix.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_M=\
	".\Release\def_gotita.mod"\
	".\Release\def_master.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_M=\
	".\Debug\def_gotita.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_memall.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_ME=\
	".\Release\def_domain.mod"\
	".\Release\def_gotita.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_ME=\
	".\Debug\def_domain.mod"\
	".\Debug\def_gotita.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_membcs.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_MEM=\
	".\Release\def_domain.mod"\
	".\Release\def_gotita.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_MEM=\
	".\Debug\def_domain.mod"\
	".\Debug\def_gotita.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_memphy.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_MEMP=\
	".\Release\def_domain.mod"\
	".\Release\def_gotita.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_MEMP=\
	".\Debug\def_domain.mod"\
	".\Debug\def_gotita.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_newmsh.f90
# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_openfi.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_O=\
	".\Release\def_gotita.mod"\
	".\Release\def_master.mod"\
	".\Release\def_postpr.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_O=\
	".\Debug\def_gotita.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_postpr.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_outerr.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_OU=\
	".\Release\def_domain.mod"\
	".\Release\def_gotita.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_OU=\
	".\Debug\def_domain.mod"\
	".\Debug\def_gotita.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_outinf.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_OUT=\
	".\Release\def_gotita.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_OUT=\
	".\Debug\def_gotita.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_outlat.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_OUTL=\
	".\Release\def_domain.mod"\
	".\Release\def_gotita.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_OUTL=\
	".\Debug\def_domain.mod"\
	".\Debug\def_gotita.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_output.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_OUTP=\
	".\Release\def_domain.mod"\
	".\Release\def_gotita.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_iofile.mod"\
	".\Release\Mod_postpr.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_OUTP=\
	".\Debug\def_domain.mod"\
	".\Debug\def_gotita.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_iofile.mod"\
	".\Debug\Mod_postpr.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_outvar.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_OUTV=\
	".\Release\def_domain.mod"\
	".\Release\def_gotita.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\mod_gradie.mod"\
	".\Release\Mod_postpr.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_OUTV=\
	".\Debug\def_domain.mod"\
	".\Debug\def_gotita.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\mod_gradie.mod"\
	".\Debug\Mod_postpr.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_parall.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_P=\
	".\Release\def_gotita.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_P=\
	".\Debug\def_gotita.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_reabcs.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_R=\
	".\Release\def_domain.mod"\
	".\Release\def_gotita.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_R=\
	".\Debug\def_domain.mod"\
	".\Debug\def_gotita.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_reanut.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_RE=\
	".\Release\def_domain.mod"\
	".\Release\def_gotita.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_RE=\
	".\Debug\def_domain.mod"\
	".\Debug\def_gotita.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_reaous.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_REA=\
	".\Release\def_gotita.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_REA=\
	".\Debug\def_gotita.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_reaphy.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_REAP=\
	".\Release\def_domain.mod"\
	".\Release\def_gotita.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_REAP=\
	".\Debug\def_domain.mod"\
	".\Debug\def_gotita.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_reapro.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_REAPR=\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_REAPR=\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_restar.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_RES=\
	".\Release\def_domain.mod"\
	".\Release\def_gotita.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_RES=\
	".\Debug\def_domain.mod"\
	".\Debug\def_gotita.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_sendat.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_S=\
	".\Release\def_domain.mod"\
	".\Release\def_gotita.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_S=\
	".\Debug\def_domain.mod"\
	".\Debug\def_gotita.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_solbgs.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_SO=\
	".\Release\def_domain.mod"\
	".\Release\def_gotita.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_SO=\
	".\Debug\def_domain.mod"\
	".\Debug\def_gotita.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_solite.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_SOL=\
	".\Release\def_domain.mod"\
	".\Release\def_gotita.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_SOL=\
	".\Debug\def_domain.mod"\
	".\Debug\def_gotita.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_solmon.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_SOLM=\
	".\Release\def_gotita.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_SOLM=\
	".\Debug\def_gotita.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_timste.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_T=\
	".\Release\def_domain.mod"\
	".\Release\def_gotita.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_T=\
	".\Debug\def_domain.mod"\
	".\Debug\def_gotita.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_tistep.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_TI=\
	".\Release\def_gotita.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_TI=\
	".\Debug\def_gotita.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_turnof.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_TU=\
	".\Release\def_gotita.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_TU=\
	".\Debug\def_gotita.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_turnon.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_TUR=\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_TUR=\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_updfix.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_U=\
	".\Release\def_domain.mod"\
	".\Release\def_gotita.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_U=\
	".\Debug\def_domain.mod"\
	".\Debug\def_gotita.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_updtss.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_UP=\
	".\Release\def_domain.mod"\
	".\Release\def_gotita.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_UP=\
	".\Debug\def_domain.mod"\
	".\Debug\def_gotita.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_updunk.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_UPD=\
	".\Release\def_domain.mod"\
	".\Release\def_gotita.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_UPD=\
	".\Debug\def_domain.mod"\
	".\Debug\def_gotita.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\got_veloci.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOT_V=\
	".\Release\def_domain.mod"\
	".\Release\def_gotita.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOT_V=\
	".\Debug\def_domain.mod"\
	".\Debug\def_gotita.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\gotita\Gotita.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_GOTIT=\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_GOTIT=\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# End Group
# Begin Group "wavequ"

# PROP Default_Filter "f90"
# Begin Source File

SOURCE=..\..\Sources\modules\wavequ\def_wavequ.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_DEF_W=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_DEF_W=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\wavequ\wav_addmas.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_WAV_A=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_WAV_A=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\wavequ\wav_assmat.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_WAV_AS=\
	".\Release\def_kintyp.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_WAV_AS=\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\wavequ\wav_begite.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_WAV_B=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_wavequ.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_WAV_B=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_wavequ.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\wavequ\wav_begste.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_WAV_BE=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_wavequ.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_WAV_BE=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_wavequ.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\wavequ\wav_bouope.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_WAV_BO=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_wavequ.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_WAV_BO=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_wavequ.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\wavequ\wav_conblk.f90
# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\wavequ\wav_concou.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_WAV_C=\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_wavequ.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_WAV_C=\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_wavequ.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\wavequ\wav_cvgunk.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_WAV_CV=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_wavequ.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_WAV_CV=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_wavequ.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\wavequ\wav_doiter.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_WAV_D=\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	".\Release\def_wavequ.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_WAV_D=\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	".\Debug\def_wavequ.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\wavequ\wav_elmope.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_WAV_E=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_wavequ.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_WAV_E=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_wavequ.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\wavequ\wav_elmsou.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_WAV_EL=\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_WAV_EL=\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\wavequ\wav_endite.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_WAV_EN=\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_wavequ.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_WAV_EN=\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_wavequ.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\wavequ\wav_endste.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_WAV_END=\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_wavequ.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_WAV_END=\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_wavequ.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\wavequ\wav_explic.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_WAV_EX=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_wavequ.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_WAV_EX=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_wavequ.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\wavequ\wav_inisol.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_WAV_I=\
	".\Release\def_solver.mod"\
	".\Release\def_wavequ.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_WAV_I=\
	".\Debug\def_solver.mod"\
	".\Debug\def_wavequ.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\wavequ\wav_iniunk.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_WAV_IN=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_wavequ.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_WAV_IN=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_wavequ.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\wavequ\wav_inivar.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_WAV_INI=\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	".\Release\def_wavequ.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_WAV_INI=\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	".\Debug\def_wavequ.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\wavequ\wav_matrix.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_WAV_M=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	".\Release\def_wavequ.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_WAV_M=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	".\Debug\def_wavequ.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\wavequ\wav_memall.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_WAV_ME=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	".\Release\def_wavequ.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_WAV_ME=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	".\Debug\def_wavequ.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\wavequ\wav_membcs.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_WAV_MEM=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_wavequ.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_WAV_MEM=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_wavequ.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\wavequ\wav_memphy.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_WAV_MEMP=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_wavequ.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_WAV_MEMP=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_wavequ.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\wavequ\wav_newmsh.f90
# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\wavequ\wav_openfi.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_WAV_O=\
	".\Release\def_master.mod"\
	".\Release\def_wavequ.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_WAV_O=\
	".\Debug\def_master.mod"\
	".\Debug\def_wavequ.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\wavequ\wav_outerr.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_WAV_OU=\
	".\Release\def_master.mod"\
	".\Release\def_wavequ.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_WAV_OU=\
	".\Debug\def_master.mod"\
	".\Debug\def_wavequ.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\wavequ\wav_outinf.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_WAV_OUT=\
	".\Release\def_master.mod"\
	".\Release\def_wavequ.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_WAV_OUT=\
	".\Debug\def_master.mod"\
	".\Debug\def_wavequ.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\wavequ\wav_output.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_WAV_OUTP=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_wavequ.mod"\
	".\Release\Mod_iofile.mod"\
	".\Release\Mod_postpr.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_WAV_OUTP=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_wavequ.mod"\
	".\Debug\Mod_iofile.mod"\
	".\Debug\Mod_postpr.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\wavequ\wav_outset.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_WAV_OUTS=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_wavequ.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_WAV_OUTS=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_wavequ.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\wavequ\wav_outvar.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_WAV_OUTV=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_wavequ.mod"\
	".\Release\Mod_memchk.mod"\
	".\Release\Mod_postpr.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_WAV_OUTV=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_wavequ.mod"\
	".\Debug\Mod_memchk.mod"\
	".\Debug\Mod_postpr.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\wavequ\wav_parall.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_WAV_P=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_wavequ.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_WAV_P=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_wavequ.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\wavequ\wav_reabcs.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_WAV_R=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_wavequ.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_WAV_R=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_wavequ.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\wavequ\wav_reanut.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_WAV_RE=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	".\Release\def_wavequ.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_WAV_RE=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	".\Debug\def_wavequ.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\wavequ\wav_reaous.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_WAV_REA=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_wavequ.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_WAV_REA=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_wavequ.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\wavequ\wav_reaphy.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_WAV_REAP=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_wavequ.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_WAV_REAP=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_wavequ.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\wavequ\wav_reapro.f90
# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\wavequ\wav_restar.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_WAV_RES=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_wavequ.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_WAV_RES=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_wavequ.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\wavequ\wav_sendat.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_WAV_S=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	".\Release\def_wavequ.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_WAV_S=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	".\Debug\def_wavequ.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\wavequ\wav_solite.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_WAV_SO=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_wavequ.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_WAV_SO=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_wavequ.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\wavequ\wav_timste.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_WAV_T=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_wavequ.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_WAV_T=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_wavequ.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\wavequ\wav_tistep.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_WAV_TI=\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_wavequ.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_WAV_TI=\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_wavequ.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\wavequ\wav_turnof.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_WAV_TU=\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	".\Release\def_wavequ.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_WAV_TU=\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	".\Debug\def_wavequ.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\wavequ\wav_turnon.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_WAV_TUR=\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_wavequ.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_WAV_TUR=\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_wavequ.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\wavequ\wav_updtss.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_WAV_U=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_wavequ.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_WAV_U=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_wavequ.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\wavequ\wav_updunk.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_WAV_UP=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_wavequ.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_WAV_UP=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_wavequ.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\wavequ\Wavequ.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_WAVEQ=\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_WAVEQ=\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# End Group
# Begin Group "levels"

# PROP Default_Filter "f90"
# Begin Source File

SOURCE=..\..\Sources\modules\levels\def_levels.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_DEF_L=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_DEF_L=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\levels\lev_assmat.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_LEV_A=\
	".\Release\def_kintyp.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_LEV_A=\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\levels\lev_begite.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_LEV_B=\
	".\Release\def_domain.mod"\
	".\Release\def_levels.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_LEV_B=\
	".\Debug\def_domain.mod"\
	".\Debug\def_levels.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\levels\lev_begste.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_LEV_BE=\
	".\Release\def_domain.mod"\
	".\Release\def_levels.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_LEV_BE=\
	".\Debug\def_domain.mod"\
	".\Debug\def_levels.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\levels\lev_calvol.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_LEV_C=\
	".\Release\def_domain.mod"\
	".\Release\def_levels.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_LEV_C=\
	".\Debug\def_domain.mod"\
	".\Debug\def_levels.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\levels\lev_conblk.f90
# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\levels\lev_concou.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_LEV_CO=\
	".\Release\def_levels.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_LEV_CO=\
	".\Debug\def_levels.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\levels\lev_conint.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_LEV_CON=\
	".\Release\def_domain.mod"\
	".\Release\def_levels.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_LEV_CON=\
	".\Debug\def_domain.mod"\
	".\Debug\def_levels.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\levels\lev_corvol.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_LEV_COR=\
	".\Release\def_domain.mod"\
	".\Release\def_levels.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_LEV_COR=\
	".\Debug\def_domain.mod"\
	".\Debug\def_levels.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\levels\lev_cvgunk.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_LEV_CV=\
	".\Release\def_domain.mod"\
	".\Release\def_levels.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_LEV_CV=\
	".\Debug\def_domain.mod"\
	".\Debug\def_levels.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\levels\lev_doiter.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_LEV_D=\
	".\Release\def_levels.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_LEV_D=\
	".\Debug\def_levels.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\levels\lev_elmdir.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_LEV_E=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_levels.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_LEV_E=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_levels.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\levels\lev_elmini.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_LEV_EL=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_LEV_EL=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\levels\lev_elmmat.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_LEV_ELM=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_LEV_ELM=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\levels\lev_elmope.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_LEV_ELMO=\
	".\Release\def_domain.mod"\
	".\Release\def_levels.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_LEV_ELMO=\
	".\Debug\def_domain.mod"\
	".\Debug\def_levels.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\levels\lev_endite.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_LEV_EN=\
	".\Release\def_levels.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_LEV_EN=\
	".\Debug\def_levels.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\levels\lev_endste.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_LEV_END=\
	".\Release\def_domain.mod"\
	".\Release\def_levels.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_LEV_END=\
	".\Debug\def_domain.mod"\
	".\Debug\def_levels.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\levels\lev_explic.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_LEV_EX=\
	".\Release\def_domain.mod"\
	".\Release\def_levels.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_LEV_EX=\
	".\Debug\def_domain.mod"\
	".\Debug\def_levels.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\levels\lev_inisol.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_LEV_I=\
	".\Release\def_domain.mod"\
	".\Release\def_levels.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_LEV_I=\
	".\Debug\def_domain.mod"\
	".\Debug\def_levels.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\levels\lev_iniunk.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_LEV_IN=\
	".\Release\def_domain.mod"\
	".\Release\def_levels.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_LEV_IN=\
	".\Debug\def_domain.mod"\
	".\Debug\def_levels.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\levels\lev_inivar.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_LEV_INI=\
	".\Release\def_levels.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_LEV_INI=\
	".\Debug\def_levels.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\levels\lev_matrix.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_LEV_M=\
	".\Release\def_domain.mod"\
	".\Release\def_levels.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_LEV_M=\
	".\Debug\def_domain.mod"\
	".\Debug\def_levels.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\levels\lev_memall.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_LEV_ME=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_levels.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_LEV_ME=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_levels.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\levels\lev_membcs.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_LEV_MEM=\
	".\Release\def_domain.mod"\
	".\Release\def_levels.mod"\
	".\Release\def_master.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_LEV_MEM=\
	".\Debug\def_domain.mod"\
	".\Debug\def_levels.mod"\
	".\Debug\def_master.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\levels\lev_newmsh.f90
# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\levels\lev_openfi.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_LEV_O=\
	".\Release\def_levels.mod"\
	".\Release\def_master.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_LEV_O=\
	".\Debug\def_levels.mod"\
	".\Debug\def_master.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\levels\lev_outerr.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_LEV_OU=\
	".\Release\def_levels.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_LEV_OU=\
	".\Debug\def_levels.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\levels\lev_outinf.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_LEV_OUT=\
	".\Release\def_levels.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_LEV_OUT=\
	".\Debug\def_levels.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\levels\lev_output.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_LEV_OUTP=\
	".\Release\def_domain.mod"\
	".\Release\def_levels.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_iofile.mod"\
	".\Release\Mod_postpr.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_LEV_OUTP=\
	".\Debug\def_domain.mod"\
	".\Debug\def_levels.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_iofile.mod"\
	".\Debug\Mod_postpr.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\levels\lev_outset.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_LEV_OUTS=\
	".\Release\def_domain.mod"\
	".\Release\def_levels.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_LEV_OUTS=\
	".\Debug\def_domain.mod"\
	".\Debug\def_levels.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\levels\lev_outvar.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_LEV_OUTV=\
	".\Release\def_domain.mod"\
	".\Release\def_levels.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_LEV_OUTV=\
	".\Debug\def_domain.mod"\
	".\Debug\def_levels.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\levels\lev_parall.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_LEV_P=\
	".\Release\def_domain.mod"\
	".\Release\def_levels.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_LEV_P=\
	".\Debug\def_domain.mod"\
	".\Debug\def_levels.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\levels\lev_reabcs.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_LEV_R=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_levels.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_LEV_R=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_levels.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\levels\lev_reanut.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_LEV_RE=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_levels.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_LEV_RE=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_levels.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\levels\lev_reaous.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_LEV_REA=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_levels.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_LEV_REA=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_levels.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\levels\lev_reaphy.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_LEV_REAP=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_levels.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_LEV_REAP=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_levels.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\levels\lev_reapro.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_LEV_REAPR=\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_LEV_REAPR=\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\levels\lev_redist.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_LEV_RED=\
	".\Release\def_domain.mod"\
	".\Release\def_levels.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_LEV_RED=\
	".\Debug\def_domain.mod"\
	".\Debug\def_levels.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\levels\lev_restar.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_LEV_RES=\
	".\Release\def_domain.mod"\
	".\Release\def_levels.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_LEV_RES=\
	".\Debug\def_domain.mod"\
	".\Debug\def_levels.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\levels\lev_sendat.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_LEV_S=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_levels.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_LEV_S=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_levels.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\levels\lev_solite.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_LEV_SO=\
	".\Release\def_domain.mod"\
	".\Release\def_levels.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_LEV_SO=\
	".\Debug\def_domain.mod"\
	".\Debug\def_levels.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\levels\lev_timste.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_LEV_T=\
	".\Release\def_domain.mod"\
	".\Release\def_levels.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_LEV_T=\
	".\Debug\def_domain.mod"\
	".\Debug\def_levels.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\levels\lev_tistep.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_LEV_TI=\
	".\Release\def_levels.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_LEV_TI=\
	".\Debug\def_levels.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\levels\lev_turnof.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_LEV_TU=\
	".\Release\def_levels.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_LEV_TU=\
	".\Debug\def_levels.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\levels\lev_turnon.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_LEV_TUR=\
	".\Release\def_levels.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_LEV_TUR=\
	".\Debug\def_levels.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\levels\lev_updtss.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_LEV_U=\
	".\Release\def_domain.mod"\
	".\Release\def_levels.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_LEV_U=\
	".\Debug\def_domain.mod"\
	".\Debug\def_levels.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\levels\lev_updunk.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_LEV_UP=\
	".\Release\def_domain.mod"\
	".\Release\def_levels.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_LEV_UP=\
	".\Debug\def_domain.mod"\
	".\Debug\def_levels.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\levels\lev_velfun.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_LEV_V=\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_LEV_V=\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\levels\Levels.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_LEVEL=\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_LEVEL=\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# End Group
# Begin Group "quanty"

# PROP Default_Filter "f90"
# Begin Source File

SOURCE=..\..\Sources\modules\quanty\def_quanty.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_DEF_Q=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_DEF_Q=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\quanty\qua_inivar.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_QUA_I=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_quanty.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_QUA_I=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_quanty.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\quanty\qua_openfi.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_QUA_O=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_quanty.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_QUA_O=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_quanty.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\quanty\qua_reaphy.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_QUA_R=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_quanty.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_QUA_R=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_quanty.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\quanty\qua_reapro.f90
# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\quanty\qua_turnon.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_QUA_T=\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_quanty.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_QUA_T=\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_quanty.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\quanty\Quanty.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_QUANT=\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_QUANT=\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# End Group
# Begin Group "partis"

# PROP Default_Filter "f90"
# Begin Source File

SOURCE=..\..\Sources\modules\partis\def_partis.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_DEF_PA=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_DEF_PA=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\partis\Partis.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PARTI=\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PARTI=\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\partis\pts_accumu.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PTS_A=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_partis.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PTS_A=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_partis.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\partis\pts_begite.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PTS_B=\
	".\Release\def_master.mod"\
	".\Release\def_partis.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PTS_B=\
	".\Debug\def_master.mod"\
	".\Debug\def_partis.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\partis\pts_begste.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PTS_BE=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_partis.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PTS_BE=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_partis.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\partis\pts_bouset.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PTS_BO=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_partis.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PTS_BO=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_partis.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\partis\pts_conblk.f90
# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\partis\pts_concou.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PTS_C=\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_partis.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PTS_C=\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_partis.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\partis\pts_cvgunk.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PTS_CV=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_partis.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PTS_CV=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_partis.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\partis\pts_doiter.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PTS_D=\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_partis.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PTS_D=\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_partis.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\partis\pts_elmdir.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PTS_E=\
	".\Release\def_kintyp.mod"\
	".\Release\def_partis.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PTS_E=\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_partis.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\partis\pts_elmgat.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PTS_EL=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_partis.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PTS_EL=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_partis.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\partis\pts_elmope.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PTS_ELM=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_partis.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PTS_ELM=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_partis.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\partis\pts_elmpro.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PTS_ELMP=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_partis.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PTS_ELMP=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_partis.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\partis\pts_endite.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PTS_EN=\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_partis.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PTS_EN=\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_partis.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\partis\pts_endste.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PTS_END=\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_partis.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PTS_END=\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_partis.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\partis\pts_hdiffu.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PTS_H=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_partis.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PTS_H=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_partis.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\partis\pts_iniunk.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PTS_I=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_partis.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PTS_I=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_partis.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\partis\pts_inivar.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PTS_IN=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_partis.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PTS_IN=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_partis.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\partis\pts_matrix.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PTS_M=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_partis.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PTS_M=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_partis.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\partis\pts_memall.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PTS_ME=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_partis.mod"\
	".\Release\def_solver.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PTS_ME=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_partis.mod"\
	".\Debug\def_solver.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\partis\pts_membcs.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PTS_MEM=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_partis.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PTS_MEM=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_partis.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\partis\pts_memose.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PTS_MEMO=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_partis.mod"\
	".\Release\Mod_iofile.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PTS_MEMO=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_partis.mod"\
	".\Debug\Mod_iofile.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\partis\pts_memphy.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PTS_MEMP=\
	".\Release\def_master.mod"\
	".\Release\def_partis.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PTS_MEMP=\
	".\Debug\def_master.mod"\
	".\Debug\def_partis.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\partis\pts_newmsh.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PTS_N=\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PTS_N=\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\partis\pts_openfi.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PTS_O=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_partis.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PTS_O=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_partis.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\partis\pts_outerr.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PTS_OU=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_partis.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PTS_OU=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_partis.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\partis\pts_outinf.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PTS_OUT=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_partis.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PTS_OUT=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_partis.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\partis\pts_output.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PTS_OUTP=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_partis.mod"\
	".\Release\Mod_iofile.mod"\
	".\Release\Mod_postpr.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PTS_OUTP=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_partis.mod"\
	".\Debug\Mod_iofile.mod"\
	".\Debug\Mod_postpr.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\partis\pts_outset.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PTS_OUTS=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_partis.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PTS_OUTS=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_partis.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\partis\pts_outvar.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PTS_OUTV=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_partis.mod"\
	".\Release\Mod_postpr.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PTS_OUTV=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_partis.mod"\
	".\Debug\Mod_postpr.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\partis\pts_outwit.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PTS_OUTW=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_partis.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PTS_OUTW=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_partis.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\partis\pts_parall.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PTS_P=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_partis.mod"\
	".\Release\def_solver.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PTS_P=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_partis.mod"\
	".\Debug\def_solver.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\partis\pts_reabcs.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PTS_R=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_partis.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PTS_R=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_partis.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\partis\pts_reamet.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PTS_RE=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_partis.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PTS_RE=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_partis.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\partis\pts_reanut.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PTS_REA=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_partis.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PTS_REA=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_partis.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\partis\pts_reaous.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PTS_REAO=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_partis.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PTS_REAO=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_partis.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\partis\pts_reaphy.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PTS_REAP=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_partis.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PTS_REAP=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_partis.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\partis\pts_reapro.f90
# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\partis\pts_reasou.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PTS_REAS=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_partis.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PTS_REAS=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_partis.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\partis\pts_restar.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PTS_RES=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_partis.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PTS_RES=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_partis.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\partis\pts_setpsi.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PTS_S=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PTS_S=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\partis\pts_solite.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PTS_SO=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_partis.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PTS_SO=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_partis.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\partis\pts_timste.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PTS_T=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_partis.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PTS_T=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_partis.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\partis\pts_tistep.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PTS_TI=\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_partis.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PTS_TI=\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_partis.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\partis\pts_turnof.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PTS_TU=\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_partis.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PTS_TU=\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_partis.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\partis\pts_turnon.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PTS_TUR=\
	".\Release\def_master.mod"\
	".\Release\def_partis.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PTS_TUR=\
	".\Debug\def_master.mod"\
	".\Debug\def_partis.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\partis\pts_updbcs.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PTS_U=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_partis.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PTS_U=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_partis.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\partis\pts_updunk.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PTS_UP=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_partis.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PTS_UP=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_partis.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\partis\pts_updvte.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PTS_UPD=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_partis.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PTS_UPD=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_partis.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\partis\pts_velfun.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PTS_V=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_partis.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PTS_V=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_partis.mod"\
	

!ENDIF 

# End Source File
# End Group
# Begin Group "nasedg"

# PROP Default_Filter "f90"
# Begin Source File

SOURCE=..\..\Sources\modules\nasedg\chksub.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CHKSU=\
	".\Release\def_kintyp.mod"\
	".\Release\def_nasedg.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CHKSU=\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_nasedg.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nasedg\comp.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_COMP_=\
	".\Release\def_kintyp.mod"\
	".\Release\def_meshin.mod"\
	".\Release\def_nasedg.mod"\
	".\Release\Mod_memchk.mod"\
	".\Release\mod_mshtol.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_COMP_=\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_meshin.mod"\
	".\Debug\def_nasedg.mod"\
	".\Debug\Mod_memchk.mod"\
	".\Debug\mod_mshtol.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nasedg\def_nasedg.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_DEF_NAS=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_DEF_NAS=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nasedg\edgtol.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_EDGTO=\
	".\Release\def_kintyp.mod"\
	".\Release\def_meshin.mod"\
	".\Release\def_nasedg.mod"\
	".\Release\Mod_memchk.mod"\
	".\Release\mod_mshtol.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_EDGTO=\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_meshin.mod"\
	".\Debug\def_nasedg.mod"\
	".\Debug\Mod_memchk.mod"\
	".\Debug\mod_mshtol.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nasedg\Nasedg.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_NASED=\
	".\Release\def_master.mod"\
	".\Release\def_nasedg.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_NASED=\
	".\Debug\def_master.mod"\
	".\Debug\def_nasedg.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nasedg\precond.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PRECO=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PRECO=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nasedg\renum.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_RENUM_=\
	".\Release\def_kintyp.mod"\
	".\Release\def_nasedg.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_RENUM_=\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_nasedg.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\nasedg\rkcomp.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_RKCOM=\
	".\Release\def_kintyp.mod"\
	".\Release\def_nasedg.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_RKCOM=\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_nasedg.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# End Group
# Begin Group "chemic"

# PROP Default_Filter "f90"
# Begin Source File

SOURCE=..\..\Sources\modules\chemic\Chemic.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CHEMI=\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CHEMI=\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\chemic\chm_begite.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CHM_B=\
	".\Release\def_chemic.mod"\
	".\Release\def_master.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CHM_B=\
	".\Debug\def_chemic.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\chemic\chm_begste.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CHM_BE=\
	".\Release\def_chemic.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CHM_BE=\
	".\Debug\def_chemic.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\chemic\chm_boumat.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CHM_BO=\
	".\Release\def_chemic.mod"\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CHM_BO=\
	".\Debug\def_chemic.mod"\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\chemic\chm_bouope.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CHM_BOU=\
	".\Release\def_chemic.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CHM_BOU=\
	".\Debug\def_chemic.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\chemic\chm_bouset.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CHM_BOUS=\
	".\Release\def_chemic.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CHM_BOUS=\
	".\Debug\def_chemic.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\chemic\chm_conblk.f90
# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\chemic\chm_concou.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CHM_C=\
	".\Release\def_chemic.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CHM_C=\
	".\Debug\def_chemic.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\chemic\chm_cvgunk.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CHM_CV=\
	".\Release\def_chemic.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CHM_CV=\
	".\Debug\def_chemic.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\chemic\chm_doiter.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CHM_D=\
	".\Release\def_chemic.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CHM_D=\
	".\Debug\def_chemic.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\chemic\chm_elmdir.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CHM_E=\
	".\Release\def_chemic.mod"\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CHM_E=\
	".\Debug\def_chemic.mod"\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\chemic\chm_elmgat.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CHM_EL=\
	".\Release\def_chemic.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CHM_EL=\
	".\Debug\def_chemic.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\chemic\chm_elmope.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CHM_ELM=\
	".\Release\def_chemic.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CHM_ELM=\
	".\Debug\def_chemic.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\chemic\chm_elmpro.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CHM_ELMP=\
	".\Release\def_chemic.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CHM_ELMP=\
	".\Debug\def_chemic.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\chemic\chm_endite.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CHM_EN=\
	".\Release\def_chemic.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CHM_EN=\
	".\Debug\def_chemic.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\chemic\chm_endste.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CHM_END=\
	".\Release\def_chemic.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CHM_END=\
	".\Debug\def_chemic.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\chemic\chm_grarea.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CHM_G=\
	".\Release\def_chemic.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CHM_G=\
	".\Debug\def_chemic.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\chemic\chm_iniunk.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CHM_I=\
	".\Release\def_chemic.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CHM_I=\
	".\Debug\def_chemic.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\chemic\chm_inivar.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CHM_IN=\
	".\Release\def_chemic.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CHM_IN=\
	".\Debug\def_chemic.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\chemic\chm_matrix.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CHM_M=\
	".\Release\def_chemic.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CHM_M=\
	".\Debug\def_chemic.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\chemic\chm_memall.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CHM_ME=\
	".\Release\def_chemic.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CHM_ME=\
	".\Debug\def_chemic.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\chemic\chm_membcs.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CHM_MEM=\
	".\Release\def_chemic.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CHM_MEM=\
	".\Debug\def_chemic.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\chemic\chm_memose.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CHM_MEMO=\
	".\Release\def_chemic.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_iofile.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CHM_MEMO=\
	".\Debug\def_chemic.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_iofile.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\chemic\chm_memphy.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CHM_MEMP=\
	".\Release\def_chemic.mod"\
	".\Release\def_master.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CHM_MEMP=\
	".\Debug\def_chemic.mod"\
	".\Debug\def_master.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\chemic\chm_newmsh.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CHM_N=\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CHM_N=\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\chemic\chm_odepro.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CHM_O=\
	".\Release\def_chemic.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CHM_O=\
	".\Debug\def_chemic.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\chemic\chm_openfi.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CHM_OP=\
	".\Release\def_chemic.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CHM_OP=\
	".\Debug\def_chemic.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\chemic\chm_outerr.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CHM_OU=\
	".\Release\def_chemic.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CHM_OU=\
	".\Debug\def_chemic.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\chemic\chm_outinf.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CHM_OUT=\
	".\Release\def_chemic.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CHM_OUT=\
	".\Debug\def_chemic.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\chemic\chm_output.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CHM_OUTP=\
	".\Release\def_chemic.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_iofile.mod"\
	".\Release\Mod_postpr.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CHM_OUTP=\
	".\Debug\def_chemic.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_iofile.mod"\
	".\Debug\Mod_postpr.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\chemic\chm_outset.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CHM_OUTS=\
	".\Release\def_chemic.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CHM_OUTS=\
	".\Debug\def_chemic.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\chemic\chm_outvar.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CHM_OUTV=\
	".\Release\def_chemic.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_postpr.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CHM_OUTV=\
	".\Debug\def_chemic.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_postpr.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\chemic\chm_outwit.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CHM_OUTW=\
	".\Release\def_chemic.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CHM_OUTW=\
	".\Debug\def_chemic.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\chemic\chm_parall.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CHM_P=\
	".\Release\def_chemic.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CHM_P=\
	".\Debug\def_chemic.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\chemic\chm_reabcs.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CHM_R=\
	".\Release\def_chemic.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CHM_R=\
	".\Debug\def_chemic.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\chemic\chm_reanut.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CHM_RE=\
	".\Release\def_chemic.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CHM_RE=\
	".\Debug\def_chemic.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\chemic\chm_reaous.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CHM_REA=\
	".\Release\def_chemic.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CHM_REA=\
	".\Debug\def_chemic.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\chemic\chm_reaphy.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CHM_REAP=\
	".\Release\def_chemic.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CHM_REAP=\
	".\Debug\def_chemic.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\chemic\chm_reapro.f90
# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\chemic\chm_restar.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CHM_RES=\
	".\Release\def_chemic.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CHM_RES=\
	".\Debug\def_chemic.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\chemic\chm_skymat.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CHM_S=\
	".\Release\def_chemic.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CHM_S=\
	".\Debug\def_chemic.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\chemic\chm_solite.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CHM_SO=\
	".\Release\def_chemic.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CHM_SO=\
	".\Debug\def_chemic.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\chemic\chm_timste.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CHM_T=\
	".\Release\def_chemic.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CHM_T=\
	".\Debug\def_chemic.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\chemic\chm_tistep.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CHM_TI=\
	".\Release\def_chemic.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CHM_TI=\
	".\Debug\def_chemic.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\chemic\chm_turnof.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CHM_TU=\
	".\Release\def_chemic.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CHM_TU=\
	".\Debug\def_chemic.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\chemic\chm_turnon.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CHM_TUR=\
	".\Release\def_chemic.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CHM_TUR=\
	".\Debug\def_chemic.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\chemic\chm_updbcs.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CHM_U=\
	".\Release\def_chemic.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CHM_U=\
	".\Debug\def_chemic.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\chemic\chm_updode.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CHM_UP=\
	".\Release\def_chemic.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CHM_UP=\
	".\Debug\def_chemic.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\chemic\chm_updtss.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CHM_UPD=\
	".\Release\def_chemic.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CHM_UPD=\
	".\Debug\def_chemic.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\chemic\chm_updunk.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CHM_UPDU=\
	".\Release\def_chemic.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CHM_UPDU=\
	".\Debug\def_chemic.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\chemic\chm_usrbcs.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CHM_US=\
	".\Release\def_chemic.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CHM_US=\
	".\Debug\def_chemic.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\chemic\chm_usrrea.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CHM_USR=\
	".\Release\def_chemic.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CHM_USR=\
	".\Debug\def_chemic.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\chemic\chm_velfun.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_CHM_V=\
	".\Release\def_chemic.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_CHM_V=\
	".\Debug\def_chemic.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\modules\chemic\def_chemic.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_DEF_CH=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_DEF_CH=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# End Group
# End Group
# Begin Group "services"

# PROP Default_Filter "f90"
# Begin Group "solmum"

# PROP Default_Filter "f90"
# Begin Source File

SOURCE=..\..\Sources\services\solmum\inimum.f90
# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\solmum\matspr.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_MATSP=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_MATSP=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\solmum\Solmum.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SOLMU=\
	"..\..\..\Mumps\include\dmumps_root.h"\
	"..\..\..\Mumps\include\dmumps_struc.h"\
	"..\..\..\Mumps\libseq\mpif.h"\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_solver.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SOLMU=\
	"..\..\..\Mumps\include\dmumps_root.h"\
	"..\..\..\Mumps\include\dmumps_struc.h"\
	"..\..\..\Mumps\libseq\mpif.h"\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_solver.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# End Group
# Begin Group "handfp"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\Sources\services\handfp\def_handfp.f90
# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\handfp\handfp.f90
# End Source File
# End Group
# Begin Group "parall"

# PROP Default_Filter "f90"
# Begin Source File

SOURCE=..\..\Sources\services\parall\def_parall.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_DEF_PAR=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_DEF_PAR=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\parall\mod_par_memchk.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_MOD_PA=\
	".\Release\def_master.mod"\
	".\Release\def_parall.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_MOD_PA=\
	".\Debug\def_master.mod"\
	".\Debug\def_parall.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\parall\par_alloca.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PAR_A=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parall.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PAR_A=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parall.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\parall\par_arrays.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PAR_AR=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parall.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PAR_AR=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parall.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\parall\par_barrie.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PAR_B=\
	"..\..\..\Mumps\libseq\mpif.h"\
	".\Release\def_master.mod"\
	".\Release\def_parall.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PAR_B=\
	"..\..\..\Mumps\libseq\mpif.h"\
	".\Debug\def_master.mod"\
	".\Debug\def_parall.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\parall\par_bounda.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PAR_BO=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_parall.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PAR_BO=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_parall.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\parall\par_broadc.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PAR_BR=\
	"..\..\..\Mumps\libseq\mpif.h"\
	".\Release\def_master.mod"\
	".\Release\def_parall.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PAR_BR=\
	"..\..\..\Mumps\libseq\mpif.h"\
	".\Debug\def_master.mod"\
	".\Debug\def_parall.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\parall\par_chkpoi.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PAR_C=\
	".\Release\def_master.mod"\
	".\Release\def_parall.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PAR_C=\
	".\Debug\def_master.mod"\
	".\Debug\def_parall.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\parall\par_colgra.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PAR_CO=\
	".\Release\def_parall.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PAR_CO=\
	".\Debug\def_parall.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\parall\par_commun.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PAR_COM=\
	".\Release\def_master.mod"\
	".\Release\def_parall.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	".\Release\mod_par_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PAR_COM=\
	".\Debug\def_master.mod"\
	".\Debug\def_parall.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	".\Debug\mod_par_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\parall\par_comset.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PAR_COMS=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parall.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PAR_COMS=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parall.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\parall\par_cregro.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PAR_CR=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parall.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	".\Release\Mod_memchk.mod"\
	".\Release\Mod_postpr.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PAR_CR=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parall.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	".\Debug\Mod_memchk.mod"\
	".\Debug\Mod_postpr.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\parall\par_deallo.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PAR_D=\
	".\Release\def_master.mod"\
	".\Release\def_parall.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PAR_D=\
	".\Debug\def_master.mod"\
	".\Debug\def_parall.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\parall\par_disbou.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PAR_DI=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_parall.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PAR_DI=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_parall.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\parall\par_domgra.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PAR_DO=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_parall.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PAR_DO=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_parall.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\parall\par_domlis.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PAR_DOM=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PAR_DOM=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\parall\par_duagra.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PAR_DU=\
	".\Release\def_parall.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PAR_DU=\
	".\Debug\def_parall.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\parall\par_elmgra.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PAR_E=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parall.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PAR_E=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parall.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\parall\par_errors.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PAR_ER=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parall.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PAR_ER=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parall.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\parall\par_exampl.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PAR_EX=\
	".\Release\def_domain.mod"\
	".\Release\def_parall.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PAR_EX=\
	".\Debug\def_domain.mod"\
	".\Debug\def_parall.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\parall\par_extnor.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PAR_EXT=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parall.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PAR_EXT=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parall.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\parall\par_filnam.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PAR_F=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parall.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PAR_F=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parall.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\parall\par_finali.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PAR_FI=\
	"..\..\..\Mumps\libseq\mpif.h"\
	".\Release\def_master.mod"\
	".\Release\def_parall.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PAR_FI=\
	"..\..\..\Mumps\libseq\mpif.h"\
	".\Debug\def_master.mod"\
	".\Debug\def_parall.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\parall\par_gather.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PAR_G=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parall.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PAR_G=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parall.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\parall\par_inidat.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PAR_I=\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parall.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PAR_I=\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parall.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\parall\par_initia.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PAR_IN=\
	"..\..\..\Mumps\libseq\mpif.h"\
	".\Release\def_master.mod"\
	".\Release\def_parall.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PAR_IN=\
	"..\..\..\Mumps\libseq\mpif.h"\
	".\Debug\def_master.mod"\
	".\Debug\def_parall.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\parall\par_livinf.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PAR_L=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parall.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PAR_L=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parall.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\parall\par_locnum.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PAR_LO=\
	"..\..\..\Mumps\libseq\mpif.h"\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parall.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PAR_LO=\
	"..\..\..\Mumps\libseq\mpif.h"\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parall.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\parall\par_memory.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PAR_M=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parall.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PAR_M=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parall.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\parall\par_memset.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PAR_ME=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parall.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PAR_ME=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parall.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\parall\par_metis.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PAR_MET=\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PAR_MET=\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\parall\par_openfi.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PAR_O=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parall.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PAR_O=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parall.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\parall\par_operat.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PAR_OP=\
	"..\..\..\Mumps\libseq\mpif.h"\
	".\Release\def_master.mod"\
	".\Release\def_parall.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PAR_OP=\
	"..\..\..\Mumps\libseq\mpif.h"\
	".\Debug\def_master.mod"\
	".\Debug\def_parall.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\parall\par_outcpu.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PAR_OU=\
	".\Release\def_master.mod"\
	".\Release\def_parall.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PAR_OU=\
	".\Debug\def_master.mod"\
	".\Debug\def_parall.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\parall\par_outinf.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PAR_OUT=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parall.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PAR_OUT=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parall.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\parall\par_outprt.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PAR_OUTP=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parall.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PAR_OUTP=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parall.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\parall\par_partit.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PAR_P=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parall.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PAR_P=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parall.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\parall\par_period.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PAR_PE=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PAR_PE=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\parall\par_pospar.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PAR_PO=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parall.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PAR_PO=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parall.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\parall\par_reapro.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PAR_R=\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parall.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PAR_R=\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parall.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\parall\par_receiv.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PAR_RE=\
	"..\..\..\Mumps\libseq\mpif.h"\
	".\Release\def_master.mod"\
	".\Release\def_parall.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PAR_RE=\
	"..\..\..\Mumps\libseq\mpif.h"\
	".\Debug\def_master.mod"\
	".\Debug\def_parall.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\parall\par_scatt2.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PAR_S=\
	"..\..\..\Mumps\libseq\mpif.h"\
	".\Release\def_master.mod"\
	".\Release\def_parall.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PAR_S=\
	"..\..\..\Mumps\libseq\mpif.h"\
	".\Debug\def_master.mod"\
	".\Debug\def_parall.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\parall\par_scatte.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PAR_SC=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parall.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PAR_SC=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parall.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\parall\par_sencom.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PAR_SE=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parall.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	".\Release\mod_par_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PAR_SE=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parall.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	".\Debug\mod_par_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\parall\par_sendat.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PAR_SEN=\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parall.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PAR_SEN=\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parall.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\parall\par_sendin.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PAR_SEND=\
	"..\..\..\Mumps\libseq\mpif.h"\
	".\Release\def_master.mod"\
	".\Release\def_parall.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PAR_SEND=\
	"..\..\..\Mumps\libseq\mpif.h"\
	".\Debug\def_master.mod"\
	".\Debug\def_parall.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\parall\par_sengeo.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PAR_SENG=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parall.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PAR_SENG=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parall.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\parall\par_senset.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PAR_SENS=\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parall.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PAR_SENS=\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parall.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\parall\par_slexch.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PAR_SL=\
	"..\..\..\Mumps\libseq\mpif.h"\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parall.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PAR_SL=\
	"..\..\..\Mumps\libseq\mpif.h"\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parall.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\parall\par_solpls.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PAR_SO=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parall.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PAR_SO=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parall.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\parall\par_subgra.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PAR_SU=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PAR_SU=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\parall\par_turnof.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PAR_T=\
	".\Release\def_master.mod"\
	".\Release\def_parall.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PAR_T=\
	".\Debug\def_master.mod"\
	".\Debug\def_parall.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\parall\par_volume.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PAR_V=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parall.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PAR_V=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parall.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\parall\Parall.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PARAL=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parall.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_postpr.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PARAL=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parall.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_postpr.mod"\
	

!ENDIF 

# End Source File
# End Group
# Begin Group "cgns24"

# PROP Default_Filter "f90;c"
# Begin Source File

SOURCE=..\..\Sources\services\cgns24\fak_cgns24.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_FAK_C=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	".\Release\Mod_postpr.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_FAK_C=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	".\Debug\Mod_postpr.mod"\
	

!ENDIF 

# End Source File
# End Group
# Begin Group "dodeme"

# PROP Default_Filter "f90"
# Begin Source File

SOURCE=..\..\Sources\services\dodeme\def_dodeme.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_DEF_DO=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_DEF_DO=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\dodeme\dod_assbcs.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_DOD_A=\
	".\Release\def_dodeme.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_DOD_A=\
	".\Debug\def_dodeme.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\dodeme\dod_assmat.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_DOD_AS=\
	".\Release\def_kintyp.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_DOD_AS=\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\dodeme\dod_bouope.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_DOD_B=\
	".\Release\def_dodeme.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_DOD_B=\
	".\Debug\def_dodeme.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\dodeme\dod_chimer.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_DOD_C=\
	".\Release\def_dodeme.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_DOD_C=\
	".\Debug\def_dodeme.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\dodeme\dod_defint.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_DOD_D=\
	".\Release\def_dodeme.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_DOD_D=\
	".\Debug\def_dodeme.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\dodeme\dod_domain.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_DOD_DO=\
	".\Release\def_dodeme.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_DOD_DO=\
	".\Debug\def_dodeme.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\dodeme\dod_domgra.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_DOD_DOM=\
	".\Release\def_dodeme.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_DOD_DOM=\
	".\Debug\def_dodeme.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\dodeme\dod_elsest.f
# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\dodeme\dod_embedd.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_DOD_E=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_DOD_E=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\dodeme\dod_finelm.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_DOD_F=\
	".\Release\def_dodeme.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_DOD_F=\
	".\Debug\def_dodeme.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\dodeme\dod_gaucon.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_DOD_G=\
	".\Release\def_dodeme.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_DOD_G=\
	".\Debug\def_dodeme.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\dodeme\dod_gra2nd.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_DOD_GR=\
	".\Release\def_dodeme.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\mod_dod_memchk.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_DOD_GR=\
	".\Debug\def_dodeme.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\mod_dod_memchk.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\dodeme\dod_initia.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_DOD_I=\
	".\Release\def_dodeme.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_DOD_I=\
	".\Debug\def_dodeme.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\dodeme\dod_inivar.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_DOD_IN=\
	".\Release\def_dodeme.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_DOD_IN=\
	".\Debug\def_dodeme.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\dodeme\dod_livinf.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_DOD_L=\
	".\Release\def_dodeme.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_DOD_L=\
	".\Debug\def_dodeme.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\dodeme\dod_memall.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_DOD_M=\
	".\Release\def_dodeme.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solver.mod"\
	".\Release\mod_dod_memchk.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_DOD_M=\
	".\Debug\def_dodeme.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solver.mod"\
	".\Debug\mod_dod_memchk.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\dodeme\dod_nodcon.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_DOD_N=\
	".\Release\def_dodeme.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_kintyp.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_DOD_N=\
	".\Debug\def_dodeme.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_kintyp.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\dodeme\dod_openfi.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_DOD_O=\
	".\Release\def_dodeme.mod"\
	".\Release\def_master.mod"\
	".\Release\def_postpr.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_DOD_O=\
	".\Debug\def_dodeme.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_postpr.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\dodeme\dod_outcpu.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_DOD_OU=\
	".\Release\def_dodeme.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_DOD_OU=\
	".\Debug\def_dodeme.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\dodeme\dod_outinf.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_DOD_OUT=\
	".\Release\def_dodeme.mod"\
	".\Release\def_master.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_DOD_OUT=\
	".\Debug\def_dodeme.mod"\
	".\Debug\def_master.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\dodeme\dod_readat.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_DOD_R=\
	".\Release\def_dodeme.mod"\
	".\Release\def_domain.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_DOD_R=\
	".\Debug\def_dodeme.mod"\
	".\Debug\def_domain.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\dodeme\dod_reapro.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_DOD_RE=\
	".\Release\def_dodeme.mod"\
	".\Release\def_inpout.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_DOD_RE=\
	".\Debug\def_dodeme.mod"\
	".\Debug\def_inpout.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\dodeme\dod_transf.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_DOD_T=\
	".\Release\def_dodeme.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_DOD_T=\
	".\Debug\def_dodeme.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\dodeme\dod_turnof.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_DOD_TU=\
	".\Release\def_dodeme.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_DOD_TU=\
	".\Debug\def_dodeme.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\dodeme\dod_turnon.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_DOD_TUR=\
	".\Release\def_dodeme.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\Mod_iofile.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_DOD_TUR=\
	".\Debug\def_dodeme.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\Mod_iofile.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\dodeme\Dodeme.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_DODEM=\
	".\Release\def_dodeme.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_DODEM=\
	".\Debug\def_dodeme.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\dodeme\mod_dod_memchk.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_MOD_D=\
	".\Release\def_dodeme.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_MOD_D=\
	".\Debug\def_dodeme.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	

!ENDIF 

# End Source File
# End Group
# Begin Group "solpls"

# PROP Default_Filter "f90"
# Begin Source File

SOURCE=..\..\Sources\services\solpls\def_solpls.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_DEF_SOL=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_DEF_SOL=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\solpls\pls_csrase.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PLS_C=\
	".\Release\def_domain.mod"\
	".\Release\def_elmope.mod"\
	".\Release\def_master.mod"\
	".\Release\def_solpls.mod"\
	".\Release\def_solver.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PLS_C=\
	".\Debug\def_domain.mod"\
	".\Debug\def_elmope.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_solpls.mod"\
	".\Debug\def_solver.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\solpls\pls_driver.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PLS_D=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_solpls.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PLS_D=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_solpls.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\solpls\pls_elmtyp.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PLS_E=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_solpls.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PLS_E=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_solpls.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\solpls\pls_graphs.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PLS_G=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_solpls.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PLS_G=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_solpls.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\solpls\pls_memory.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_PLS_M=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solpls.mod"\
	".\Release\Mod_memchk.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_PLS_M=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solpls.mod"\
	".\Debug\Mod_memchk.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\services\solpls\Solpls.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_SOLPL=\
	".\Release\def_domain.mod"\
	".\Release\def_master.mod"\
	".\Release\def_parame.mod"\
	".\Release\def_solpls.mod"\
	".\Release\Mod_postpr.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_SOLPL=\
	".\Debug\def_domain.mod"\
	".\Debug\def_master.mod"\
	".\Debug\def_parame.mod"\
	".\Debug\def_solpls.mod"\
	".\Debug\Mod_postpr.mod"\
	

!ENDIF 

# End Source File
# End Group
# Begin Group "adapti"

# PROP Default_Filter "f90"
# End Group
# End Group
# Begin Group "Automatic"

# PROP Default_Filter "f90"
# Begin Group "fakes"

# PROP Default_Filter "f90"
# Begin Source File

SOURCE=..\..\Sources\automatic\fakes\fak_adapti.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_FAK_A=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_FAK_A=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\Sources\automatic\fakes\fak_latbol.f90

!IF  "$(CFG)" == "Alya - Win32 Release"

DEP_F90_FAK_L=\
	".\Release\def_kintyp.mod"\
	

!ELSEIF  "$(CFG)" == "Alya - Win32 Debug"

DEP_F90_FAK_L=\
	".\Debug\def_kintyp.mod"\
	

!ENDIF 

# End Source File
# End Group
# Begin Group "info"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\Sources\automatic\info\infmak.f90
# End Source File
# Begin Source File

SOURCE=..\..\Sources\automatic\info\infsvn.f90
# End Source File
# End Group
# End Group
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl;fi;fd"
# End Group
# Begin Group "Resource Files"

# PROP Default_Filter "ico;cur;bmp;dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe"
# End Group
# End Target
# End Project
