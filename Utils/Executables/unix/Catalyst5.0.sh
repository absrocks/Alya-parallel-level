#!/bin/bash
module purge
module load MKL/11.0.1

gcc  \
-I/apps/PARAVIEW/SRC/5.0.0_catalyst/ParaViewSuperbuild-build/install/include \
-I/apps/PARAVIEW/SRC/5.0.0_catalyst/ParaViewSuperbuild-build/paraview/src/paraview/CoProcessing/Catalyst \
-I/apps/PARAVIEW/SRC/5.0.0_catalyst/ParaViewSuperbuild-build/paraview/src/paraview/CoProcessing/PythonCatalyst \
-I/apps/PARAVIEW/SRC/5.0.0_catalyst/ParaViewSuperbuild-build/paraview/src/paraview/VTK/Common/Core \
-I/apps/PARAVIEW/SRC/5.0.0_catalyst/ParaViewSuperbuild-build/paraview/src/paraview/VTK/Utilities/KWIML \
-I/apps/PARAVIEW/SRC/5.0.0_catalyst/ParaViewSuperbuild-build/paraview/src/paraview/VTK/Common/DataModel \
-I/apps/PARAVIEW/SRC/5.0.0_catalyst/ParaViewSuperbuild-build/paraview/src/paraview/VTK/Filters/Extraction \
-I/apps/PARAVIEW/SRC/5.0.0_catalyst/ParaViewSuperbuild-build/paraview/src/paraview/VTK/Filters/Core \
-I/apps/PARAVIEW/SRC/5.0.0_catalyst/ParaViewSuperbuild-build/paraview/src/paraview/VTK/Filters/General \
-I/apps/PARAVIEW/SRC/5.0.0_catalyst/ParaViewSuperbuild-build/paraview/src/paraview/VTK/Filters/Statistics \
-I/apps/PARAVIEW/SRC/5.0.0_catalyst/ParaViewSuperbuild-build/paraview/src/paraview/VTK/Filters/Parallel \
-I/apps/PARAVIEW/SRC/5.0.0_catalyst/ParaViewSuperbuild-build/paraview/src/paraview/VTK/Filters/Geometry \
-I/apps/PARAVIEW/SRC/5.0.0_catalyst/ParaViewSuperbuild-build/paraview/src/paraview/VTK/Filters/Modeling \
-I/apps/PARAVIEW/SRC/5.0.0_catalyst/ParaViewSuperbuild-build/paraview/src/paraview/VTK/Filters/Sources \
-I/apps/PARAVIEW/SRC/5.0.0_catalyst/ParaViewSuperbuild-build/paraview/src/paraview/VTK/Rendering/Core \
-I/apps/PARAVIEW/SRC/5.0.0_catalyst/ParaViewSuperbuild-build/paraview/src/paraview/ParaViewCore/ServerManager/SMApplication \
-I/apps/PARAVIEW/SRC/5.0.0_catalyst/ParaViewSuperbuild-build/paraview/src/paraview/ParaViewCore/ServerManager/Core \
-I/apps/PARAVIEW/SRC/5.0.0_catalyst/ParaViewSuperbuild-build/paraview/src/paraview/ParaViewCore/ServerImplementation/Core \
-I/apps/PARAVIEW/SRC/5.0.0_catalyst/ParaViewSuperbuild-build/paraview/src/paraview/ParaViewCore/ClientServerCore/Core \
-I/apps/PARAVIEW/SRC/5.0.0_catalyst/ParaViewSuperbuild-build/paraview/src/paraview/ParaViewCore/VTKExtensions/Core \
-I/apps/PARAVIEW/SRC/5.0.0_catalyst/ParaViewSuperbuild-build/paraview/src/paraview-build/CoProcessing/Catalyst \
-I/apps/PARAVIEW/SRC/5.0.0_catalyst/ParaViewSuperbuild-build/paraview/src/paraview-build/CoProcessing/PythonCatalyst \
-I/apps/PARAVIEW/SRC/5.0.0_catalyst/ParaViewSuperbuild-build/paraview/src/paraview-build/VTK/Common/Core \
-I/apps/PARAVIEW/SRC/5.0.0_catalyst/ParaViewSuperbuild-build/paraview/src/paraview-build/VTK/Utilities/KWIML \
-I/apps/PARAVIEW/SRC/5.0.0_catalyst/ParaViewSuperbuild-build/paraview/src/paraview-build/VTK/Common/DataModel \
-I/apps/PARAVIEW/SRC/5.0.0_catalyst/ParaViewSuperbuild-build/paraview/src/paraview-build/VTK/Filters/Extraction \
-I/apps/PARAVIEW/SRC/5.0.0_catalyst/ParaViewSuperbuild-build/paraview/src/paraview-build/VTK/Filters/Core \
-I/apps/PARAVIEW/SRC/5.0.0_catalyst/ParaViewSuperbuild-build/paraview/src/paraview-build/VTK/Filters/General \
-I/apps/PARAVIEW/SRC/5.0.0_catalyst/ParaViewSuperbuild-build/paraview/src/paraview-build/VTK/Filters/Statistics \
-I/apps/PARAVIEW/SRC/5.0.0_catalyst/ParaViewSuperbuild-build/paraview/src/paraview-build/VTK/Filters/Parallel \
-I/apps/PARAVIEW/SRC/5.0.0_catalyst/ParaViewSuperbuild-build/paraview/src/paraview-build/VTK/Filters/Geometry \
-I/apps/PARAVIEW/SRC/5.0.0_catalyst/ParaViewSuperbuild-build/paraview/src/paraview-build/VTK/Filters/Modeling \
-I/apps/PARAVIEW/SRC/5.0.0_catalyst/ParaViewSuperbuild-build/paraview/src/paraview-build/VTK/Filters/Sources \
-I/apps/PARAVIEW/SRC/5.0.0_catalyst/ParaViewSuperbuild-build/paraview/src/paraview-build/VTK/Rendering/Core \
-I/apps/PARAVIEW/SRC/5.0.0_catalyst/ParaViewSuperbuild-build/paraview/src/paraview-build/ParaViewCore/ServerManager/SMApplication \
-I/apps/PARAVIEW/SRC/5.0.0_catalyst/ParaViewSuperbuild-build/paraview/src/paraview-build/ParaViewCore/ServerManager/Core \
-I/apps/PARAVIEW/SRC/5.0.0_catalyst/ParaViewSuperbuild-build/paraview/src/paraview-build/ParaViewCore/ServerImplementation/Core \
-I/apps/PARAVIEW/SRC/5.0.0_catalyst/ParaViewSuperbuild-build/paraview/src/paraview-build/ParaViewCore/ClientServerCore/Core \
-I/apps/PARAVIEW/SRC/5.0.0_catalyst/ParaViewSuperbuild-build/paraview/src/paraview-build/ParaViewCore/VTKExtensions/Core \
-c ../../Thirdparties/Catalyst/FECxxAdaptor.cxx -o  ../../Thirdparties/Catalyst/FECxxAdaptor.o 

module load intel/15.0.2
ifort -c -O3 ../../Thirdparties/Catalyst/FEFortranAdaptor.f90 -I ../../Executables/unix/Objects_x/ -o ../../Thirdparties/Catalyst/FEFortranAdaptor.o
