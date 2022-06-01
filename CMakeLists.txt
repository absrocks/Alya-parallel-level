cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
find_package(PkgConfig REQUIRED)

project(alya-project)

set(CMAKE_VERBOSE_MAKEFILE OFF)

#----------Options---------

set(CONF "auto" CACHE STRING "Configuration preset loaded by cmake")
set_property(CACHE CONF PROPERTY STRINGS none custom auto)

#Integer size
set(INTEGER_SIZE "4" CACHE STRING "Integer size (in bytes)")
set_property(CACHE INTEGER_SIZE PROPERTY STRINGS 4 8)

#Optimization level
set(OPTIMIZATION_LEVEL "3" CACHE STRING "Optimization level")
set_property(CACHE OPTIMIZATION_LEVEL PROPERTY STRINGS 0 1 2 3)

#NDIMEPAR
option(WITH_NDIMEPAR "Enable NDIMEPAR" OFF)

#STD2008
option(WITH_STD2008 "Enable fortran 2008 standard compilation" ON)

#IPO
option(WITH_IPO "Enable ipo flag" OFF)

#VECTOR_SIZE
set(VECTOR_SIZE "16" CACHE STRING "Vector size")

#Metis
set(WITH_METIS "auto" CACHE STRING "Enable metis")
set_property(CACHE WITH_METIS PROPERTY STRINGS auto OFF)

#COMMDOM
set(WITH_COMMDOM "OFF" CACHE STRING "Enable commdom")
set_property(CACHE WITH_COMMDOM PROPERTY STRINGS auto OFF)

#GMSH
set(WITH_GMSH "OFF" CACHE STRING "Enable bridge to gmsh")
set_property(CACHE WITH_GMSH PROPERTY STRINGS auto OFF)

#CANTERA
set(WITH_CANTERA "OFF" CACHE STRING "Enable bridge to cantera")
set_property(CACHE WITH_CANTERA PROPERTY STRINGS external OFF)

set(CANTERA_DIR "" CACHE PATH "Cantera installation directory")

#TALP
set(WITH_DLB_TALP "OFF" CACHE STRING "Enable TALP")
set_property(CACHE WITH_DLB_TALP PROPERTY STRINGS external OFF)

set(DLB_DIR "" CACHE PATH "DLB installation directory")

#Shared library
option(SHARED "Build shared libraries" OFF)

#MPI
option(WITH_MPI "Enable MPI" ON)

#OPENMP
option(WITH_OPENMP "Enable OPENMP" OFF)

#GPU
option(WITH_GPU "Enable compilation for GPU" OFF)

#Code coverage
option(WITH_CODE_COVERAGE "Compile using code coverage flags" OFF)

#Include directory
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

#Modules
include(subdir)
SUBDIRLIST(MODULES "${CMAKE_SOURCE_DIR}/src/modules")
foreach(MOD ${MODULES})
  string(TOUPPER ${MOD} UPMOD)
  option(WITH_MODULE_${UPMOD} "Enable ${MOD} module" ON)
endforeach()

#----------CMake configuration---------

#Build type
include(buildType)

if(CONF MATCHES "none")
#nothing
elseif(CONF MATCHES "custom")
  include("${CMAKE_SOURCE_DIR}/config.cmake")
else()
  if (CMAKE_C_COMPILER_ID STREQUAL "GNU")
    # using GNU
    include("${CMAKE_SOURCE_DIR}/config/gnu.cmake")
  elseif (CMAKE_C_COMPILER_ID STREQUAL "Clang")
    # using Clang
    message(WARNING "Default compiler is Clang, but no Clang configuration is available. Switching to GNU configuration. If you have some issues, try to define manually CMAKE_C_COMPILER and CMAKE_CXX_COMPILER to gcc and g++ paths")
    include("${CMAKE_SOURCE_DIR}/config/gnu.cmake")
  elseif (CMAKE_C_COMPILER_ID STREQUAL "Intel")
    # using Intel
    include("${CMAKE_SOURCE_DIR}/config/intel.cmake")
  elseif (CMAKE_C_COMPILER_ID STREQUAL "PGI")
    # using PGI
    if (NOT WITH_GPU)
      include("${CMAKE_SOURCE_DIR}/config/pgi.cmake")
    else()
      include("${CMAKE_SOURCE_DIR}/config/pgi-gpu.cmake")
    endif()
  elseif (CMAKE_C_COMPILER_ID STREQUAL "XL")
    # using XL
    include("${CMAKE_SOURCE_DIR}/config/xl.cmake")
  else()
    message(FATAL_ERROR "No preset found for your compiler, use a custom configuration!")
  endif()
endif()

#Metis management
include(metis)
init_metis()

#Cantera management
include(cantera)
init_cantera()

#TALP management
include(talp)
init_talp()

#MPI Management
include(mpi)
init_mpi()

if(SHARED)
  set(LIBRARY_TYPE SHARED)
else()
  set(LIBRARY_TYPE STATIC)
endif(SHARED)

#if (GPROF_OPT)
#    set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-pg")
#    set(CMAKE_C_FLAGS ${CMAKE_C_FLAGS} "-pg")
#endif()

# use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# the RPATH to be used when installing, but only if it's not a system directory
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
     set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif("${isSystemDir}" STREQUAL "-1")

enable_testing()

ADD_SUBDIRECTORY(deps)
ADD_SUBDIRECTORY(external)
ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(unitt)
ADD_SUBDIRECTORY(tests)

set(DOCUMENTATION_DIR "${CMAKE_SOURCE_DIR}/doc/markdown")
set(DOCUMENTATION_BIN_DIR "${CMAKE_BINARY_DIR}/doc")
set(MARKDOWN_BIN_DIR "${DOCUMENTATION_BIN_DIR}/md")
set(WIKI_BIN_DIR "${DOCUMENTATION_BIN_DIR}/html")
set(DOCUMENTATION_SRC_DIR "${CMAKE_SOURCE_DIR}/Sources")

add_custom_target(uninstall "${CMAKE_COMMAND}" -P "${CMAKE_MODULE_PATH}/uninstall.cmake")
add_custom_target(documentation ./generate-md "${DOCUMENTATION_SRC_DIR}" "${MARKDOWN_BIN_DIR}" WORKING_DIRECTORY "${DOCUMENTATION_DIR}")
add_custom_target(wiki ./generate-html "${MARKDOWN_BIN_DIR}" "${WIKI_BIN_DIR}" WORKING_DIRECTORY "${DOCUMENTATION_DIR}")
add_dependencies(wiki documentation)

