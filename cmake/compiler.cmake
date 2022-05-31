function(set_flags)
  #Fortran Management

  if (DEFINED CUSTOM_Fortran_FLAGS)
    set(APPEND_Fortran_FLAGS "${CUSTOM_Fortran_FLAGS}")
  endif()
  
  #Standard 2008
  if (WITH_STD2008)
    set(APPEND_Fortran_FLAGS "${APPEND_Fortran_FLAGS} ${CUSTOM_Fortran_FLAGS_STD2008}")
  endif()
 
  #Code coverage
  if (WITH_CODE_COVERAGE)
    set(APPEND_Fortran_FLAGS "${APPEND_Fortran_FLAGS} ${CUSTOM_Fortran_FLAGS_CODE_COVERAGE}")
  endif()

  if (DEFINED APPEND_Fortran_FLAGS)
    #Fortran compilation flag management
    string(REPLACE " " ";" REPLACED_Fortran_FLAGS ${APPEND_Fortran_FLAGS})
    target_compile_options(${PROJECT_NAME} PUBLIC $<$<COMPILE_LANGUAGE:Fortran>:${REPLACED_Fortran_FLAGS}>)
  endif()

  #Release flag management
  if (DEFINED CUSTOM_Fortran_FLAGS_RELEASE)
    set(CMAKE_Fortran_FLAGS_RELEASE "${CUSTOM_Fortran_FLAGS_RELEASE}" PARENT_SCOPE)
  endif()

  #Optimizations
  if (DEFINED CUSTOM_Fortran_FLAGS_O${OPTIMIZATION_LEVEL})
    set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${CUSTOM_Fortran_FLAGS_O${OPTIMIZATION_LEVEL}}" PARENT_SCOPE)
  endif()
  
  #Architecture
  if (DEFINED CUSTOM_Fortran_FLAGS_ARCHITECTURE)
    set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${CUSTOM_Fortran_FLAGS_ARCHITECTURE}" PARENT_SCOPE)
  endif()
  
  if (WITH_IPO)
    #Architecture
    if (DEFINED CUSTOM_Fortran_FLAGS_IPO)
      set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${CUSTOM_Fortran_FLAGS_IPO}" PARENT_SCOPE)
    endif()
  endif()
  
  #Debug flag management
  if (DEFINED CUSTOM_Fortran_FLAGS_DEBUG)
    set(CMAKE_Fortran_FLAGS_DEBUG "${CUSTOM_Fortran_FLAGS_DEBUG}" PARENT_SCOPE)
  endif()

  #C Management
  
  if (DEFINED CUSTOM_C_FLAGS)
    set(APPEND_C_FLAGS "${CUSTOM_C_FLAGS}")
    string(REPLACE " " ";" REPLACED_C_FLAGS ${APPEND_C_FLAGS})
    target_compile_options(${PROJECT_NAME} PUBLIC $<$<COMPILE_LANGUAGE:C>:${REPLACED_C_FLAGS}>)
  endif()

  #C++ Management

  if (DEFINED CUSTOM_CXX_FLAGS)
    set(APPEND_CXX_FLAGS "${CUSTOM_CXX_FLAGS}")
    string(REPLACE " " ";" REPLACED_CXX_FLAGS ${APPEND_CXX_FLAGS})
    target_compile_options(${PROJECT_NAME} PUBLIC $<$<COMPILE_LANGUAGE:CXX>:${REPLACED_CXX_FLAGS}>)
  endif()

endfunction()

function(set_defs)
  #CMake flag
  add_definitions(-DCMAKE)
  
  #MPI flag
  if (NOT WITH_MPI)
    add_definitions(-DMPI_OFF)
  endif()

  #NDIMEPAR flag
  if (WITH_NDIMEPAR)
    add_definitions(-DNDIMEPAR)
  endif()

  #Vector size flag
  if (VECTOR_SIZE MATCHES "0")
  #nothing
  else()
    add_definitions(-DVECTOR_SIZE=${VECTOR_SIZE})
  endif()

  #Integer management
  if (INTEGER_SIZE MATCHES "4")
    set(CMAKE_Fortran_FLAGS "${CUSTOM_Fortran_FLAGS_I4}" PARENT_SCOPE)
  elseif (INTEGER_SIZE MATCHES "8")
    add_definitions(-DI8)
    set(CMAKE_Fortran_FLAGS "${CUSTOM_Fortran_FLAGS_I8}" PARENT_SCOPE)
  endif()
endfunction()

function(set_openmpi)
 include(mpi)
 include_directories_mpi()
 target_link_libraries_mpi()
 include(openmp)
 target_link_libraries_openmp()
endfunction()

function(set_mpi)
 include(mpi)
 include_directories_mpi()
 target_link_libraries_mpi()
endfunction()

function(set_mpi_c)
 include(mpi)
 include_directories_mpi_c()
 target_link_libraries_mpi_c()
endfunction()

function(set_mpi_cxx)
 include(mpi)
 include_directories_mpi_cxx()
 target_link_libraries_mpi_cxx()
endfunction()
