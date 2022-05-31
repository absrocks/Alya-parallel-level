function(target_link_libraries_openmp)
  if (WITH_OPENMP)
    FIND_PACKAGE(OpenMP_Fortran)
    SET_TARGET_PROPERTIES(${PROJECT_NAME} PROPERTIES
                        COMPILE_FLAGS "${OpenMP_Fortran_FLAGS}"
                        LINK_FLAGS "${OpenMP_Fortran_FLAGS}")
  endif()
endfunction()
