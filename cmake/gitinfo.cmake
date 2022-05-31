# cmake/gitversion.cmake
message(STATUS "Resolving GIT Version")
 
set(_git_remote "unknown")
set(_git_revision "unknown")
set(_git_branch "unknown")
set(_git_dirty "unknown")

message(STATUS "Looking for ${CMAKE_SOURCE_DIR}.git")

if(EXISTS ${root_dir}/.git)

  find_package(Git)

  if(GIT_FOUND)
    execute_process(
      COMMAND ${GIT_EXECUTABLE} remote get-url origin
      WORKING_DIRECTORY "${local_dir}"
      OUTPUT_VARIABLE _git_remote
      OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    
    execute_process(
      COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
      WORKING_DIRECTORY "${local_dir}"
      OUTPUT_VARIABLE _git_revision
      OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    
    execute_process(
      COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
      WORKING_DIRECTORY "${local_dir}"
      OUTPUT_VARIABLE _git_branch
      OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    
  #  execute_process(
  #    COMMAND ${GIT_EXECUTABLE} diff HEAD --quiet && echo $?
  #    WORKING_DIRECTORY "${local_dir}"
  #    OUTPUT_VARIABLE _git_dirty
  #    OUTPUT_STRIP_TRAILING_WHITESPACE
  #  )
    message(STATUS "GIT remote: ${_git_remote}")
    message(STATUS "GIT revision: ${_git_revision}")
    message(STATUS "GIT branch: ${_git_branch}")
  #  message(STATUS "GIT dirty: ${_git_dirty}")
  else()
    message(STATUS "GIT not found")
  endif()
else()
  message(STATUS "This is not a git repository")
endif() 
string(TIMESTAMP _time_stamp)

configure_file(${local_dir}/gitinfo.f90.in ${output_dir}/gitinfo.f90 @ONLY)

