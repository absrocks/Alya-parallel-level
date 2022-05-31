# Add this directory for internal users.
include_directories(.)
# Find sources.
file(GLOB metis_sources *.c)
add_definitions(-DIDXTYPEWIDTH=${METIS_INTEGER})
add_definitions(-DREALTYPEWIDTH=${METIS_REAL})

# Build libmetis.
add_library(metis ${METIS_LIBRARY_TYPE} ${GKlib_sources} ${metis_sources})
if(UNIX)
  target_link_libraries(metis m)
endif()

if(METIS_INSTALL)
  install(TARGETS metis
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION lib
    ARCHIVE DESTINATION lib)
endif()
