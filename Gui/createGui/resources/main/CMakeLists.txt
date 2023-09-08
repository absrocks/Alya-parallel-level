
cmake_minimum_required(VERSION 2.6)

project(ParaViewAlya)

FIND_PACKAGE(ParaView REQUIRED)
FIND_PACKAGE(Qt4)
INCLUDE(${PARAVIEW_USE_FILE})

INCLUDE(${ParaView_CMAKE_DIR}/ParaViewBranding.cmake)
INCLUDE(${ParaView_CMAKE_DIR}/ParaViewPlugins.cmake)

#------------------------------------------------------------------------------
# Add extra library containing custom code for the client.
QT4_WRAP_CPP(MOC_BUILT_SOURCES
  myMainWindow.h
  alyawidget.h
  <<hModuleFiles>>)

QT4_WRAP_UI(UI_BUILT_SOURCES
  myMainWindow.ui
  alyawidget.ui
  <<uiModuleFiles>>
)

##########################################################################
## Add resources (icons etc.)
SET (UI_RESOURCES
  alya.qrc
)

QT4_ADD_RESOURCES(RCS_BUILT_SOURCES
  ${UI_RESOURCES}
)

SOURCE_GROUP("Resources" FILES
  ${UI_RESOURCES})

SOURCE_GROUP("Generated" FILES
  ${MOC_BUILT_SOURCES}
  ${UI_BUILT_SOURCES}
)

SET(ParaView_SOURCE_FILES
  myMainWindow.cxx
  myMainWindow.h
  module.h
  alyawidget.cpp
  alyawidget.h
<<sources>>
  ${MOC_BUILT_SOURCES}
  ${UI_BUILT_SOURCES})

INCLUDE_DIRECTORIES(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR})

build_paraview_client(paraviewAlya
  TITLE "ParaView ALYA"
  ORGANIZATION  "BSC-CNS"
  VERSION_MAJOR 1 
  VERSION_MINOR 1
  VERSION_PATCH 1
#  SPLASH_IMAGE "${CMAKE_CURRENT_SOURCE_DIR}/RSplash.png"
  PVMAIN_WINDOW myMainWindow
  PVMAIN_WINDOW_INCLUDE myMainWindow.h
  GUI_CONFIGURATION_XMLS
    ${CMAKE_CURRENT_SOURCE_DIR}/ParaViewSources.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/ParaViewFilters.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/ParaViewReaders.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/ParaViewWriters.xml
  SOURCES ${ParaView_SOURCE_FILES} ${UI_BUILT_SOURCES} ${RCS_BUILT_SOURCES}
)

TARGET_LINK_LIBRARIES(paraviewAlya ${QT_QTCORE_LIBRARY} ${QT_QTGUI_LIBRARY}  vtkCommon vtkPVCommon pqCore pqComponents)
