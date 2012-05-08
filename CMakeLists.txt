########################################################
# cmake file for building Marlin example Package
# @author Jan Engels, Desy IT
CMAKE_MINIMUM_REQUIRED(VERSION 2.6 FATAL_ERROR)
########################################################


# project name
PROJECT( CEDViewer )

# project version
SET( ${PROJECT_NAME}_VERSION_MAJOR 1 )
SET( ${PROJECT_NAME}_VERSION_MINOR 5 )
SET( ${PROJECT_NAME}_VERSION_PATCH 1 )

cmake_policy(SET CMP0008 NEW)

### SETTINGS ################################################################

FIND_PACKAGE( ILCUTIL COMPONENTS ILCSOFT_CMAKE_MODULES REQUIRED )

# load default settings from ILCSOFT_CMAKE_MODULES
INCLUDE( ilcsoft_default_settings )



### DEPENDENCIES ############################################################

FIND_PACKAGE( Marlin 1.0 REQUIRED )
FIND_PACKAGE( MarlinUtil 1.2 REQUIRED )

FIND_PACKAGE( LCCD ) # optional - needed for TrackerRawViewer
FIND_PACKAGE( MarlinTPC ) # optional - needed for TrackerRawViewer
#--- workaround for MarlinTPC (no cmake config files yet...)
#IF( MarlinTPC_DIR )
#  SET( MarlinTPC_INCLUDE_DIRS ${MarlinTPC_DIR}/tpcconddata/include )
#  IF( APPLE ) 
#    SET( MarlinTPC_LIBRARIES    ${MarlinTPC_DIR}/lib/libtpcconddata.dylib )
#  ELSE()
#    SET( MarlinTPC_LIBRARIES    ${MarlinTPC_DIR}/lib/libtpcconddata.so )
#  ENDIF()
#  SET( MarlinTPC_FOUND TRUE)
#ENDIF()


FOREACH( pkg Marlin MarlinUtil LCCD MarlinTPC)
    IF( ${pkg}_FOUND )
        INCLUDE_DIRECTORIES( ${${pkg}_INCLUDE_DIRS} )
        LINK_LIBRARIES( ${${pkg}_LIBRARIES} )
        ADD_DEFINITIONS ( ${${pkg}_DEFINITIONS} )
    ENDIF()
ENDFOREACH()







### DOCUMENTATION ###########################################################

OPTION( INSTALL_DOC "Set to OFF to skip build/install Documentation" OFF )

IF( INSTALL_DOC AND EXISTS "${PROJECT_SOURCE_DIR}/doc/CMakeLists.txt" )
    ADD_SUBDIRECTORY( doc )
ELSE()
    MESSAGE( STATUS "INSTALL_DOC set to OFF" )
    SET( INSTALL_DOC OFF )
ENDIF()



### LIBRARY AND TOOLS #######################################################

# library name
SET( libname ${PROJECT_NAME} )

# include directories
INCLUDE_DIRECTORIES( ./include )

# definitions to pass to the compiler
ADD_DEFINITIONS( "-Wall -ansi -pedantic" )
ADD_DEFINITIONS( "-Wno-long-long" )

# add library
AUX_SOURCE_DIRECTORY( ./src library_sources )

LIST( REMOVE_ITEM library_sources ./src/TrackerRawViewer.cc )


#don't build TrackerRawViewer if LCCD not available
IF( LCCD_FOUND AND MarlinTPC_FOUND )
  LIST( APPEND library_sources ./src/TrackerRawViewer.cc )
ELSE() 
  MESSAGE( "\n ***** TrackerRawViewer will not be built - it needs MarlinTPC and LCCD  \n"
           " ***** to include TrackerRawData specify -DLCCD_DIR=... and -DMarlinTPC_DIR=...   !! \n"  ) 
ENDIF()




ADD_SHARED_LIBRARY( ${libname} ${library_sources} )
INSTALL_SHARED_LIBRARY( ${libname} DESTINATION lib )

# ced2go
ADD_EXECUTABLE( extractdetector ced2go/extractdetector.cc )
INSTALL( TARGETS extractdetector DESTINATION bin )
INSTALL( FILES ./ced2go/ced2go-template.xml DESTINATION bin )
INSTALL( PROGRAMS ./ced2go/ced2go DESTINATION bin )


# display some variables and write them to cache
DISPLAY_STD_VARIABLES()

