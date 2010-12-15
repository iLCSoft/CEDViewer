########################################################
# cmake file for building Marlin example Package
# @author Jan Engels, Desy IT
CMAKE_MINIMUM_REQUIRED(VERSION 2.6 FATAL_ERROR)
########################################################


# project name
PROJECT( CEDViewer )

# project version
SET( ${PROJECT_NAME}_VERSION_MAJOR 1 )
SET( ${PROJECT_NAME}_VERSION_MINOR 1 )
SET( ${PROJECT_NAME}_VERSION_PATCH 0 )



### SETTINGS ################################################################

# ilcsoft cmake macros / default settings
INCLUDE( ilcsoft_default_settings )
#INCLUDE( build_32bit_compatible )
INCLUDE( MacroAddSharedLibrary )



### DEPENDENCIES ############################################################

# additionally search for cmake macros/moules defined by this project
LIST( APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake )

# ------- required packages -------------------
SET( required_packages LCIO GEAR CED Marlin MarlinUtil streamlog CLHEP GSL )

# ------- optional packages -------------------
#SET( optional_packages AIDA )


# ------- do the search -----------------------

# fix for backwards compatibility
FOREACH( pkg ${required_packages} ${optional_packages} )
    IF( NOT ${pkg}_DIR AND DEFINED ${pkg}_HOME )
        SET( ${pkg}_DIR ${${pkg}_HOME} )
    ENDIF()
ENDFOREACH()

IF( NOT streamlog_DIR AND DEFINED Marlin_DIR )
    SET( streamlog_DIR "${Marlin_DIR}" CACHE PATH "Path to streamlog" FORCE )
ENDIF()

FOREACH( pkg ${required_packages} )
    FIND_PACKAGE( ${pkg} REQUIRED )
ENDFOREACH()

FOREACH( pkg ${optional_packages} )
    FIND_PACKAGE( ${pkg} QUIET )
ENDFOREACH()

FOREACH( pkg ${required_packages} ${optional_packages} )
    IF( ${pkg}_FOUND )
        INCLUDE_DIRECTORIES( ${${pkg}_INCLUDE_DIRS} )
        LINK_LIBRARIES( ${${pkg}_LIBRARIES} )
        ADD_DEFINITIONS ( ${${pkg}_DEFINITIONS} )
    ELSE()
        MESSAGE( STATUS "${pkg} not found" )
    ENDIF()
ENDFOREACH()



### DOCUMENTATION ###########################################################

OPTION( INSTALL_DOC "Set to OFF to skip build/install Documentation" ON )

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
INCLUDE_DIRECTORIES( "include" )
INSTALL( DIRECTORY include DESTINATION . FILES_MATCHING PATTERN "*.h" )

# definitions to pass to the compiler
ADD_DEFINITIONS( "-Wall -ansi -pedantic" )
ADD_DEFINITIONS( "-Wno-long-long" )

# add library
AUX_SOURCE_DIRECTORY( src library_sources )
ADD_SHARED_LIBRARY( ${libname} ${library_sources} )

# ced2go
ADD_EXECUTABLE( extractdetector ced2go/extractdetector.cc )
INSTALL( TARGETS extractdetector DESTINATION bin )
INSTALL( FILES ced2go/ced2go-template.xml DESTINATION bin )
INSTALL( FILES ced2go/ced2go DESTINATION bin PERMISSIONS
        OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ GROUP_EXECUTE
        WORLD_READ WORLD_EXECUTE
)

# display some variables and write them to cache
DISPLAY_STD_VARIABLES()

