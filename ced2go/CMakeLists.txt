########################################################
CMAKE_MINIMUM_REQUIRED( VERSION 2.6 )
########################################################

# project name
#PROJECT( extractdetector )

#ENABLE_TESTING()
#INCLUDE(CTest)

# append link pathes to rpath list
#SET( CMAKE_INSTALL_RPATH_USE_LINK_PATH 1 )
#MARK_AS_ADVANCED( CMAKE_INSTALL_RPATH_USE_LINK_PATH )

FIND_PACKAGE( LCIO REQUIRED )
INCLUDE_DIRECTORIES( ${LCIO_INCLUDE_DIRS} )

ADD_EXECUTABLE( extractdetector extractdetector )
TARGET_LINK_LIBRARIES( extractdetector ${LCIO_LIBRARIES} )

