# #require a certain version
CMAKE_MINIMUM_REQUIRED( VERSION 2.8 )

# declare the project name
PROJECT(genfit)

# set project version
SET( ${PROJECT_NAME}_VERSION_MAJOR 2 )
SET( ${PROJECT_NAME}_VERSION_MINOR 0 )
SET( ${PROJECT_NAME}_VERSION_PATCH 0 )

# define path to macros and settings
LIST( APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake )
# include the macros and settings
INCLUDE(settings)
INCLUDE(macros)

FIND_PACKAGE(ROOT REQUIRED)


#search for boost, provide some information to user
FIND_PACKAGE(Boost)
IF (Boost_FOUND)
    INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
    MESSAGE("-- Boost is found at: ${Boost_INCLUDE_DIR}")
ELSE()
    MESSAGE( " Boost was not found on your system, using specified path in CMakeLists.txt file.")
    SET(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} "/usr/include")
ENDIF()


# require proper c++
ADD_DEFINITIONS( "-O2 -Wall -g -ansi -pedantic -Wno-long-long -Wshadow -pedantic -ftree-vectorize -ftree-vectorizer-verbose=0 "   )

ADD_DEFINITIONS("-std=c++03")
ADD_DEFINITIONS(-D_GFOVERRIDE=)
ADD_DEFINITIONS(-D_GFNULLPTR=NULL)

# include directories -- for now both core and RKTrackRep
INCLUDE_DIRECTORIES( BEFORE ./core/include  ${ROOT_INCLUDE_DIRS})
INSTALL( DIRECTORY ./core/include/ ./RKTrackRep/include/  DESTINATION ./include/genfit PATTERN ".svn" EXCLUDE )
#INCLUDE_DIRECTORIES( BEFORE ./core/include ./RKTrackRep/include ${ROOT_INCLUDE_DIRS})
#INSTALL( DIRECTORY ./core/include ./RKTrackRep/include  DESTINATION . PATTERN ".svn" EXCLUDE )


## normally this would then stand inside an if clause
#INCLUDE_DIRECTORIES( BEFORE  ./RKTrackRep/include)
#INSTALL( DIRECTORY ./RKTrackRep/include  DESTINATION . PATTERN ".svn" EXCLUDE )
## and this as well
#INCLUDE_DIRECTORIES( BEFORE  ./eventDisplay/include)
#INSTALL( DIRECTORY ./eventDisplay/include  DESTINATION . PATTERN ".svn" EXCLUDE )


# declare the stupid ROOT library path
LINK_DIRECTORIES(${ROOT_LIBRARY_DIR})

# add library
AUX_SOURCE_DIRECTORY( ./core/src  library_sources )

## and again if-clauses:
AUX_SOURCE_DIRECTORY( ./RKTrackRep/src  library_sources )

ADD_SHARED_LIBRARY( ${PROJECT_NAME} ${library_sources} )
INSTALL( TARGETS ${PROJECT_NAME} DESTINATION lib )

TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${ROOT_LIBS})

# the examples, build with "make examples"
#ADD_SUBDIRECTORY( examples/ )

# generate and install following configuration files
#GENERATE_PACKAGE_CONFIGURATION_FILES( genfitConfig.cmake genfitConfigVersion.cmake )
