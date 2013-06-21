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

# require proper c++
ADD_DEFINITIONS( "-O2 -Wall -g -ansi -pedantic -Wno-long-long -std=c++11 -Wshadow -pedantic -ftree-vectorize -ftree-vectorizer-verbose=0 "   )

# include directories -- for now both core and RKTrackRep
INCLUDE_DIRECTORIES( BEFORE ./core/include ./RKTrackRep/include ${ROOT_INCLUDE_DIRS})
INSTALL( DIRECTORY ./core/include/ ./RKTrackRep/include/  DESTINATION ./include/genfit PATTERN ".svn" EXCLUDE )


## normally this would then stand inside an if clause
#INCLUDE_DIRECTORIES( BEFORE  ./RKTrackRep/include)
#INSTALL( DIRECTORY ./RKTrackRep/include  DESTINATION . PATTERN ".svn" EXCLUDE )
## and this as well
#INCLUDE_DIRECTORIES( BEFORE  ./eventDisplay/include)
#INSTALL( DIRECTORY ./eventDisplay/include  DESTINATION . PATTERN ".svn" EXCLUDE )

# add a target to generate API documentation with Doxygen
FIND_PACKAGE(Doxygen)

if(DOXYGEN_FOUND)
    
    CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile @ONLY)

    ADD_CUSTOM_TARGET(doc
        ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/doc
        COMMENT "Generating API documentation with Doxygen" VERBATIM
    )

ELSE()
    MESSAGE( STATUS "Doxygen not installed or found -- API documentation can't be generatred!" )

ENDIF(DOXYGEN_FOUND)


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
