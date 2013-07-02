# #require a certain version
CMAKE_MINIMUM_REQUIRED( VERSION 2.6 )

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
    INSTALL( DIRECTORY ./core/include/ ./RKTrackRep/include/  DESTINATION ./include/genfit PATTERN ".svn" EXCLUDE )
ELSE()
    MESSAGE( " Boost was not found on your system, using the one shipped with GENFIT.")
    INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/extern/ )
    INSTALL( DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/extern/ ./core/include/ ./RKTrackRep/include/  DESTINATION ./include/genfit PATTERN ".svn" EXCLUDE )
ENDIF()

# require proper c++, standard is C++03, flag reads: c++98
ADD_DEFINITIONS( " -Wall -g -ansi -pedantic -Wno-long-long -Wshadow -pedantic -ftree-vectorize -ftree-vectorizer-verbose=0 "   )
ADD_DEFINITIONS("-std=c++03")
ADD_DEFINITIONS(-D_GFNULLPTR=NULL)

# include directories -- for now both core and RKTrackRep
INCLUDE_DIRECTORIES( BEFORE ./core/include ./RKTrackRep/include ${ROOT_INCLUDE_DIRS})
# and this is the rest:
#INCLUDE_DIRECTORIES( BEFORE  ./eventDisplay/include)
#INSTALL( DIRECTORY ./eventDisplay/include  DESTINATION . PATTERN ".svn" EXCLUDE )

# add a target to generate API documentation with Doxygen
FIND_PACKAGE(Doxygen)

if(DOXYGEN_FOUND)
    CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile @ONLY)

    ADD_CUSTOM_TARGET(doc
        ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/doc
        COMMENT "Generating API documentation with Doxygen" VERBATIM)
ELSE()
    MESSAGE( STATUS "Doxygen not installed or found -- API documentation can't be generatred!" )
ENDIF(DOXYGEN_FOUND)


# declare the stupid ROOT library path
LINK_DIRECTORIES(${ROOT_LIBRARY_DIR})

# add the list of source code files -- careful, when a new file is added!
AUX_SOURCE_DIRECTORY( ./core/src  library_sources )
AUX_SOURCE_DIRECTORY( ./RKTrackRep/src  library_sources )

ADD_SHARED_LIBRARY( ${PROJECT_NAME} ${library_sources} )
INSTALL( TARGETS ${PROJECT_NAME} DESTINATION lib )

# wtf I hate root
TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${ROOT_LIBS} -lGeom)

# the examples, build with "make tests"
#ADD_CUSTOM_TARGET( tests )
#MACRO( ADD_GENFIT_TEST _testname )
#   INCLUDE_DIRECTORIES( BEFORE eventDisplay/include)
#   ADD_EXECUTABLE( ${_testname} ${ARGN} )
#   ADD_DEPENDENCIES( tests  ${_testname} )
#   TARGET_LINK_LIBRARIES( ${_testname} ${PROJECT_NAME}  ${ROOT_LIBS} -lGeom )
    #INSTALL( TARGETS ${_testname} DESTINATION ${EXECUTABLE_INSTALL_DIR})
#ENDMACRO( ADD_GENFIT_TEST )

#ADD_GENFIT_TEST( fitterTests test/fitterTests/main.cpp eventDisplay/src/EventDisplay.cpp)


# generate and install following configuration files
#GENERATE_PACKAGE_CONFIGURATION_FILES( genfitConfig.cmake genfitConfigVersion.cmake )
