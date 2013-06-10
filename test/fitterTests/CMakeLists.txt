cmake_minimum_required(VERSION 2.6)

IF(NOT DEFINED ENV{ROOTSYS})
   MESSAGE(FATAL_ERROR "You did not set the environment variable ROOTSYS. Please check your ROOT configuration.") 
ENDIF(NOT DEFINED ENV{ROOTSYS})

IF(NOT DEFINED ENV{GENFIT})
   MESSAGE(FATAL_ERROR "You did not set the environment variable GENFIT.") 
ENDIF(NOT DEFINED ENV{GENFIT})

FIND_PROGRAM(ROOT_CONFIG root-config $ENV{ROOTSYS}/bin)

IF(ROOT_CONFIG MATCHES "ROOT_CONFIG-NOTFOUND")
  MESSAGE(FATAL_ERROR "Could not find root-config executable.")
ENDIF(ROOT_CONFIG MATCHES "ROOT_CONFIG-NOTFOUND")

Execute_process(COMMAND ${ROOT_CONFIG} --cflags
	OUTPUT_VARIABLE ROOT_CFLAGS
    	OUTPUT_STRIP_TRAILING_WHITESPACE)
    	
Execute_process(COMMAND ${ROOT_CONFIG} --incdir
	OUTPUT_VARIABLE ROOT_INCDIR
    	OUTPUT_STRIP_TRAILING_WHITESPACE)
    	
SET(CMAKE_CXX_FLAGS "-g -O2 -Wall -Wshadow -std=c++11 -pedantic ${ROOT_CFLAGS}")
#SET(CMAKE_VERBOSE_MAKEFILE ON)

Execute_process(COMMAND ${ROOT_CONFIG} --glibs --ldflags --evelibs
	OUTPUT_VARIABLE ROOT_LIBS
    	OUTPUT_STRIP_TRAILING_WHITESPACE)

INCLUDE_DIRECTORIES(
	${ROOT_INCDIR}
	$ENV{GENFIT}/core/include
	$ENV{GENFIT}/RKTrackRep/include
	$ENV{GENFIT}/eventDisplay/include
)

ADD_EXECUTABLE(main 
	main.cpp
	$ENV{GENFIT}/eventDisplay/include/EventDisplay.h
	$ENV{GENFIT}/eventDisplay/src/EventDisplay.cpp
)

TARGET_LINK_LIBRARIES(main 
	${ROOT_LIBS} 
	-L$ENV{GENFIT}/lib
	-lgenfit2 
)

