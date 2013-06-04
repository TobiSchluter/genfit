cmake_minimum_required(VERSION 2.8)

IF(NOT DEFINED ENV{ROOTSYS})
   MESSAGE(FATAL_ERROR "You did not set the environment variable ROOTSYS. Please check your ROOT configuration.") 
ENDIF(NOT DEFINED ENV{ROOTSYS})

IF(NOT DEFINED ENV{GENFIT})
   MESSAGE(FATAL_ERROR "You did not set the environment variable GENFIT.") 
ENDIF(NOT DEFINED ENV{GENFIT})

#IF(NOT DEFINED ENV{GLOGPATH})
#   MESSAGE(FATAL_ERROR "You did not set the environment variable GLOGPATH. Please install the Google Logging Library (glog) and set GLOGPATH.") 
#ENDIF(NOT DEFINED ENV{GLOGPATH})

FIND_PROGRAM(ROOT_CONFIG root-config $ENV{ROOTSYS}/bin)

IF(ROOT_CONFIG MATCHES "ROOT_CONFIG-NOTFOUND")
  MESSAGE(FATAL_ERROR "Could not find root-config executable.")
ENDIF(ROOT_CONFIG MATCHES "ROOT_CONFIG-NOTFOUND")

Execute_process(COMMAND ${ROOT_CONFIG} --cflags
	OUTPUT_VARIABLE ROOT_CFLAGS
    	OUTPUT_STRIP_TRAILING_WHITESPACE)

Execute_process(COMMAND ${ROOT_CONFIG} --glibs --ldflags --evelibs
	OUTPUT_VARIABLE ROOT_LIBS
    	OUTPUT_STRIP_TRAILING_WHITESPACE)

Execute_process(COMMAND ${ROOT_CONFIG} --incdir
	OUTPUT_VARIABLE ROOT_INCDIR
    	OUTPUT_STRIP_TRAILING_WHITESPACE)

#include(macros.cmake)

#SET(CMAKE_VERBOSE_MAKEFILE ON)

SET(CMAKE_CXX_FLAGS "-g -O2 -Wall -Wshadow -std=c++11 -pedantic -ftree-vectorize -ftree-vectorizer-verbose=0 ${ROOT_CFLAGS}")


include(CheckCXXCompilerFlag)
check_cxx_compiler_flag(-march=native HAS_MARCH_NATIVE)
IF(HAS_MARCH_NATIVE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
ENDIF(HAS_MARCH_NATIVE)

IF(NOT DEFINED ENV{GENFIT})
  SET(ENV{GENFIT} "${CMAKE_SOURCE_DIR}")
ENDIF(NOT DEFINED ENV{GENFIT})
get_filename_component(GENFIT $ENV{GENFIT} REALPATH)

SET(GENFITLIBDIR ${GENFIT}/lib)
SET(GENFITLIBINSTALLDIR ${GENFITLIBDIR})
SET(LIBRARY_OUTPUT_PATH ${GENFITLIBDIR})

SET(GENFITINCDIR ${GENFIT}/include/genfit2)
SET(GENFITINCINSTALLDIR ${GENFITINCDIR})
SET(INCLUDE_OUTPUT_PATH ${GENFITINCDIR})

INCLUDE_DIRECTORIES(
	${GENFIT}/core/include
	${GENFIT}/RKTrackRep/include
	${GENFIT}/MeasurementExamples/include
	${ROOT_INCDIR}
	$ENV{GLOGPATH}/include
)
add_subdirectory(core/src)
add_subdirectory(RKTrackRep/src)
add_subdirectory(core/include)
add_subdirectory(RKTrackRep/include)
#add_subdirectory(GFRave)
#add_subdirectory(MeasurementExamples/src)
#add_subdirectory(MeasurementExamples/include)

ADD_LIBRARY(genfit2 SHARED ${GENFIT_SRCS})
TARGET_LINK_LIBRARIES(genfit2 ${ROOT_LIBS})
INSTALL(TARGETS genfit2 DESTINATION ${GENFITLIBINSTALLDIR})

