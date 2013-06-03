cmake_minimum_required(VERSION 2.8)

#check for ENV variables that are needed
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

#SET(CMAKE_VERBOSE_MAKEFILE ON)


include(CheckCXXCompilerFlag)
check_cxx_compiler_flag(-march=native HAS_MARCH_NATIVE)
IF(HAS_MARCH_NATIVE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
ENDIF(HAS_MARCH_NATIVE)

INCLUDE_DIRECTORIES(
	${GENFIT}/core/include
	${GENFIT}/RKTrackRep/include
	${GENFIT}/MeasurementExamples/include
	${ROOT_INCDIR}
)


AUX_SOURCE_DIRECTORY(${GENFIT}/MeasurementExamples/src MEASUREMENTEXAMPLES_SRCS)

SET(MEASUREMENTEXAMPLES_HEADERS 
  PixMeasurement.h
  SpacepointMeasurement.h 
  ProlateSpacepointMeasurement.h
  StripMeasurement.h
  WireMeasurement.h
  WirePointMeasurement.h
)


################ install ###################

#add ROOT lib
Execute_process(COMMAND root-config --glibs --cflags
	OUTPUT_VARIABLE ROOT_LIBS
    	OUTPUT_STRIP_TRAILING_WHITESPACE)


ADD_LIBRARY(MeasurementExamples SHARED ${MEASUREMENTEXAMPLES_SRCS})
TARGET_LINK_LIBRARIES(MeasurementExamples ${ROOT_LIBS} -L${GENFITLIBDIR} -lgenfit2 -lgenfit2RK  -lGeom -lVMC -lEG)
INSTALL(TARGETS MeasurementExamples DESTINATION ${GENFITLIBINSTALLDIR} )

