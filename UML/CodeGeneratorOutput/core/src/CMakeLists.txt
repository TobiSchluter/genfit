cmake_minimum_required(VERSION 2.8)

#check for ENV variables that are needed
FIND_PROGRAM(ROOT_CONFIG root-config $ENV{ROOTSYS}/bin)

IF(ROOT_CONFIG MATCHES "ROOT_CONFIG-NOTFOUND")
  MESSAGE(FATAL_ERROR "Could not find root-config executable.")
ENDIF(ROOT_CONFIG MATCHES "ROOT_CONFIG-NOTFOUND")

Execute_process(COMMAND ${ROOT_CONFIG} --cflags
	OUTPUT_VARIABLE ROOT_CFLAGS
    	OUTPUT_STRIP_TRAILING_WHITESPACE)

#SET(CMAKE_VERBOSE_MAKEFILE ON)


include(CheckCXXCompilerFlag)
check_cxx_compiler_flag(-march=native HAS_MARCH_NATIVE)
IF(HAS_MARCH_NATIVE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
ENDIF(HAS_MARCH_NATIVE)

INCLUDE_DIRECTORIES(
	${GENFIT}/core/include
)

set(GENFIT_SRCS 
	AbsFinitePlane.cpp
	AbsFitterInfo.cpp
	AbsMeasurement.cpp
	AbsTrackRep.cpp
	DetPlane.cpp
	KalmanFittedStateOnPlane.cpp
	KalmanFitterInfo.cpp
	MaterialInfo.cpp
	MeasuredStateOnPlane.cpp
	MeasurementOnPlane.cpp
	PlanarMeasurement.cpp
	PlanarPixelMeasurement.cpp
	PlanarStripMeasurement.cpp
	ProlateSpacePointMeasurement.cpp
	RectangularFinitePlane.cpp
	ReferenceStateOnPlane.cpp
	SpacePointMeasurement.cpp
	StateOnPlane.cpp
	TrackCand.cpp
	TrackCandHit.cpp
	Track.cpp
	TrackPoint.cpp
	WireMeasurement.cpp
	WirePointMeasurement.cpp
)

set(GENFIT_HEADERS 
	${CMAKE_CURRENT_SOURCE_DIR}/../include/AbsFinitePlane.h
	${CMAKE_CURRENT_SOURCE_DIR}/../include/AbsFitterInfo.h
	${CMAKE_CURRENT_SOURCE_DIR}/../include/AbsMeasurement.h
	${CMAKE_CURRENT_SOURCE_DIR}/../include/AbsTrackRep.h
	${CMAKE_CURRENT_SOURCE_DIR}/../include/DetPlane.h
	${CMAKE_CURRENT_SOURCE_DIR}/../include/KalmanFittedStateOnPlane.h
	${CMAKE_CURRENT_SOURCE_DIR}/../include/KalmanFitterInfo.h
	${CMAKE_CURRENT_SOURCE_DIR}/../include/MaterialInfo.h
	${CMAKE_CURRENT_SOURCE_DIR}/../include/MeasuredStateOnPlane.h
	${CMAKE_CURRENT_SOURCE_DIR}/../include/MeasurementOnPlane.h
	${CMAKE_CURRENT_SOURCE_DIR}/../include/PlanarMeasurement.h
	${CMAKE_CURRENT_SOURCE_DIR}/../include/PlanarPixelMeasurement.h
	${CMAKE_CURRENT_SOURCE_DIR}/../include/PlanarStripMeasurement.h
	${CMAKE_CURRENT_SOURCE_DIR}/../include/ProlateSpacePointMeasurement.h
	${CMAKE_CURRENT_SOURCE_DIR}/../include/RectangularFinitePlane.h
	${CMAKE_CURRENT_SOURCE_DIR}/../include/ReferenceStateOnPlane.h
	${CMAKE_CURRENT_SOURCE_DIR}/../include/SpacePointMeasurement.h
	${CMAKE_CURRENT_SOURCE_DIR}/../include/StateOnPlane.h
	${CMAKE_CURRENT_SOURCE_DIR}/../include/TrackCand.h
	${CMAKE_CURRENT_SOURCE_DIR}/../include/TrackCandHit.h
	${CMAKE_CURRENT_SOURCE_DIR}/../include/Track.h
	${CMAKE_CURRENT_SOURCE_DIR}/../include/TrackPoint.h
	${CMAKE_CURRENT_SOURCE_DIR}/../include/WireMeasurement.h
	${CMAKE_CURRENT_SOURCE_DIR}/../include/WirePointMeasurement.h
)

#string (REPLACE ".cpp" ".h" GENFIT_HEADERS "${GENFIT_SRCS}")


set(GENFIT_LINKDEF ../include/genfitLinkDef.h)
set(GENFIT_DICTIONARY genfitDict.cpp) 

FIND_PROGRAM(ROOT_CINT_EXECUTABLE
	NAMES rootcint
	PATHS $ENV{ROOTSYS}/bin)

ADD_CUSTOM_COMMAND(OUTPUT ${CMAKE_BINARY_DIR}/${GENFIT_DICTIONARY}
	COMMAND ${ROOT_CINT_EXECUTABLE}
      	ARGS -f ${CMAKE_BINARY_DIR}/${GENFIT_DICTIONARY} -c -p -DHAVE_CONFIG_H
      	-I${GENFIT}/core/include ${GENFIT_HEADERS} ${GENFIT_LINKDEF}
      	DEPENDS ${GENFIT_HEADERS} ${GENFIT_LINKDEF})
      	
SET(GENFIT_SRCS ${GENFIT_SRCS} ${CMAKE_BINARY_DIR}/${GENFIT_DICTIONARY})

################ install ###################

#add ROOT lib
Execute_process(COMMAND root-config --glibs --cflags
	OUTPUT_VARIABLE ROOT_LIBS
    	OUTPUT_STRIP_TRAILING_WHITESPACE)


ADD_LIBRARY(genfit SHARED ${GENFIT_SRCS})

TARGET_LINK_LIBRARIES(genfit ${ROOT_LIBS})

INSTALL(TARGETS genfit DESTINATION ${GENFITLIBINSTALLDIR})
