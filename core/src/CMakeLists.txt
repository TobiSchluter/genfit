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
	${ROOT_INCDIR}
	$ENV{GLOGPATH}/include
)

set(GENFIT_SRCS 
	AbsFinitePlane.cpp
	AbsFitterInfo.cpp
	AbsMeasurement.cpp
	AbsTrackRep.cpp
	ConstField.cpp
	DetPlane.cpp
	Exception.cpp
	FieldManager.cpp
	KalmanFittedStateOnPlane.cpp
	KalmanFitterInfo.cpp
	KalmanFitter.cpp
	MaterialInfo.cpp
	MaterialProperties.cpp
	MeasuredStateOnPlane.cpp
	MeasurementOnPlane.cpp
	PlanarMeasurement.cpp
	PlanarPixelMeasurement.cpp
	PlanarStripMeasurement.cpp
	ProlateSpacePointMeasurement.cpp
	RectangularFinitePlane.cpp
	ReferenceStateOnPlane.cpp
	SimpleKalmanFitterInfo.cpp
	SpacePointMeasurement.cpp
	StateOnPlane.cpp
	Tools.cpp
	TrackCand.cpp
	TrackCandHit.cpp
	Track.cpp
	TrackPoint.cpp
	WireMeasurement.cpp
	WirePointMeasurement.cpp
)


################ install ###################

#add ROOT lib
Execute_process(COMMAND ${ROOT_CONFIG} --glibs --ldflags --evelibs
	OUTPUT_VARIABLE ROOT_LIBS
    	OUTPUT_STRIP_TRAILING_WHITESPACE)

ADD_LIBRARY(genfit2 SHARED ${GENFIT_SRCS})

TARGET_LINK_LIBRARIES(genfit2 ${ROOT_LIBS})

INSTALL(TARGETS genfit2 DESTINATION ${GENFITLIBINSTALLDIR})

