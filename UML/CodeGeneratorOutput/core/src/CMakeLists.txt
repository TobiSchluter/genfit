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
	KalmanFitterInfo.cpp        
	ProlateSpacePointMeasurement.cpp  
	TrackCandHit.cpp
	AbsFitterInfo.cpp             
	MaterialInfo.cpp            
	RectangularFinitePlane.cpp        
	Track.cpp
	AbsMeasurement.cpp            
	MeasuredStateOnPlane.cpp    	
	ReferenceStateOnPlane.cpp         
	TrackPoint.cpp
	AbsTrackRep.cpp               
	PlanarMeasurement.cpp       
	SpacePointMeasurement.cpp         
	WireMeasurement.cpp
	DetPlane.cpp                  
	PlanarPixelMeasurement.cpp  
	StateOnPlane.cpp                  
	WirePointMeasurement.cpp
	KalmanFittedStateOnPlane.cpp  
	PlanarStripMeasurement.cpp  
	TrackCand.cpp
)

set(GENFIT_HEADERS 
	../include/AbsFinitePlane.h            
	../include/KalmanFitterInfo.h        
	../include/ProlateSpacePointMeasurement.h  
	../include/TrackCandHit.h
	../include/AbsFitterInfo.h             
	../include/MaterialInfo.h            
	../include/RectangularFinitePlane.h        
	../include/Track.h
	../include/AbsMeasurement.h            
	../include/MeasuredStateOnPlane.h    	
	../include/ReferenceStateOnPlane.h         
	../include/TrackPoint.h
	../include/AbsTrackRep.h               
	../include/PlanarMeasurement.h       
	../include/SpacePointMeasurement.h         
	../include/WireMeasurement.h
	../include/DetPlane.h                  
	../include/PlanarPixelMeasurement.h  
	../include/StateOnPlane.h                  
	../include/WirePointMeasurement.h
	../include/KalmanFittedStateOnPlane.h  
	../include/PlanarStripMeasurement.h  
	../include/TrackCand.h
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
      	-I${GENFIT}/core ${GENFIT_HEADERS} ${GENFIT_LINKDEF}
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
