cmake_minimum_required(VERSION 2.8)

set(GENFIT_SRCS
        ${GENFIT_SRCS}
	${CMAKE_CURRENT_SOURCE_DIR}/PixMeasurement.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/SpacepointMeasurement.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/ProlateSpacepointMeasurement.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/StripMeasurement.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/WireMeasurement.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/WirePointMeasurement.cpp
 PARENT_SCOPE
)

