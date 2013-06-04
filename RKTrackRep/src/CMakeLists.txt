cmake_minimum_required(VERSION 2.8)

set(GENFIT_SRCS
        ${GENFIT_SRCS}
	${CMAKE_CURRENT_SOURCE_DIR}/MaterialEffects.cpp  
	${CMAKE_CURRENT_SOURCE_DIR}/RKTools.cpp  
	${CMAKE_CURRENT_SOURCE_DIR}/RKTrackRep.cpp  
	${CMAKE_CURRENT_SOURCE_DIR}/StepLimits.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/TGeoMaterialInterface.cpp
 PARENT_SCOPE
)
