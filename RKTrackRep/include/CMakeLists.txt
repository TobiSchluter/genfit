cmake_minimum_required(VERSION 2.8)

set(GENFIT_INCS RKTrackRep.h AbsMaterialInterface.h StepLimits.h TGeoMaterialInterface.h MaterialEffects.h)

INSTALL(FILES ${GENFIT_INCS} DESTINATION ${GENFITINCINSTALLDIR})