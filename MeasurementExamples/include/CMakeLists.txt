cmake_minimum_required(VERSION 2.8)

set(GENFIT_INCS
  PixMeasurement.h
  SpacepointMeasurement.h
  WireMeasurement.h
  ProlateSpacepointMeasurement.h
  StripMeasurement.h
  WirePointMeasurement.h
)

INSTALL(FILES ${GENFIT_INCS} DESTINATION ${GENFITINCINSTALLDIR})
