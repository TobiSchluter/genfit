cmake_minimum_required(VERSION 2.8)

set(GENFIT_INCS
  ${CMAKE_CURRENT_SOURCE_DIR}/AbsBField.h
  ${CMAKE_CURRENT_SOURCE_DIR}/AbsFinitePlane.h
  ${CMAKE_CURRENT_SOURCE_DIR}/AbsFitter.h
  ${CMAKE_CURRENT_SOURCE_DIR}/AbsFitterInfo.h
  ${CMAKE_CURRENT_SOURCE_DIR}/AbsMeasurement.h
  ${CMAKE_CURRENT_SOURCE_DIR}/AbsTrackRep.h
  ${CMAKE_CURRENT_SOURCE_DIR}/ConstField.h
  ${CMAKE_CURRENT_SOURCE_DIR}/DetPlane.h
  ${CMAKE_CURRENT_SOURCE_DIR}/Exception.h
  ${CMAKE_CURRENT_SOURCE_DIR}/FieldManager.h
  ${CMAKE_CURRENT_SOURCE_DIR}/genfitLinkDef.h
  ${CMAKE_CURRENT_SOURCE_DIR}/KalmanFittedStateOnPlane.h
  ${CMAKE_CURRENT_SOURCE_DIR}/KalmanFitter.h
  ${CMAKE_CURRENT_SOURCE_DIR}/KalmanFitterRefTrack.h
  ${CMAKE_CURRENT_SOURCE_DIR}/KalmanFitterInfo.h
  ${CMAKE_CURRENT_SOURCE_DIR}/MaterialInfo.h
  ${CMAKE_CURRENT_SOURCE_DIR}/MaterialProperties.h
  ${CMAKE_CURRENT_SOURCE_DIR}/MeasuredStateOnPlane.h
  ${CMAKE_CURRENT_SOURCE_DIR}/MeasurementOnPlane.h
  ${CMAKE_CURRENT_SOURCE_DIR}/PlanarMeasurement.h
  ${CMAKE_CURRENT_SOURCE_DIR}/ProlateSpacepointMeasurement.h
  ${CMAKE_CURRENT_SOURCE_DIR}/RectangularFinitePlane.h
  ${CMAKE_CURRENT_SOURCE_DIR}/ReferenceStateOnPlane.h
  ${CMAKE_CURRENT_SOURCE_DIR}/SharedMaterialPropertiesPtr.h
  ${CMAKE_CURRENT_SOURCE_DIR}/SharedPlanePtr.h
  ${CMAKE_CURRENT_SOURCE_DIR}/SimpleKalmanFitterInfo.h
  ${CMAKE_CURRENT_SOURCE_DIR}/SpacepointMeasurement.h
  ${CMAKE_CURRENT_SOURCE_DIR}/StateOnPlane.h
  ${CMAKE_CURRENT_SOURCE_DIR}/Tools.h
  ${CMAKE_CURRENT_SOURCE_DIR}/TrackCand.h
  ${CMAKE_CURRENT_SOURCE_DIR}/TrackCandHit.h
  ${CMAKE_CURRENT_SOURCE_DIR}/Track.h
  ${CMAKE_CURRENT_SOURCE_DIR}/TrackPoint.h
  ${CMAKE_CURRENT_SOURCE_DIR}/WireMeasurement.h
  ${CMAKE_CURRENT_SOURCE_DIR}/WirePointMeasurement.h)

INSTALL(FILES ${GENFIT_INCS} DESTINATION ${GENFITINCINSTALLDIR})
