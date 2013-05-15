cmake_minimum_required(VERSION 2.8)

set(GENFIT_INCS
  AbsBField.h
  AbsFinitePlane.h
  AbsFitter.h
  AbsFitterInfo.h
  AbsMeasurement.h
  AbsTrackRep.h
  ConstField.h
  DetPlane.h
  Exception.h
  FieldManager.h
  genfitLinkDef.h
  KalmanFittedStateOnPlane.h
  KalmanFitter.h
  KalmanFitterInfo.h
  MaterialInfo.h
  MaterialProperties.h
  MeasuredStateOnPlane.h
  MeasurementOnPlane.h
  PlanarMeasurement.h
  PlanarPixelMeasurement.h
  PlanarStripMeasurement.h
  ProlateSpacePointMeasurement.h
  RectangularFinitePlane.h
  ReferenceStateOnPlane.h
  SharedMaterialPropertiesPtr.h
  SharedPlanePtr.h
  SpacePointMeasurement.h
  StateOnPlane.h
  Tools.h
  TrackCand.h
  TrackCandHit.h
  Track.h
  TrackPoint.h
  WireMeasurement.h
  WirePointMeasurement.h)

INSTALL(FILES ${GENFIT_INCS} DESTINATION ${GENFITINCINSTALLDIR})