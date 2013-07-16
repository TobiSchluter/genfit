#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

// These need no special treatment.
#pragma link C++ class genfit::AbsFinitePlane;
#pragma link C++ class genfit::RectangularFinitePlane;
#pragma link C++ class genfit::AbsTrackRep;
#pragma link C++ class genfit::MeasurementOnPlane;
#pragma link C++ class genfit::MeasuredStateOnPlane;
#pragma link C++ class genfit::FitStatus;

// These need their owners fixed up after reading.
#pragma link C++ class genfit::AbsMeasurement; // trackPoint_
#pragma link C++ class genfit::AbsFitterInfo; // trackPoint_, rep_

// These cannot be dealt with by default streamers because of
// shared_ptrs<> or scoped_ptrs<>.  Additionally, they may need their
// owners fixed up.
#pragma link C++ class genfit::DetPlane-;  // scoped_ptr<> finitePlane_
#pragma link C++ class genfit::StateOnPlane-;  // rep_, sharedPlanePtr
#pragma link C++ class genfit::Track-;
#pragma link C++ class genfit::TrackPoint-; // track_, fixup the map

#endif
