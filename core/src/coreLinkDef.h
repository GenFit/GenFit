#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

// These need no special treatment.
#pragma link C++ class genfit::AbsFinitePlane+;
#pragma link C++ class genfit::AbsHMatrix+;
#pragma link C++ class genfit::RectangularFinitePlane+;
#pragma link C++ class genfit::FitStatus+;
#pragma link C++ class genfit::MaterialProperties+;
#pragma link C++ class genfit::TrackCand+;
#pragma link C++ class genfit::TrackCandHit+;

// These inherit from classes with custom streamers, or reference share_ptrs in their interfaces.
#pragma link C++ class genfit::AbsTrackRep;
#pragma link C++ class genfit::MeasuredStateOnPlane;

// These need their owners fixed up after reading.
#pragma link C++ class genfit::AbsMeasurement; // trackPoint_

// These cannot be dealt with by default streamers because of
// shared_ptrs<> or scoped_ptrs<>.  Additionally, they may need their
// owners fixed up.
#pragma link C++ class genfit::AbsFitterInfo-; // trackPoint_, rep_, sharedPlanePtr
#pragma link C++ class genfit::DetPlane-;  // scoped_ptr<> finitePlane_
#pragma link C++ class genfit::MeasurementOnPlane-; // scoped_ptr<> hMatrix_
#pragma link C++ class genfit::StateOnPlane-;  // rep_, sharedPlanePtr
#pragma link C++ class genfit::ThinScatterer-; // sharedPlanePtr
#pragma link C++ class genfit::Track-;
#pragma link C++ class genfit::TrackPoint-; // track_, fixup the map

#endif
