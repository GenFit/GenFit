#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

// these need no special tratment
#pragma link C++ class genfit::AbsKalmanFitter+;
#pragma link C++ class genfit::KalmanFitStatus;
#pragma link C++ class genfit::KalmanFitterRefTrack+;

// these inherit from classes that need custom streamers
#pragma link C++ class genfit::KalmanFittedStateOnPlane+;
#pragma link C++ class genfit::ReferenceStateOnPlane+;

// Classes that needed manually written streamers:
#pragma link C++ class genfit::KalmanFitter-;
#pragma link C++ class genfit::KalmanFitterInfo-;
#pragma link C++ class genfit::DAF-;

#endif
