#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#ifdef __CLING__
#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;
#pragma link C++ namespace genfit;
#endif

#pragma link C++ class genfit::GFGbl+;
#pragma link C++ class genfit::GblFitter+;
#pragma link C++ class genfit::ICalibrationParametersDerivatives+;
#pragma link C++ class genfit::GblFitStatus+;
#pragma link C++ class genfit::GblFitterInfo+;
#pragma link C++ class genfit::GblTrackSegmentController+;

#endif
