#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#ifdef __CLING__
#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;

#pragma link C++ namespace genfit;
#endif

// These need no special treatment.
#pragma link C++ class genfit::EventDisplay+;

#endif
