#include <iostream>
#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>

#include <AbsFinitePlane.h>
#include <AbsFitterInfo.h>
#include <AbsMeasurement.h>
#include <AbsTrackRep.h>
#include <ConstField.h>
#include <DetPlane.h>
#include <Exception.h>
#include <FieldManager.h>
#include <KalmanFittedStateOnPlane.h>
#include <KalmanFitterInfo.h>
#include <MeasuredStateOnPlane.h>
#include <MeasurementOnPlane.h>
#include <PlanarMeasurement.h>
#include <ProlateSpacepointMeasurement.h>
#include <RectangularFinitePlane.h>
#include <ReferenceStateOnPlane.h>
#include <SharedPlanePtr.h>
#include <SpacepointMeasurement.h>
#include <StateOnPlane.h>
#include <Tools.h>
#include <TrackCand.h>
#include <TrackCandHit.h>
#include <Track.h>
#include <TrackPoint.h>
#include <WireMeasurement.h>
#include <WirePointMeasurement.h>

#include <MaterialEffects.h>
#include <RKTools.h>
#include <RKTrackRep.h>
#include <StepLimits.h>
#include <TGeoMaterialInterface.h>

#include <TApplication.h>
#include <TCanvas.h>
#include <TDatabasePDG.h>
#include <TEveManager.h>
#include <TGeoManager.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TRandom.h>
#include <TStyle.h>
#include <TVector3.h>

#include <vector>
#include <map>
#include <stdio.h>
#include <stdlib.h>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include "TDatabasePDG.h"
#include <TMath.h>
#include <TString.h>

//#include <callgrind/callgrind.h>


enum e_testStatus {
  kPassed,
  kFailed,
  kException
};

constexpr bool verbose = false;

void handler(int sig) {
  void *array[10];
  size_t size;

  // get void*'s for all entries on the stack
  size = backtrace(array, 10);

  // print out all the frames to stderr
  fprintf(stderr, "Error: signal %d:\n", sig);
  backtrace_symbols_fd(array, size, 2);
  exit(1);
}

int randomPdg() {
  int pdg;

  switch(int(gRandom->Uniform(8))) {
  case 1:
    pdg = -11; break;
  case 2:
    pdg = 11; break;
  case 3:
    pdg = 13; break;
  case 4:
    pdg = -13; break;
  case 5:
    pdg = 211; break;
  case 6:
    pdg = -211; break;
  case 7:
    pdg = 2212; break;
  default:
    pdg = 211;
  }

  return pdg;
}


int randomSign() {
  if (gRandom->Uniform(1) > 0.5)
    return 1;
  return -1;
}


e_testStatus compareMatrices(const TMatrixTBase<double>& A, const TMatrixTBase<double>& B, double maxRelErr) {
  e_testStatus retVal = kPassed;

  // search max abs value
  double max(0);
  for (int i=0; i<A.GetNrows(); ++i) {
    for (int j=0; j<A.GetNcols(); ++j) {
      if (fabs(A(i,j)) > max)
        max = fabs(A(i,j));
    }
  }

  double maxAbsErr = maxRelErr*max;

  for (int i=0; i<A.GetNrows(); ++i) {
    for (int j=0; j<A.GetNcols(); ++j) {
      double absErr = A(i,j) - B(i,j);
      if ( fabs(absErr) > maxAbsErr ) {
        double relErr = A(i,j)/B(i,j) - 1;
        if ( fabs(relErr) > maxRelErr ) {
          if (verbose)  {
            std::cout << "compareMatrices: A("<<i<<","<<j<<") = " << A(i,j) << "  B("<<i<<","<<j<<") = " << B(i,j) << "     absErr = " << absErr << "    relErr = " << relErr << "\n";
          }
          retVal = kFailed;
        }
      }
    }
  }
  return retVal;
}

e_testStatus isCovMatrix(TMatrixTBase<double>& cov) {

  if (!(cov.IsSymmetric())) {
    if (verbose) {
      std::cout << "isCovMatrix: not symmetric\n";
    }
    return kFailed;
  }

  for (int i=0; i<cov.GetNrows(); ++i) {
    for (int j=0; j<cov.GetNcols(); ++j) {
       if (std::isnan(cov(i,j))) {
         if (verbose) {
           std::cout << "isCovMatrix: element isnan\n";
         }
         return kFailed;
       }
       if (i==j && cov(i,j) < 0) {
         if (verbose) {
           std::cout << "isCovMatrix: negative diagonal element\n";
         }
         return kFailed;
       }
    }
  }

  return kPassed;
}



e_testStatus checkSetGetPosMom(bool writeHisto = false) {

  if (writeHisto)
    return kPassed;

  double epsilonLen = 1.E-10;
  double epsilonMom = 1.E-10;

  int pdg = randomPdg();
  genfit::AbsTrackRep* rep;
  rep = new genfit::RKTrackRep(pdg);

  //TVector3 pos(0,0,0);
  TVector3 pos(gRandom->Gaus(0,0.1),gRandom->Gaus(0,0.1),gRandom->Gaus(0,0.1));
  TVector3 mom(0,0.5,gRandom->Gaus(0,0.3));
  mom.SetMag(0.5);
  mom *= randomSign();


  genfit::StateOnPlane state(rep);
  rep->setPosMom(state, pos, mom);

  // check if we can set another position in the same plane
  if (randomSign() == 1) {
    genfit::SharedPlanePtr plane = state.getPlane();
    const TVector3& u = plane->getU();
    const TVector3& v = plane->getV();

    // random position on plane
    pos += gRandom->Gaus() * u;
    pos += gRandom->Gaus() * v;

    // new random momentum
    mom.SetXYZ(0,0.5,gRandom->Gaus(0,0.3));
    mom.SetMag(0.5);
    mom *= randomSign();

    rep->setPosMom(state, pos, mom);

    // check if plane has changed
    if (state.getPlane() != plane) {
        if (verbose) {
          std::cout << "plane has changed unexpectedly! \n";
        }
      delete rep;
      return kFailed;
    }
  }


  // compare
  if ((pos - rep->getPos(state)).Mag() > epsilonLen ||
      (mom - rep->getMom(state)).Mag() > epsilonMom) {

    if (verbose) {
      state.Print();
      std::cout << "pos difference = " << (pos - rep->getPos(state)).Mag() << "\n";
      std::cout << "mom difference = " << (mom - rep->getMom(state)).Mag() << "\n";

      std::cout << std::endl;
    }
    delete rep;
    return kFailed;
  }

  delete rep;
  return kPassed;

}


e_testStatus compareForthBackExtrapolation(bool writeHisto = false) {
  static std::map<int, std::vector<TH2D*> > histoMap;

  static const int nPDGs(5);
  int pdgs[nPDGs]={0, 211,13,11,2212};
  static bool fill(true);
  if (fill) {
    for (int i = 0; i<nPDGs; ++i) {
      int pdg = pdgs[i];

      std::string Result;//string which will contain the result
      std::stringstream convert; // stringstream used for the conversion
      convert << pdg;//add the value of Number to the characters in the stream
      Result = convert.str();//set Result to the content of the stream

      histoMap[pdg].push_back(new TH2D((std::string("deviationRel_")+Result).c_str(), "log(betaGamma) vs relative deviation", 100000, -1.e-2, 1.e-2, 50, -4, 8));
      histoMap[pdg].push_back(new TH2D((std::string("deviationAbs_")+Result).c_str(), "log(betaGamma) vs absolute deviation; deviation (keV)", 100000, -90.0, 10.0, 50, -4, 8));
      histoMap[pdg].push_back(new TH2D((std::string("ExtrapLen_")+Result).c_str(), "delta ExtrapLen vs relative deviation", 50000, -5.e-2, 5.e-2, 400, -0.1, 0.1));
    }
    fill = false;
  }


  if (writeHisto) {
    TFile outfile("deviation.root", "recreate");
    outfile.cd();

    for (int i = 0; i<nPDGs; ++i) {
      int pdg = pdgs[i];
      histoMap[pdg][0]->Write();
      histoMap[pdg][1]->Write();
      histoMap[pdg][2]->Write();
    }

    outfile.Close();

    return kPassed;
  }


  double epsilonLen = 5.E-5; // 0.5 mu
  double epsilonMom = 1.E-6; // 1 keV

  int pdg = randomPdg();
  //pdg = 211;
  genfit::AbsTrackRep* rep;
  rep = new genfit::RKTrackRep(pdg);

  //rep->setDebugLvl(1);
  //genfit::MaterialEffects::getInstance()->setDebugLvl(1);



  TParticlePDG* part = TDatabasePDG::Instance()->GetParticle(pdg);
  double mass = part->Mass(); // GeV

  //TVector3 pos(0,0,0);
  TVector3 pos(gRandom->Gaus(0,0.1),gRandom->Gaus(0,0.1),gRandom->Gaus(0,0.1));
  TVector3 mom(0,0.5,gRandom->Gaus(0,0.3));
  mom.SetMag( exp(gRandom->Uniform(-4, 8)) * mass );
  mom *= randomSign();


  double betaGamma = log(mom.Mag()/mass);

  //mom.Print();

  genfit::StateOnPlane state(rep);
  rep->setPosMom(state, pos, mom);

  genfit::SharedPlanePtr origPlane = state.getPlane();
  genfit::SharedPlanePtr plane(new genfit::DetPlane(TVector3(0,randomSign()*10,0), TVector3(0,randomSign()*1,0)));

  genfit::StateOnPlane origState(state);

  TVector3 mom2;
  double momLoss1, momLoss2;

  // forth
  double extrapLen(0);
  try {
    extrapLen = rep->extrapolateToPlane(state, plane);

    mom2 = state.getMom();
    momLoss1 = mom.Mag()-mom2.Mag();

    //mom2.Print();
    //std::cout << "MomLoss = " << momLoss1 << "\n";
  }
  catch (genfit::Exception& e) {
    if (verbose) {
      std::cerr << "Exception in forth Extrapolation. PDG = " << pdg << "; mom: \n";
      mom.Print();

      std::cerr << e.what();
    }
    delete rep;
    return kException;
  }



  // back
  double backExtrapLen(0);
  try {
    backExtrapLen = rep->extrapolateToPlane(state, origPlane);

    momLoss2 = mom2.Mag()-state.getMom().Mag();

    //state.getMom().Print();
    //std::cout << "MomLoss = " << momLoss2 << "\n";

    double deviation = 1. + momLoss1/momLoss2;
    histoMap[abs(pdg)][0]->Fill(deviation, betaGamma);
    histoMap[abs(pdg)][1]->Fill((mom.Mag() - state.getMom().Mag())*1e6, betaGamma);
    histoMap[abs(pdg)][2]->Fill(deviation, extrapLen+backExtrapLen);

    histoMap[0][0]->Fill(deviation, betaGamma);
    histoMap[0][1]->Fill((mom.Mag() - state.getMom().Mag())*1e6, betaGamma);
    histoMap[0][2]->Fill(deviation, extrapLen+backExtrapLen);

    //std::cout << "deviation = " << deviation << "\n";
  }
  catch (genfit::Exception& e) {
    if (verbose) {
      std::cerr << "Exception in back Extrapolation. PDG = " << pdg << "; mom:  \n";
      mom.Print();
      std::cerr << "mom2:  \n";
      mom2.Print();
    }
    std::cerr << e.what();

    delete rep;
    return kException;
  }

  // compare
  if ((rep->getPos(origState) - rep->getPos(state)).Mag() > epsilonLen ||
      (rep->getMom(origState) - rep->getMom(state)).Mag() > epsilonMom ||
      fabs(extrapLen + backExtrapLen) > epsilonLen) {

    if (verbose) {
        origState.Print();
        state.Print();

        std::cout << "pos difference = " << (rep->getPos(origState) - rep->getPos(state)).Mag() << "\n";
        std::cout << "mom difference = " << (rep->getMom(origState) - rep->getMom(state)).Mag() << "\n";
        std::cout << "len difference = " << extrapLen + backExtrapLen << "\n";

        std::cout << std::endl;
    }
    delete rep;
    return kFailed;
  }

  delete rep;
  return kPassed;

}


e_testStatus compareForthBackJacNoise(bool writeHisto = false) {

  if (writeHisto)
    return kPassed;

  e_testStatus retVal(kPassed);

  bool fx( randomSign() > 0);
  genfit::MaterialEffects::getInstance()->setNoEffects(!fx);

  double deltaJac = 0.005; // relative
  double deltaNoise = 0.00001;
  double deltaState = 3.E-6; // absolute

  if (fx) {
    deltaJac = 0.1; // relative
    deltaNoise = 0.1;
    deltaState = 5.E-4; // absolute
  }

  int pdg = randomPdg();
  genfit::AbsTrackRep* rep;
  rep = new genfit::RKTrackRep(pdg);

  //TVector3 pos(0,0,0);
  //TVector3 mom(0,1,2);
  TVector3 pos(gRandom->Gaus(0,0.1),gRandom->Gaus(0,0.1),gRandom->Gaus(0,0.1));
  TVector3 mom(0, 0.5, gRandom->Gaus(0, 1));
  mom *= randomSign();
  mom.SetMag(gRandom->Uniform(2)+0.3);
  //mom.SetMag(3);

  TMatrixD jac_f, jac_fi, jac_b, jac_bi;
  TMatrixDSym noise_f, noise_fi, noise_b, noise_bi;
  TVectorD c_f, c_fi, c_b, c_bi;
  TVectorD state_man, stateOrig_man;

  // original state and plane
  genfit::MeasuredStateOnPlane state(rep);
  rep->setPosMom(state, pos, mom);

  static const double smear = 0.2;
  TVector3 normal(mom);
  normal.SetMag(1);
  normal.SetXYZ(gRandom->Gaus(normal.X(), smear),
      gRandom->Gaus(normal.Y(), smear),
      gRandom->Gaus(normal.Z(), smear));
  genfit::DetPlane* origPlanePtr = new genfit::DetPlane (pos, normal);
  //genfit::DetPlane* origPlanePtr = new genfit::DetPlane (pos, TVector3(1,0,0), TVector3(0,0,1));
  double rotAngleOrig = gRandom->Uniform(2.*TMath::Pi());
  origPlanePtr->rotate(rotAngleOrig);
  genfit::SharedPlanePtr origPlane(origPlanePtr);
  //genfit::SharedPlanePtr origPlane = state.getPlane();
  rep->extrapolateToPlane(state, origPlane);

  const genfit::StateOnPlane origState(state);


  // dest plane
  normal = mom;
  normal.SetMag(1);
  normal.SetXYZ(gRandom->Gaus(normal.X(), smear),
      gRandom->Gaus(normal.Y(), smear),
      gRandom->Gaus(normal.Z(), smear));
  TVector3 dest(mom);
  dest.SetMag(10);
  genfit::DetPlane* planePtr = new genfit::DetPlane (dest, normal);
  //genfit::DetPlane* planePtr = new genfit::DetPlane (dest, TVector3(1,0,0), TVector3(0,0,1));
  double rotAngle = gRandom->Uniform(2.*TMath::Pi());
  planePtr->rotate(rotAngle);
  genfit::SharedPlanePtr plane(planePtr);

 /* genfit::DetPlane* planePtr = new genfit::DetPlane (*origPlane);
  planePtr->setO(TVector3(0,randomSign()*10,0));
  //planePtr->rotate(rotAngle);
  genfit::SharedPlanePtr plane(planePtr);
*/

  // numerical calculation
  TMatrixD jac_f_num;
  rep->calcJacobianNumerically(origState, plane, jac_f_num);


  // forth
  genfit::StateOnPlane extrapolatedState;
  try {
    //std::cout << "DO FORTH EXTRAPOLATION \n";
    rep->extrapolateToPlane(state, plane);
    //std::cout << "GET INFO FOR FORTH EXTRAPOLATION \n";
    extrapolatedState = state;
    rep->getForwardJacobianAndNoise(jac_f, noise_f, c_f);
    rep->getBackwardJacobianAndNoise(jac_fi, noise_fi, c_fi);
  }
  catch (genfit::Exception& e) {
    std::cerr << e.what();

    delete rep;
    genfit::MaterialEffects::getInstance()->setNoEffects(false);
    return kException;
  }

  // back
  try {
    //std::cout << "DO BACK EXTRAPOLATION \n";
    rep->extrapolateToPlane(state, origPlane);
    //std::cout << "GET INFO FOR BACK EXTRAPOLATION \n";
    rep->getForwardJacobianAndNoise(jac_b, noise_b, c_b);
    rep->getBackwardJacobianAndNoise(jac_bi, noise_bi, c_bi);
  }
  catch (genfit::Exception& e) {
    std::cerr << e.what();

    delete rep;
    genfit::MaterialEffects::getInstance()->setNoEffects(false);
    return kException;
  }


  // compare
  if (!((origState.getState() - state.getState()).Abs()  < deltaState) ) {
    if (verbose) {
      std::cout << "(origState.getState() - state.getState()) ";
      (origState.getState() - state.getState()).Print();
    }

    retVal = kFailed;
  }

  // check c
  if (!((jac_f * origState.getState() + c_f  -  extrapolatedState.getState()).Abs() < deltaState) ||
      !((jac_bi * origState.getState() + c_bi  -  extrapolatedState.getState()).Abs() < deltaState) ||
      !((jac_b * extrapolatedState.getState() + c_b  -  origState.getState()).Abs() < deltaState) ||
      !((jac_fi * extrapolatedState.getState() + c_fi  -  origState.getState()).Abs() < deltaState)   ) {

    if (verbose) {
      std::cout << "(jac_f * origState.getState() + c_f  -  extrapolatedState.getState()) ";
      (jac_f * origState.getState() + c_f - extrapolatedState.getState()).Print();
      std::cout << "(jac_bi * origState.getState() + c_bi  -  extrapolatedState.getState()) ";
      (jac_bi * origState.getState() + c_bi - extrapolatedState.getState()).Print();
      std::cout << "(jac_b * extrapolatedState.getState() + c_b  -  origState.getState()) ";
      (jac_b * extrapolatedState.getState() + c_b - origState.getState()).Print();
      std::cout << "(jac_fi * extrapolatedState.getState() + c_fi  -  origState.getState()) ";
      (jac_fi * extrapolatedState.getState() + c_fi - origState.getState()).Print();
    }
    retVal = kFailed;
  }

  if (isCovMatrix(state.getCov()) == kFailed) {
    retVal = kFailed;
  }


  // compare
  if (compareMatrices(jac_f, jac_bi, deltaJac) == kFailed) {
    if (verbose) {
      std::cout << "jac_f = ";
      jac_f.Print();
      std::cout << "jac_bi = ";
      jac_bi.Print();
      std::cout << std::endl;
    }
    retVal = kFailed;
  }

  // compare
  if (compareMatrices(jac_b, jac_fi, deltaJac) == kFailed) {
    if (verbose) {
      std::cout << "jac_b = ";
      jac_b.Print();
      std::cout << "jac_fi = ";
      jac_fi.Print();
      std::cout << std::endl;
    }
    retVal = kFailed;
  }

  // compare
  if (compareMatrices(noise_f, noise_bi, deltaNoise) == kFailed) {
    if (verbose) {
      std::cout << "noise_f = ";
      noise_f.Print();
      std::cout << "noise_bi = ";
      noise_bi.Print();
      std::cout << std::endl;
    }
    retVal = kFailed;
  }

  // compare
  if (compareMatrices(noise_b, noise_fi, deltaNoise) == kFailed) {
    if (verbose) {
      std::cout << "noise_b = ";
      noise_b.Print();
      std::cout << "noise_fi = ";
      noise_fi.Print();
      std::cout << std::endl;
    }
    retVal = kFailed;
  }


  if (!fx) {
    // compare
    if (compareMatrices(jac_f, jac_f_num, deltaJac) == kFailed) {
      if (verbose) {
        std::cout << "jac_f = ";
        jac_f.Print();
        std::cout << "jac_f_num = ";
        jac_f_num.Print();
        std::cout << std::endl;
      }
      retVal = kFailed;
    }
  }

  delete rep;
  genfit::MaterialEffects::getInstance()->setNoEffects(false);

  return retVal;
}


e_testStatus checkStopAtBoundary(bool writeHisto = false) {

  if (writeHisto)
    return kPassed;

  double epsilonLen = 1.E-4; // 1 mu
  //double epsilonMom = 1.E-5; // 10 keV

  double matRadius(1.);

  int pdg = randomPdg();
  genfit::AbsTrackRep* rep;
  rep = new genfit::RKTrackRep(pdg);

  //TVector3 pos(0,0,0);
  TVector3 pos(gRandom->Gaus(0,0.1),gRandom->Gaus(0,0.1),gRandom->Gaus(0,0.1));
  TVector3 mom(0,0.5,gRandom->Gaus(0,0.1));
  mom *= randomSign();


  genfit::StateOnPlane state(rep);
  rep->setPosMom(state, pos, mom);

  genfit::SharedPlanePtr origPlane = state.getPlane();
  genfit::SharedPlanePtr plane(new genfit::DetPlane(TVector3(0,randomSign()*10,0), TVector3(0,randomSign()*1,0)));

  genfit::StateOnPlane origState(state);

  // forth
  try {
    rep->extrapolateToPlane(state, plane, true);
  }
  catch (genfit::Exception& e) {
    std::cerr << e.what();

    delete rep;
    return kException;
  }


  // compare
  if (fabs(rep->getPos(state).Perp() - matRadius) > epsilonLen) {
      if (verbose) {
        origState.Print();
        state.Print();

        std::cerr << "radius difference = " << rep->getPos(state).Perp() - matRadius << "\n";

        std::cerr << std::endl;
      }
      delete rep;
      return kFailed;
    }

    delete rep;
    return kPassed;

}


e_testStatus checkErrorPropagation(bool writeHisto = false) {

  if (writeHisto)
    return kPassed;

  //double epsilonLen = 1.E-4; // 1 mu
  //double epsilonMom = 1.E-5; // 10 keV

  int pdg = randomPdg();
  genfit::AbsTrackRep* rep;
  rep = new genfit::RKTrackRep(pdg);

  //TVector3 pos(0,0,0);
  TVector3 pos(gRandom->Gaus(0,0.1),gRandom->Gaus(0,0.1),gRandom->Gaus(0,0.1));
  TVector3 mom(0,0.5,gRandom->Gaus(0,0.1));
  mom *= randomSign();


  genfit::MeasuredStateOnPlane state(rep);
  rep->setPosMom(state, pos, mom);

  genfit::SharedPlanePtr origPlane = state.getPlane();
  genfit::SharedPlanePtr plane(new genfit::DetPlane(TVector3(0,randomSign()*50,0), TVector3(0,randomSign()*1,0)));

  genfit::MeasuredStateOnPlane origState(state);

  // forth
  try {
    rep->extrapolateToPlane(state, plane);
  }
  catch (genfit::Exception& e) {
    std::cerr << e.what();

    delete rep;
    return kException;
  }


  // check
  if (isCovMatrix(state.getCov()) == kFailed) {
    if (verbose) {
      origState.Print();
      state.Print();
    }
    delete rep;
    return kFailed;
  }

  delete rep;
  return kPassed;

}


e_testStatus checkExtrapolateToLine(bool writeHisto = false) {

  if (writeHisto)
    return kPassed;

  double epsilonLen = 1.E-4; // 1 mu
  double epsilonMom = 1.E-5; // 10 keV

  int pdg = randomPdg();
  genfit::AbsTrackRep* rep;
  rep = new genfit::RKTrackRep(pdg);

  //TVector3 pos(0,0,0);
  TVector3 pos(0+gRandom->Gaus(0,0.1),0+gRandom->Gaus(0,0.1),0+gRandom->Gaus(0,0.1));
  TVector3 mom(0,0.5,0.);
  mom *= randomSign();


  genfit::StateOnPlane state(rep);
  rep->setPosMom(state, pos, mom);

  genfit::SharedPlanePtr origPlane = state.getPlane();
  genfit::StateOnPlane origState(state);

  TVector3 linePoint(gRandom->Gaus(),randomSign()*10+gRandom->Gaus(),gRandom->Gaus());
  TVector3 lineDirection(gRandom->Gaus(),gRandom->Gaus(),randomSign()*10+gRandom->Gaus());

  // forth
  try {
    rep->extrapolateToLine(state, linePoint, lineDirection, false);
  }
  catch (genfit::Exception& e) {
    std::cerr << e.what();

    delete rep;
    return kException;
  }


  // compare
  if (fabs(state.getPlane()->distance(linePoint)) > epsilonLen ||
      fabs(state.getPlane()->distance(linePoint+lineDirection)) > epsilonLen ||
      (rep->getMom(state).Unit() * state.getPlane()->getNormal()) > epsilonMom) {

      if (verbose) {
        origState.Print();
        state.Print();

        std::cout << "distance of linePoint to plane = " << state.getPlane()->distance(linePoint) << "\n";
        std::cout << "distance of linePoint+lineDirection to plane = "
                  << state.getPlane()->distance(linePoint + lineDirection) << "\n";
        std::cout << "direction * plane normal = " << rep->getMom(state).Unit() * state.getPlane()->getNormal() << "\n";
      }
      delete rep;
      return kFailed;
    }

    delete rep;
    return kPassed;

}


e_testStatus checkExtrapolateToPoint(bool writeHisto = false) {

  if (writeHisto)
    return kPassed;

  double epsilonLen = 1.E-4; // 1 mu
  double epsilonMom = 1.E-5; // 10 keV

  int pdg = randomPdg();
  genfit::AbsTrackRep* rep;
  rep = new genfit::RKTrackRep(pdg);

  //TVector3 pos(0,0,0);
  TVector3 pos(gRandom->Gaus(0,0.1),gRandom->Gaus(0,0.1),gRandom->Gaus(0,0.1));
  TVector3 mom(0,0.5,gRandom->Gaus(0,0.1));
  mom *= randomSign();


  genfit::StateOnPlane state(rep);
  rep->setPosMom(state, pos, mom);

  genfit::SharedPlanePtr origPlane = state.getPlane();
  genfit::StateOnPlane origState(state);

  TVector3 point(gRandom->Gaus(),randomSign()*10+gRandom->Gaus(),gRandom->Gaus());

  // forth
  try {
    rep->extrapolateToPoint(state, point, false);
  }
  catch (genfit::Exception& e) {
    std::cerr << e.what();

    delete rep;
    return kException;
  }


  // compare
  if (fabs(state.getPlane()->distance(point)) > epsilonLen ||
      fabs((rep->getMom(state).Unit() * state.getPlane()->getNormal())) - 1 > epsilonMom) {
      if (verbose) {
        origState.Print();
        state.Print();

        std::cout << "distance of point to plane = " << state.getPlane()->distance(point) << "\n";
        std::cout << "direction * plane normal = " << rep->getMom(state).Unit() * state.getPlane()->getNormal() << "\n";
      }
      delete rep;
      return kFailed;
    }

    delete rep;
    return kPassed;

}


e_testStatus checkExtrapolateToCylinder(bool writeHisto = false) {

  if (writeHisto)
    return kPassed;

  double epsilonLen = 1.E-4; // 1 mu
  //double epsilonMom = 1.E-5; // 10 keV

  int pdg = randomPdg();
  genfit::AbsTrackRep* rep;
  rep = new genfit::RKTrackRep(pdg);

  //TVector3 pos(0,0,0);
  TVector3 pos(0+gRandom->Gaus(0,0.1),0+gRandom->Gaus(0,0.1),0+gRandom->Gaus(0,0.1));
  TVector3 mom(0,0.5,0.);
  mom *= randomSign();


  genfit::StateOnPlane state(rep);
  rep->setPosMom(state, pos, mom);

  genfit::SharedPlanePtr origPlane = state.getPlane();
  genfit::StateOnPlane origState(state);

  const TVector3 linePoint(gRandom->Gaus(0,5), gRandom->Gaus(0,5), gRandom->Gaus(0,5));
  const TVector3 lineDirection(gRandom->Gaus(),gRandom->Gaus(),2+gRandom->Gaus());
  const double radius = gRandom->Uniform(10);

  // forth
  try {
    rep->extrapolateToCylinder(state, radius, linePoint, lineDirection, false);
  }
  catch (genfit::Exception& e) {
    std::cerr << e.what();

    delete rep;

    return kException;
  }

  TVector3 pocaOnLine(lineDirection);
  double t = 1./(pocaOnLine.Mag2()) * ((rep->getPos(state)*pocaOnLine) - (linePoint*pocaOnLine));
  pocaOnLine *= t;
  pocaOnLine += linePoint;

  TVector3 radiusVec = rep->getPos(state) - pocaOnLine;

  // compare
  if (fabs(state.getPlane()->getNormal()*radiusVec.Unit())-1 > epsilonLen ||
      fabs(lineDirection*radiusVec) > epsilonLen ||
      fabs(radiusVec.Mag()-radius) > epsilonLen) {
      if (verbose) {
        origState.Print();
        state.Print();

        std::cout << "lineDirection*radiusVec = " << lineDirection * radiusVec << "\n";
        std::cout << "radiusVec.Mag()-radius = " << radiusVec.Mag() - radius << "\n";
      }
      delete rep;
      return kFailed;
    }

    delete rep;
    return kPassed;

}


e_testStatus checkExtrapolateToSphere(bool writeHisto = false) {

  if (writeHisto)
    return kPassed;

  double epsilonLen = 1.E-4; // 1 mu
  //double epsilonMom = 1.E-5; // 10 keV

  int pdg = randomPdg();
  genfit::AbsTrackRep* rep;
  rep = new genfit::RKTrackRep(pdg);

  //TVector3 pos(0,0,0);
  TVector3 pos(0+gRandom->Gaus(0,0.1),0+gRandom->Gaus(0,0.1),0+gRandom->Gaus(0,0.1));
  TVector3 mom(0,0.5,0.);
  mom *= randomSign();


  genfit::StateOnPlane state(rep);
  rep->setPosMom(state, pos, mom);

  genfit::SharedPlanePtr origPlane = state.getPlane();
  genfit::StateOnPlane origState(state);

  const TVector3 centerPoint(gRandom->Gaus(0,10), gRandom->Gaus(0,10), gRandom->Gaus(0,10));
  const double radius = gRandom->Uniform(10);

  // forth
  try {
    rep->extrapolateToSphere(state, radius, centerPoint, false);
  }
  catch (genfit::Exception& e) {
    std::cerr << e.what();

    delete rep;

    return kException;
  }


  TVector3 radiusVec = rep->getPos(state) - centerPoint;

  // compare
  if (fabs(state.getPlane()->getNormal()*radiusVec.Unit())-1 > epsilonLen ||
      fabs(radiusVec.Mag()-radius) > epsilonLen) {
      if (verbose) {
        origState.Print();
        state.Print();

        std::cout << "state.getPlane()->getNormal()*radiusVec = " << state.getPlane()->getNormal() * radiusVec << "\n";
        std::cout << "radiusVec.Mag()-radius = " << radiusVec.Mag() - radius << "\n";
      }
      delete rep;
      return kFailed;
    }

    delete rep;
    return kPassed;

}


e_testStatus checkExtrapolateBy(bool writeHisto = false) {

  if (writeHisto)
    return kPassed;

  double epsilonLen = 1.E-3; // 10 mu

  int pdg = randomPdg();
  genfit::AbsTrackRep* rep;
  rep = new genfit::RKTrackRep(pdg);

  //TVector3 pos(0,0,0);
  TVector3 pos(0+gRandom->Gaus(0,0.1),0+gRandom->Gaus(0,0.1),0+gRandom->Gaus(0,0.1));
  TVector3 mom(0,0.5,0.);
  mom *= randomSign();


  genfit::StateOnPlane state(rep);
  rep->setPosMom(state, pos, mom);

  TVector3 posOrig(state.getPos());

  genfit::SharedPlanePtr origPlane = state.getPlane();
  genfit::StateOnPlane origState(state);

  double step = gRandom->Uniform(-15.,15.);
  double extrapolatedLen(0);

  // forth
  try {
    extrapolatedLen = rep->extrapolateBy(state, step, false);
  }
  catch (genfit::Exception& e) {
    return kException;
  }

  TVector3 posExt(state.getPos());




  // compare
  if (fabs(extrapolatedLen-step) > epsilonLen ||
      (posOrig - posExt).Mag() > fabs(step)) {
      if (verbose) {
        origState.Print();
        state.Print();

        std::cout << "extrapolatedLen-step = " << extrapolatedLen - step << "\n";
        std::cout << "started extrapolation from: ";
        posOrig.Print();
        std::cout << "extrapolated to ";
        posExt.Print();
        std::cout << "difference = " << (posOrig - posExt).Mag() << "; step = " << step << "; delta = "
                  << (posOrig - posExt).Mag() - fabs(step) << "\n";
      }
      delete rep;
      return kFailed;
    }

    delete rep;
    return kPassed;

}
//=====================================================================================================================
//=====================================================================================================================
//=====================================================================================================================
//=====================================================================================================================

struct TestCase {
  TestCase(std::string name, e_testStatus(*function)(bool)) : name_(name), function_(function), nPassed_(0), nFailed_(0), nException_(0) {;}
  void Print() {std::cout << name_ << " \t" << nPassed_ << " \t" << nFailed_ << " \t" << nException_ << "\n";}

  std::string name_;
  e_testStatus(*function_)(bool);
  unsigned int nPassed_;
  unsigned int nFailed_;
  unsigned int nException_;
};


int main() {

  const double BField = 15.;       // kGauss
  //const bool debug = true;

  gRandom->SetSeed(10);
  signal(SIGSEGV, handler);   // install our handler

  // init geometry and mag. field
  new TGeoManager("Geometry", "Geane geometry");
  TGeoManager::Import("genfitGeom.root");
  genfit::FieldManager::getInstance()->init(new genfit::ConstField(0.,0.,BField));
  genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());

  /*genfit::MaterialEffects::getInstance()->drawdEdx(2212);
  genfit::MaterialEffects::getInstance()->drawdEdx(11);
  genfit::MaterialEffects::getInstance()->drawdEdx(211);
  genfit::MaterialEffects::getInstance()->drawdEdx(13);
  return 0;*/

  TDatabasePDG::Instance()->GetParticle(211);

  const unsigned int nTests(1000);

  std::vector<TestCase> testCases;
  testCases.push_back(TestCase(std::string("checkSetGetPosMom()            "), &checkSetGetPosMom));
  testCases.push_back(TestCase(std::string("compareForthBackExtrapolation()"), &compareForthBackExtrapolation));
  testCases.push_back(TestCase(std::string("checkStopAtBoundary()          "), &checkStopAtBoundary));
  testCases.push_back(TestCase(std::string("checkErrorPropagation()        "), &checkErrorPropagation));
  testCases.push_back(TestCase(std::string("checkExtrapolateToLine()       "), &checkExtrapolateToLine));
  testCases.push_back(TestCase(std::string("checkExtrapolateToPoint()      "), &checkExtrapolateToPoint));
  testCases.push_back(TestCase(std::string("checkExtrapolateToCylinder()   "), &checkExtrapolateToCylinder));
  testCases.push_back(TestCase(std::string("checkExtrapolateToSphere()     "), &checkExtrapolateToSphere));
  testCases.push_back(TestCase(std::string("checkExtrapolateBy()           "), &checkExtrapolateBy));
  testCases.push_back(TestCase(std::string("compareForthBackJacNoise()     "), &compareForthBackJacNoise));


  for (unsigned int i=0; i<nTests; ++i) {

    for (unsigned int j=0; j<testCases.size(); ++j) {
      e_testStatus status = testCases[j].function_(false);
      switch (status) {
      case kPassed:
        testCases[j].nPassed_++;
        break;
      case kFailed:
        testCases[j].nFailed_++;
        std::cout << "failed " << testCases[j].name_ << " nr " << i << "\n";
        break;
      case kException:
        testCases[j].nException_++;
        std::cout << "exception at " << testCases[j].name_ << " nr " << i << "\n";
      }
    }

  }

  //CALLGRIND_STOP_INSTRUMENTATION;
  std::cout << "name                           " << " \t" << "pass" << " \t" << "fail" << " \t" << "exception" << "\n";
  for (unsigned int j=0; j<testCases.size(); ++j) {
    testCases[j].Print();
  }

  for (unsigned int j=0; j<testCases.size(); ++j) {
    testCases[j].function_(true);
  }

  return 0;
}


