//-*-mode: C++; c-basic-offset: 2; -*-
/* Copyright 2013-2014
 *  Authors: Sergey Yashchenko and Tadeas Bilka
 *
 *  This is an interface to General Broken Lines
 * 
 *  CHK Feb 2023: try to implement GBLv3 (with thick scatterers) and ambiguities (CDC)
 *  - GBL extentions:
 *    + GblData::GblData()
 *    + std::vector<GblData> GblTrajectory::getData()
 *  - not yet working:
 *    + debug: debugLvl_ not set (from steering), ostream debugOut not appearing in cout
 * 
 *  Version: 6 --------------------------------------------------------------
 *  - complete rewrite to GblFitter using GblFitterInfo and GblFitStatus
 *  - mathematics should be the same except for additional iterations 
 *    (not possible before)
 *  - track is populated with scatterers + additional points with 
 *    scatterers and no measurement (optional)
 *  - final track contains GblFitStatus and fitted states from GBL prediction
 *  - 1D/2D hits supported (pixel, single strip, combined strips(2D), wire)
 *  - At point: Only the very first raw measurement is used and from
 *    that, constructed MeasurementOnPlane with heighest weight
 *  Version: 5 --------------------------------------------------------------
 *  - several bug-fixes:
 *    - Scatterers at bad points
 *    - Jacobians at a point before they should be (code reorganized)q
 *    - Change of sign of residuals
 *  Version: 4 --------------------------------------------------------------
 *    Fixed calculation of equvivalent scatterers (solution by C. Kleinwort)
 *    Now a scatterer is inserted at each measurement (except last) and
 *    between each two measurements.
 *    TrueHits/Clusters. Ghost (1D) hits ignored. With or
 *    without magnetic field.
 *  Version: 3 --------------------------------------------------------------
 *    This version now supports both TrueHits and Clusters for VXD.
 *    It can be used for arbitrary material distribution between
 *    measurements. Moments of scattering distribution are computed
 *    and translated into two equivalent thin GBL scatterers placed
 *    at computed positions between measurement points.
 *  Version: 2 ... never published -----------------------------------------
 *    Scatterer at each boundary (tooo many scatterers). TrueHits/Clusters.
 *    Without global der.&MP2 output.
 *  Version: 1 --------------------------------------------------------------
 *    Scatterers at measurement planes. TrueHits
 *  Version: 0 --------------------------------------------------------------
 *    Without scatterers. Genfit 1.
 *  -------------------------------------------------------------------------
 *
 *  This file is part of GENFIT.
 *
 *  GENFIT is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published
 *  by the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  GENFIT is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with GENFIT.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "GblFitter2.h"
#include "../include/GblFitStatus.h"
#include "GblFitterInfo2.h"
#include "GblTrajectory.h"
#include "GblPoint.h"
#include "ICalibrationParametersDerivatives.h"

#include "Track.h"
#include <TFile.h>
#include <TH1F.h>
#include <TTree.h>
#include <string>
#include <list>
#include <FieldManager.h>
#include <HMatrixU.h>
#include <HMatrixV.h>
#include <Math/SMatrix.h>
#include <TMatrixD.h>
#include <TVectorDfwd.h>
#include <TMatrixT.h>
#include <TVector3.h>
#include "IO.h"

using namespace gbl;
using namespace std;
using namespace genfit;

/**
 * Destructor
 */
GblFitter2::~GblFitter2() {
}

void GblFitter2::processTrackWithRep(Track* trk, const AbsTrackRep* rep, bool resortHits)
{
  cleanGblInfo(trk, rep);
  
  if (resortHits)
    sortHits(trk, rep);
  
  // This flag enables/disables fitting of q/p parameter in GBL
  // It is switched off automatically if no B-field at (0,0,0) is detected.
  bool fitQoverP = true;
  //TODO: Use clever way to determine zero B-field
  double Bfield = genfit::FieldManager::getInstance()->getFieldVal(TVector3(0., 0., 0.)).Mag();
  if (!(Bfield > 1.e-16))
    fitQoverP = false;
  // degrees of freedom after fit
  int Ndf = 0;
  // Chi2 after fit
  double Chi2 = 0.;
  //FIXME: d-w's not used so far...
  double lostWeight = 0.;  

  // Preparation of points (+add reference states) for GBL fit
  // -----------------------------------------------------------------
  genfit::GblFitStatus* gblfs = new genfit::GblFitStatus();
  trk->setFitStatus(gblfs, rep);
  gblfs->setCurvature(fitQoverP);
  //
  // Propagate reference seed, create scattering points, calc Jacobians
  // and store everything in fitter infos. (ready to collect points and fit)
  //
  //
  gblfs->setTrackLen(
  //
    constructGblInfo(trk, rep)
  //
  );
  //
  //
  gblfs->setIsFittedWithReferenceTrack(true);
  gblfs->setNumIterations(0); //default value, still valid, No GBL iteration
  if (m_externalIterations < 1)
    return;
  // -----------------------------------------------------------------
  
  // cppcheck-suppress unreadVariable
  unsigned int nFailed = 0;
  // cppcheck-suppress unreadVariable
  int fitRes = 0;
  std::vector<std::string> gblIterations;
  gblIterations.push_back(m_gblInternalIterations);

  // Iterations and updates of fitter infos and fit status
  // ------------------------------------------------------------------- 
  for (unsigned int iIter = 0; iIter < m_externalIterations; iIter++) {
    // GBL refit (1st of reference, then refit of GBL trajectory itself)
 
    if (iIter > 0) {
      // This repropagates to get new Jacobians and Noise,
      // if planes changed, predictions are extrapolated to new planes
      GblFitterInfo2* prevFitterInfo = 0;
      GblFitterInfo2* currFitterInfo = 0;
      for (unsigned int ip = 0; ip < trk->getNumPoints(); ip++) {
        if (trk->getPoint(ip)->hasFitterInfo(rep) && 
          (currFitterInfo = dynamic_cast<GblFitterInfo2*>(trk->getPoint(ip)->getFitterInfo(rep)))) {
          currFitterInfo->recalculateJacobian(prevFitterInfo);
          prevFitterInfo = currFitterInfo;
        }
      }
    }
    
    int nscat = 0, nmeas = 0, ndummy = 0, nall = 0;
    std::vector<gbl::GblPoint> points = collectGblPoints(trk, rep, (iIter < m_resolveAmbiguties));
    for(unsigned int ip = 0;ip<points.size(); ip++) {
      GblPoint & p = points.at(ip);
      if (p.getScatDim())
        nscat++;
      if (p.numMeasurements()) {
        nmeas++;
        nall += p.numMeasurements(); }
      if(!p.numMeasurements()&&!p.getScatDim())
        ndummy++;
    }
    if (debugLvl_ >= 10)
      /*debugOut*/ cout << " GblTrajectory " << points.size() << " " << nscat << " " << nmeas << " " << nall << " " << ndummy << std::endl;

    //
    // Ambiguities ?
    //
    if (nall > nmeas) {
      gbl::GblTrajectory trajAmbig(points, gblfs->hasCurvature());
      fitRes = trajAmbig.fit(Chi2, Ndf, lostWeight, "hh"); // 2 Huber iterations
      if (debugLvl_ >= 10)
        /*debugOut*/ cout << " GBLfitRes ambig. " << fitRes << "; " << Chi2 << " " << Ndf << " " << lostWeight << std::endl;
      if (!fitRes) {
        //
        // resolve (1D) ambiguities (and update down-weights)
        //
        unsigned int nResolved = 0;
        unsigned int nSwapped =0;
        // Update points in fitterInfo(rep) and GblPoint
        for (unsigned int ip = 0; ip < trk->getNumPoints(); ip++) {      
          GblFitterInfo2 * gblfi = dynamic_cast<GblFitterInfo2*>(trk->getPoint(ip)->getFitterInfo(rep));
          if (gblfi) gblfi->resolveAmbiguities(trajAmbig, points.at(ip), nResolved, nSwapped);    
        } 
        if (debugLvl_ >= 15)
          /*debugOut*/ cout << " ambitities resolved " << nResolved << ", swapped " << nSwapped << std::endl;
        // last iteration with resolving ambiguities: update down-weights
        if (iIter == m_resolveAmbiguties-1) updateGblDownweights(trajAmbig, trk, rep);
      }
    }
    
    gbl::GblTrajectory traj(points, gblfs->hasCurvature());
    //traj.printTrajectory(1);
    //traj.printPoints(1);
    //traj.printData();
    
    fitRes = traj.fit(Chi2, Ndf, lostWeight, m_gblInternalIterations);
//traj.printTrajectory(1000);
//traj.printPoints(1000);
    //if (debugLvl_ >= 1)
      /*debugOut*/ cout << " GBLfitRes-4 " << iIter << " " << fitRes << "; " << Chi2 << " " << Ndf << " " << lostWeight << std::endl; 
    // Update down-weights
    if (m_gblInternalIterations.size() > 0) updateGblDownweights(traj, trk, rep);

    // Update fit results in fitterinfos
    updateGblInfo(traj, trk, rep);
        
    gblfs->setIsFitted(true);
    gblfs->setIsFitConvergedPartially(fitRes == 0);
    nFailed = trk->getNumPointsWithMeasurement() - nmeas;
    gblfs->setNFailedPoints(nFailed);
    gblfs->setIsFitConvergedFully(fitRes == 0 && nFailed == 0);
    gblfs->setNumIterations(iIter + 1);
    gblfs->setChi2(Chi2);    
    gblfs->setNdf(Ndf);
    gblfs->setCharge(trk->getFittedState().getCharge());
    //std::cout << " gblfs " << fitRes << " " << nFailed << " " << iIter << " " << Chi2 << " " << Ndf << std::endl; 
    
//    #ifdef
    int npoints_meas = trk->getNumPointsWithMeasurement();  
    int npoints_all = trk->getNumPoints();
        
    cout << "-------------------------------------------------------" << endl;
    cout << "               GBL processed genfit::Track            " << endl;
    cout << "-------------------------------------------------------" << endl;
    cout << " # Track Points       :  " << npoints_all  << endl;
    cout << " # Meas. Points       :  " << npoints_meas << endl;
    cout << " # GBL points all     :  " << traj.getNumPoints();
    if (ndummy)
      cout << " (" << ndummy << " dummy) ";
    cout << endl;
    cout << " # GBL points meas    :  " << nmeas << endl;
    cout << " # GBL points scat    :  " << nscat << endl;    
    cout << "-------------- GBL Fit Results ----------- Iteration  " << iIter+1 << " " << ((iIter < gblIterations.size()) ? gblIterations[iIter] : "") << endl;
    cout << " Fit q/p parameter    :  " << (gblfs->hasCurvature() ? ("True") : ("False")) << endl;
    cout << " Valid trajectory     :  " << ((traj.isValid()) ? ("True") : ("False")) << endl;
    cout << " Fit result           :  " << fitRes << "    (0 for success)" << endl;
    cout << " GBL track NDF        :  " << Ndf << "    (-1 for failure)" << endl;
    cout << " GBL track Chi2       :  " << Chi2 << endl;
    cout << " GBL track P-value    :  " << TMath::Prob(Chi2, Ndf) << endl;
    cout << "-------------------------------------------------------" << endl;
//    #endif
    
  }  
  // -------------------------------------------------------------------

}

void GblFitter2::cleanGblInfo(Track* trk, const AbsTrackRep* rep) const {
  
  for (int ip = trk->getNumPoints() - 1; ip >=0; ip--) {
    trk->getPoint(ip)->setScatterer(nullptr); 
    trk->getPoint(ip)->deleteFitterInfo(rep);
    //TODO
    if (!trk->getPoint(ip)->hasRawMeasurements())
      trk->deletePoint(ip);
  }
}

void GblFitter2::sortHits(Track* trk, const AbsTrackRep* rep) const { 
  // All measurement points in ref. track
  int npoints_meas = trk->getNumPointsWithMeasurement();  
  // Prepare state for extrapolation of track seed
  StateOnPlane reference(rep);
  rep->setTime(reference, trk->getTimeSeed());
  rep->setPosMom(reference, trk->getStateSeed());
  // Take the state to first plane
  SharedPlanePtr firstPlane(trk->getPointWithMeasurement(0)->getRawMeasurement(0)->constructPlane(reference));
  reference.extrapolateToPlane(firstPlane);
  //1st point is at arc-len=0
  double arcLenPos = 0;
  
  // Loop only between meas. points 
  for (int ipoint_meas = 0; ipoint_meas < npoints_meas - 1; ipoint_meas++) {
    // current measurement point
    TrackPoint* point_meas = trk->getPointWithMeasurement(ipoint_meas);    
    // Current detector plane
    SharedPlanePtr plane = point_meas->getRawMeasurement(0)->constructPlane(reference);    
    // Get the next plane
    SharedPlanePtr nextPlane(trk->getPointWithMeasurement(ipoint_meas + 1)->getRawMeasurement(0)->constructPlane(reference));    
    
    point_meas->setSortingParameter(arcLenPos);
    arcLenPos += reference.extrapolateToPlane(nextPlane);
    
  } // end of loop over track points with measurement
  trk->getPointWithMeasurement(npoints_meas - 1)->setSortingParameter(arcLenPos);
  trk->sort();
}

std::vector<gbl::GblPoint> GblFitter2::collectGblPoints(genfit::Track* trk, const genfit::AbsTrackRep* rep, bool allowAmbiguities) {
  //TODO store collected points in in fit status? need streamer for GblPoint (or something like that)
  std::vector<gbl::GblPoint> thePoints;
  thePoints.clear();
  
  // Collect points from track and fitterInfo(rep)
  for (unsigned int ip = 0; ip < trk->getNumPoints(); ip++) {   
    GblFitterInfo2 * gblfi = dynamic_cast<GblFitterInfo2*>(trk->getPoint(ip)->getFitterInfo(rep));
    if (!gblfi)
      continue;
    gblfi->setLabel(ip+1);
    thePoints.push_back(gblfi->constructGblPoint(ip, allowAmbiguities));      
  }  
  return thePoints;
}

void GblFitter2::updateGblDownweights(gbl::GblTrajectory& traj, genfit::Track* trk, const genfit::AbsTrackRep* rep) {
  //FIXME
  if (!traj.isValid())
    return;
  
  // Update fitterInfo(rep)
  for (unsigned int ip = 0; ip < trk->getNumPoints(); ip++) {      
    GblFitterInfo2 * gblfi = dynamic_cast<GblFitterInfo2*>(trk->getPoint(ip)->getFitterInfo(rep));
    if (gblfi) gblfi->updateDownweights(traj);    
  }
}

void GblFitter2::updateGblInfo(gbl::GblTrajectory& traj, genfit::Track* trk, const genfit::AbsTrackRep* rep) {
  //FIXME
  if (!traj.isValid())
    return;
  
  // Update points in track and fitterInfo(rep)
  for (unsigned int ip = 0; ip < trk->getNumPoints(); ip++) {      
    GblFitterInfo2 * gblfi = dynamic_cast<GblFitterInfo2*>(trk->getPoint(ip)->getFitterInfo(rep));
    if (!gblfi)
      continue;
    
    // The point will calculate its position on the track
    gblfi->updateFitResults(traj);
    
  }
}

double GblFitter2::constructGblInfo(Track* trk, const AbsTrackRep* rep)
{ 
  // All measurement points in ref. track
  int npoints_meas = trk->getNumPointsWithMeasurement();  
  
  // Prepare state for extrapolation of track seed
  // Take the state to first plane
  StateOnPlane reference(rep);
  rep->setTime(reference, trk->getTimeSeed());
  rep->setPosMom(reference, trk->getStateSeed());

  SharedPlanePtr firstPlane(trk->getPointWithMeasurement(0)->getRawMeasurement(0)->constructPlane(reference));  
  reference.extrapolateToPlane(firstPlane);
  
  double sumTrackLen = 0;
  
  // Loop only between meas. points 
  for (int ipoint_meas = 0; ipoint_meas < npoints_meas; ipoint_meas++) {
    // current measurement point
    TrackPoint* point_meas = trk->getPointWithMeasurement(ipoint_meas);     
    
    // First (starting) point?
    if (ipoint_meas == 0) {
      //std::cout << " pointMeas " << ipoint_meas << " " << sumTrackLen << std::endl;
      // should define jacobian=1, noise=0
      GblFitterInfo2* gblfimeas(new GblFitterInfo2(point_meas, rep, reference));
      point_meas->setFitterInfo(gblfimeas);
      continue;
    }
        
    // Current detector plane
    SharedPlanePtr plane = point_meas->getRawMeasurement(0)->constructPlane(reference);   
    
    // Extrapolate from previous measurement  
    double nextStep = rep->extrapolateToPlane(reference, plane, false, true);
    
    if (0. > nextStep) {
      //if (debugLvl_ > 0)
        /*debugOut*/ cout << " ERROR: The extrapolation to measurement point " << ipoint_meas << " stepped back by " << nextStep << "cm !!! Track will be cut before this point." << endl;
      // stop trajectory construction here
      break;
    }          
   
    // Track length up to  point
    sumTrackLen += nextStep; 
    
    if (debugLvl_ >= 100)
      /*debugOut*/ cout << " pointMeas " << ipoint_meas << " " << sumTrackLen << std::endl;
    GblFitterInfo2* gblfimeas(new GblFitterInfo2(point_meas, rep, reference));
    point_meas->setFitterInfo(gblfimeas);   

  } // end of loop over track points with measurement
  return sumTrackLen;
}
