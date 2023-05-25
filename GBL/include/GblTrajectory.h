/*
 * GblTrajectory.h
 *
 *  Created on: Aug 18, 2011
 *      Author: kleinwrt
 */

/** \file
 *  GblTrajectory definition.
 *
 *  \author Claus Kleinwort, DESY, 2011 (Claus.Kleinwort@desy.de)
 *
 *  \copyright
 *  Copyright (c) 2011 - 2016 Deutsches Elektronen-Synchroton,
 *  Member of the Helmholtz Association, (DESY), HAMBURG, GERMANY \n\n
 *  This library is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU Library General Public License as
 *  published by the Free Software Foundation; either version 2 of the
 *  License, or (at your option) any later version. \n\n
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details. \n\n
 *  You should have received a copy of the GNU Library General Public
 *  License along with this program (see the file COPYING.LIB for more
 *  details); if not, write to the Free Software Foundation, Inc.,
 *  675 Mass Ave, Cambridge, MA 02139, USA.
 */

#ifndef GBLTRAJECTORY_H_
#define GBLTRAJECTORY_H_

#include "GblPoint.h"
#include "GblData.h"
#include "GblPoint.h"
#include "BorderedBandMatrix.h"
#include "MilleBinary.h"
#include "TMatrixDSymEigen.h"

//! Namespace for the general broken lines package
namespace gbl {

/// GBL trajectory.
  /**
   * List of GblPoints ordered by arc length.
   * Can be fitted and optionally written to MP-II binary file.
   */
  class GblTrajectory {
  public:
    explicit GblTrajectory(const std::vector<GblPoint>& aPointList, bool flagCurv = true,
                           bool flagU1dir = true, bool flagU2dir = true);
    GblTrajectory(const std::vector<GblPoint>& aPointList, unsigned int aLabel,
                  const TMatrixDSym& aSeed, bool flagCurv = true, bool flagU1dir =
                    true, bool flagU2dir = true);
    explicit GblTrajectory(
      const std::vector<std::pair<std::vector<GblPoint>, TMatrixD> >& aPointaAndTransList);
    GblTrajectory(
      const std::vector<std::pair<std::vector<GblPoint>, TMatrixD> >& aPointaAndTransList,
      const TMatrixD& extDerivatives, const TVectorD& extMeasurements,
      const TVectorD& extPrecisions);
    GblTrajectory(
      const std::vector<std::pair<std::vector<GblPoint>, TMatrixD> >& aPointaAndTransList,
      const TMatrixD& extDerivatives, const TVectorD& extMeasurements,
      const TMatrixDSym& extPrecisions);
    virtual ~GblTrajectory();
    bool isValid() const;
    unsigned int getNumPoints() const;
    unsigned int getResults(int aSignedLabel, TVectorD& localPar,
                            TMatrixDSym& localCov) const;
    unsigned int getMeasResults(unsigned int aLabel, unsigned int& numRes,
                                TVectorD& aResiduals, TVectorD& aMeasErrors, TVectorD& aResErrors,
                                TVectorD& aDownWeights);
    unsigned int getScatResults(unsigned int aLabel, unsigned int& numRes,
                                TVectorD& aResiduals, TVectorD& aMeasErrors, TVectorD& aResErrors,
                                TVectorD& aDownWeights);
    unsigned int getLabels(std::vector<unsigned int>& aLabelList);
    unsigned int getLabels(std::vector<std::vector<unsigned int> >& aLabelList);
    unsigned int fit(double& Chi2, int& Ndf, double& lostWeight,
                     std::string optionList = "");
    void milleOut(MilleBinary& aMille);
    void printTrajectory(unsigned int level = 0);
    void printPoints(unsigned int level = 0);
    void printData();
    std::vector<GblData> getData() {return theData;}

  private:
    unsigned int numAllPoints; ///< Number of all points on trajectory
    std::vector<unsigned int> numPoints; ///< Number of points on (sub)trajectory
    unsigned int numTrajectories; ///< Number of trajectories (in composed trajectory)
    unsigned int numOffsets; ///< Number of (points with) offsets on trajectory
    unsigned int numInnerTrans; ///< Number of inner transformations to external parameters
    unsigned int numCurvature; ///< Number of curvature parameters (0 or 1) or external parameters
    unsigned int numParameters; ///< Number of fit parameters
    unsigned int numLocals; ///< Total number of (additional) local parameters
    unsigned int numMeasurements; ///< Total number of measurements
    unsigned int externalPoint; ///< Label of external point (or 0)
    bool constructOK; ///< Trajectory has been successfully constructed (ready for fit/output)
    bool fitOK; ///< Trajectory has been successfully fitted (results are valid)
    std::vector<unsigned int> theDimension; ///< List of active dimensions (0=u1, 1=u2) in fit
    std::vector<std::vector<GblPoint> > thePoints; ///< (list of) List of points on trajectory
    std::vector<GblData> theData; ///< List of data blocks
    std::vector<unsigned int> measDataIndex; ///< mapping points to data blocks from measurements
    std::vector<unsigned int> scatDataIndex; ///< mapping points to data blocks from scatterers
    TMatrixDSym externalSeed; ///< Precision (inverse covariance matrix) of external seed
    std::vector<TMatrixD> innerTransformations; ///< Transformations at innermost points of
    // composed trajectory (from common external parameters)
    TMatrixD externalDerivatives; // Derivatives for external measurements of composed trajectory
    TVectorD externalMeasurements; // Residuals for external measurements of composed trajectory
    TVectorD externalPrecisions; // Precisions for external measurements of composed trajectory
    VVector theVector; ///< Vector of linear equation system
    BorderedBandMatrix theMatrix; ///< (Bordered band) matrix of linear equation system

    std::pair<std::vector<unsigned int>, TMatrixD> getJacobian(
      int aSignedLabel) const;
    void getFitToLocalJacobian(std::vector<unsigned int>& anIndex,
                               SMatrix55& aJacobian, const GblPoint& aPoint, unsigned int measDim,
                               unsigned int nJacobian = 1) const;
    void getFitToKinkJacobian(std::vector<unsigned int>& anIndex,
                              SMatrix27& aJacobian, const GblPoint& aPoint) const;
    void construct();
    void defineOffsets();
    void calcJacobians();
    void prepare();
    void buildLinearEquationSystem();
    void predict();
    double downWeight(unsigned int aMethod);
    void getResAndErr(unsigned int aData, double& aResidual,
                      double& aMeadsError, double& aResError, double& aDownWeight);
  };
}
#endif /* GBLTRAJECTORY_H_ */
