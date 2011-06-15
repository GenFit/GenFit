
#include <GFException.h>
#include <GFDetPlane.h>
#include <GFTools.h>
#include <GFTrack.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TFitResult.h>
#include <TH1D.h>
#include <TMatrixT.h>
#include <TROOT.h>
#include <TStyle.h>
#include <vector>

#ifndef SMOOTHINGPULLS_H
#define SMOOTHINGPULLS_H

namespace SmoothingPulls {

	void displayPulls(std::vector<std::vector<GFTrack*>*> events, TApplication* rootapp);

}

#endif
