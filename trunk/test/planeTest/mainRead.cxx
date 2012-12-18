#include <iostream>

#include "GFDetPlane.h"
#include "TFile.h"
#include "TROOT.h"

int main(){
  TFile *f = TFile::Open("plane.root");

  GFDetPlane *pl = (GFDetPlane*) gROOT->FindObject("DP");
  pl->Print();
  
  for(int i=-3;i<=3;++i){
    for(int j=-3;j<=3;++j){
      std::cout << i+0.1 << " " << j+0.1 << " " << pl->inActive(i+0.1,j+0.1) << std::endl;
    }
  }
  


}


