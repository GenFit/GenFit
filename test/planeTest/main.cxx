#include<iostream>

#include "GFDetPlane.h"
#include "GFRectFinitePlane.h"
#include"TFile.h"

int main(){
  GFDetPlane d(TVector3(0.0,0.0,1.),TVector3(1.,1.,0.),TVector3(-1.,1.,0.),
	       new GFRectFinitePlane(-2,2,-2,2));

  for(int i=-3;i<=3;++i){
    for(int j=-3;j<=3;++j){
      std::cout << i+0.1 << " " << j+0.1 << " " << d.inActive(i+0.1,j+0.1) << std::endl;
    }
  }

  TVector3 pos(1.,0.,0.);
  TVector3 mom(1.,0.,0.);
  d.straightLineToPlane(pos,mom).Print();

  TFile *f = TFile::Open("plane.root","RECREATE");

  d.Print();
  d.Write("DP");
  f->Close();


}


