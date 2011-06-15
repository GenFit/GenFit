#define pulls_cxx
#include "pulls.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void pulls::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L pulls.C
//      Root > pulls t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   gROOT->SetStyle("Plain");
   /*
     gStyle->SetTextFont(62);
     gStyle->SetTextSize(16);
     gStyle->SetLabelFont(62,"xyz");
     gStyle->SetLabelSize(16,"xyz");
     gStyle->SetTitleFont(62,"xyz");
     gStyle->SetTitleSize(16,"xyz");
   */

   Long64_t nentries = fChain->GetEntriesFast();

   TH1D *hmomPu = new TH1D("hmomPu","mom pull",100,-6.,6.);
   TH1D *hqopPu = new TH1D("hqopPu","q/p pull",100,-6.,6.);
   TH1D *hupPu = new TH1D("hupPu","u pull",100,-6.,6.);
   TH1D *hvpPu = new TH1D("hvpPu","v pull",100,-6.,6.);
   TH1D *huPu = new TH1D("huPu","u' pull",100,-6.,6.);
   TH1D *hvPu = new TH1D("hvPu","v' pull",100,-6.,6.);
   double momTruth = 0.5;
   TH1D *hmomRe = new TH1D("hmomRe","hmomRe",100,0.85*momTruth,1.15*momTruth) ;
   TH1D *hfail = new TH1D("hfail","hfail",20,0.,20.);
   TH1D *hchi2 = new TH1D("chi2","chi2",100,0.,6.);
   Long64_t nbytes = 0, nb = 0;

   double yMax(100);
   hmomPu->GetYaxis()->SetRangeUser(0.,yMax);
   hqopPu->GetYaxis()->SetRangeUser(0.,yMax);
   huPu->GetYaxis()->SetRangeUser(0.,yMax);
   hvPu->GetYaxis()->SetRangeUser(0.,yMax);
   hupPu->GetYaxis()->SetRangeUser(0.,yMax);
   hvpPu->GetYaxis()->SetRangeUser(0.,yMax);

   const int iqop=0;
   const int iup=1;
   const int ivp=2;
   const int iu=3;
   const int iv=4;
   const int iqopT=0;
   const int iupT=1;
   const int ivpT=2;
   const int iuT=3;
   const int ivT=4;

   TFile *outfile = TFile::Open("histo.root","RECREATE");

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      //      stMCT->Print();
      //      stREC->Print();
      double invmom = (*stREC)[iqop][0];
      double sigmasqustate = (*covREC)[iqop][iqop];
      if(sigmasqustate<1.E-16) continue;
      double sigma_p = 1/pow(invmom,4.) * sigmasqustate;
      double momSi=TMath::Sqrt(sigma_p);
      double momRe=fabs(1./((*(stREC))[iqop][0]));
      double momTr=fabs(1./((*(stMCT))[iqopT][0]));
      double momPu=(momRe-momTr)/momSi;
      
      double qopSi=TMath::Sqrt((*covREC)[iqop][iqop]);
      double qopRe=(*stREC)[iqop][0];
      double qopTr=(*stMCT)[iqopT][0];
      double qopPu=(qopRe-qopTr)/qopSi;

      double upSi=TMath::Sqrt((*covREC)[iup][iup]);
      double upRe=(*stREC)[iup][0];
      double upTr=(*stMCT)[iupT][0];
      double upPu=(upRe-upTr)/upSi;
      double vpSi=TMath::Sqrt((*covREC)[ivp][ivp]);
      double vpRe=(*stREC)[ivp][0];
      double vpTr=(*stMCT)[ivpT][0];
      double vpPu=(vpRe-vpTr)/vpSi;
      
      double uSi=TMath::Sqrt((*covREC)[iu][iu]);
      double uRe=(*stREC)[iu][0];
      double uTr=(*stMCT)[iuT][0];
      double uPu=(uRe-uTr)/uSi;
      double vSi=TMath::Sqrt((*covREC)[iv][iv]);
      double vRe=(*stREC)[iv][0];
      double vTr=(*stMCT)[ivT][0];
      double vPu=(vRe-vTr)/vSi;
      
      //      std::cout << "momRe " << momRe << std::endl;
      hmomRe->Fill(momRe);
      hmomPu->Fill(momPu);
      hqopPu->Fill(qopPu);
      hupPu->Fill(upPu);
      hvpPu->Fill(vpPu);
      huPu->Fill(uPu);
      hvPu->Fill(vPu);
      hfail->Fill(nfail);
      hchi2->Fill(chi2);
   }

   hmomPu->Write();
   hqopPu->Write();
   huPu->Write();
   hvPu->Write();
   hupPu->Write();
   hvpPu->Write();
   hchi2->Write();
   outfile->Close();

   gStyle->SetOptFit(1111);
   TCanvas *cpulls = new TCanvas("cpulls","cpulls",1100,700);
   cpulls->Divide(3,2);
   cpulls->cd(1);
   hqopPu->Fit("gaus");
   hqopPu->Draw();
   cpulls->cd(2);
   hupPu->Fit("gaus");
   hupPu->Draw();
   cpulls->cd(3);
   hvpPu->Fit("gaus");
   hvpPu->Draw();
   cpulls->cd(4);
   huPu->Fit("gaus");
   huPu->Draw();
   cpulls->cd(5);
   hvPu->Fit("gaus");
   hvPu->Draw();
   cpulls->cd(6);
   hmomRe->Fit("gaus");
   hmomRe->Draw();
   //cpulls->cd(7);
   //hfail->Draw();
   //cpulls->cd(8);
   //hchi2->Draw();
}
