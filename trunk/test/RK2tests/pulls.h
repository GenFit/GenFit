//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jan 21 16:07:18 2010 by ROOT version 5.24/00
// from TTree t/example output
// found on file: out.root
//////////////////////////////////////////////////////////

#ifndef pulls_h
#define pulls_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class pulls {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   TMatrixT<double> *stMCT;
   TMatrixT<double> *covMCT;
   TMatrixT<double> *stREC;
   TMatrixT<double> *covREC;
   Double_t        chi2;
   Int_t           ndf;
   Int_t           nfail;

   // List of branches
   TBranch        *b_stMCT;   //!
   TBranch        *b_covMCT;   //!
   TBranch        *b_stREC;   //!
   TBranch        *b_covREC;   //!
   TBranch        *b_chi2;   //!
   TBranch        *b_ndf;   //!
   TBranch        *b_nfail;   //!

   pulls(TTree *tree=0);
   virtual ~pulls();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef pulls_cxx
pulls::pulls(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("out.root");
      if (!f) {
         f = new TFile("out.root");
      }
      tree = (TTree*)gDirectory->Get("t");

   }
   Init(tree);
}

pulls::~pulls()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t pulls::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t pulls::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void pulls::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   stMCT = 0;
   covMCT = 0;
   stREC = 0;
   covREC = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("stMCT", &stMCT, &b_stMCT);
   fChain->SetBranchAddress("covMCT", &covMCT, &b_covMCT);
   fChain->SetBranchAddress("stREC", &stREC, &b_stREC);
   fChain->SetBranchAddress("covREC", &covREC, &b_covREC);
   fChain->SetBranchAddress("chi2", &chi2, &b_chi2);
   fChain->SetBranchAddress("ndf", &ndf, &b_ndf);
   fChain->SetBranchAddress("nfail", &nfail, &b_nfail);
   Notify();
}

Bool_t pulls::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void pulls::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t pulls::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef pulls_cxx
