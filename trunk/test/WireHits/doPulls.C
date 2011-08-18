void doPulls(TString infile,TString pdf){
  TFile::Open(infile);
  gROOT->ProcessLine(".L pullsOld.C");
  pulls p(t);
  p.Loop(pdf);
}
