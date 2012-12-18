rootlogon()
{  
  if(isLibrary("libgenfit"))gSystem->Load("libgenfit");
  if(isLibrary("libgenfitRK"))gSystem->Load("libgenfitRK");
  if(isLibrary("libgenfitGeane"))gSystem->Load("libgenfitGeane");
  if(isLibrary("libmyRecoHitExamples"))gSystem->Load("libmyRecoHitExamples");
  if(isLibrary("libgfrave"))gSystem->Load("libgfrave");
}

Bool_t isLibrary(const char* libName)
{
  if (TString(gSystem->DynamicPathName(libName, kTRUE)) != TString(""))
    return kTRUE;
  else  
    return kFALSE;
}    

