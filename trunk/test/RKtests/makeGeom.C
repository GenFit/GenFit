void makeGeom()
{
  //--- Definition of a simple geometry
  //   gSystem->Load("libGeom");
   new TGeoManager("genfitGeom", "GENFIT geometry");
   gROOT->Macro("../../geometry/media.C");


   TGeoMedium *air = gGeoManager->GetMedium("air");
   assert(air!=NULL);
   TGeoMedium *vacuum = gGeoManager->GetMedium("vacuum");
   assert(vacuum!=NULL);
   TGeoMedium *sil = gGeoManager->GetMedium("silicon");
   assert(sil!=NULL);

   TGeoVolume *top = gGeoManager->MakeBox("TOPPER", vacuum, 500., 500., 500.);
   gGeoManager->SetTopVolume(top); // mandatory !

   //--- close the geometry
   gGeoManager->CloseGeometry();

   //--- draw the ROOT box
   gGeoManager->SetVisLevel(10);
   //top->Draw("ogl");
   TFile *outfile = TFile::Open("genfitGeom.root","RECREATE");
   gGeoManager->Write();
   outfile->Close();
}
