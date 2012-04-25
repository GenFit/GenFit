void makeGeom()
{
  //--- Definition of a simple geometry
  //   gSystem->Load("libGeom");
   new TGeoManager("genfitGeom", "GENFIT geometry");
   gROOT->Macro("../../geometry/media.C");


   TGeoMedium *vacuum = gGeoManager->GetMedium("vacuum");
   assert(vacuum!=NULL);
   TGeoMedium *air = gGeoManager->GetMedium("air");
   assert(air!=NULL);
   TGeoMedium *sil = gGeoManager->GetMedium("silicon");
   assert(sil!=NULL);

   TGeoVolume *top = gGeoManager->MakeBox("TOPPER", vacuum, 1000., 1000., 1000.);
   gGeoManager->SetTopVolume(top); // mandatory !

   /*TGeoVolume *redBullCan = gGeoManager->MakeTube("redBullCan", sil, 3.-0.01, 3., 10.);//, 90., 270.);
   redBullCan->SetLineColor(kRed);
   top->AddNode(redBullCan, 1, gGeoIdentity);

   TGeoVolume *redBullCan2 = gGeoManager->MakeTube("redBullCan2", sil, 4.-0.01, 4., 10.);//, 90., 270.);
   redBullCan2->SetLineColor(kRed);
   top->AddNode(redBullCan2, 1, gGeoIdentity);*/

   //--- close the geometry
   gGeoManager->CloseGeometry();

   //--- draw the ROOT box
   gGeoManager->SetVisLevel(10);
   //top->Draw("ogl");
   TFile *outfile = TFile::Open("genfitGeom.root","RECREATE");
   gGeoManager->Write();
   outfile->Close();
}
