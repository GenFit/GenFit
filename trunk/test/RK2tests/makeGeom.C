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

   TGeoVolume *top = gGeoManager->MakeBox("TOPPER", air, 500., 500., 500.);
   gGeoManager->SetTopVolume(top); // mandatory !
   
   TGeoVolume *slab = gGeoManager->MakeBox("slab", sil, 5., 5., .012);
   top->AddNode(slab, 1, new TGeoTranslation(0.,0., 0.3));
   top->AddNode(slab, 2, new TGeoTranslation(0.,0., 0.6));
   top->AddNode(slab, 3, new TGeoTranslation(0.,0., 0.9));
   top->AddNode(slab, 4, new TGeoTranslation(0.,0., 1.2));
   top->AddNode(slab, 5, new TGeoTranslation(0.,0., 1.5));
   //   top->AddNode(slab, 6, new TGeoTranslation(0.,0., 1.8));


   //--- close the geometry
   gGeoManager->CloseGeometry();

   //--- draw the ROOT box
   gGeoManager->SetVisLevel(10);
   //top->Draw("ogl");
   TFile *outfile = TFile::Open("genfitGeom.root","RECREATE");
   gGeoManager->Write();
   outfile->Close();
}
