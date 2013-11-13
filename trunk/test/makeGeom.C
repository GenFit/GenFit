void makeGeom()
{
  //--- Definition of a simple geometry
  //   gSystem->Load("libGeom");
   new TGeoManager("genfitGeom", "GENFIT geometry");


   unsigned int medInd(0);
   Double_t mPar[10];
   //TAG sollte wieder 0 werden sens flag
   mPar[0]=0.;//sensitive volume flag
   mPar[1]=1.;//magnetic field flag
   mPar[2]=30.;//max fiel in kGauss
   mPar[3]=0.1;//maximal angular dev. due to field
   mPar[4]=0.01;//max step allowed (in cm)
   mPar[5]=1.e-5;//max fractional energy loss
   mPar[6]=1.e-3;//boundary crossing accuracy
   mPar[7]=1.e-5;//minimum step
   mPar[8]=0.;//not defined
   mPar[9]=0.;//not defined

   TGeoMaterial *_siliconMat = new TGeoMaterial("siliconMat",28.0855,14.,2.329);
   _siliconMat->SetRadLen(1.);//calc automatically, need this for elemental mats.
   TGeoMedium *_silicon = new TGeoMedium("silicon",medInd++,_siliconMat,mPar);

   TGeoMixture *_airMat = new TGeoMixture("airMat",3);
   _airMat->AddElement(14.01,7.,.78);
   _airMat->AddElement(16.00,8.,.21);
   _airMat->AddElement(39.95,18.,.01);
   _airMat->SetDensity(1.2e-3);
   TGeoMedium *_air = new TGeoMedium("air",medInd++,_airMat,mPar);

   TGeoMixture *_vacuumMat = new TGeoMixture("vacuumMat",3);
   _vacuumMat->AddElement(14.01,7.,.78);
   _vacuumMat->AddElement(16.00,8.,.21);
   _vacuumMat->AddElement(39.95,18.,.01);
   _vacuumMat->SetDensity(1.2e-15);
   TGeoMedium *_vacuum = new TGeoMedium("vacuum",medInd++,_vacuumMat,mPar);




   TGeoMedium *vacuum = gGeoManager->GetMedium("vacuum");
   assert(vacuum!=NULL);
   TGeoMedium *air = gGeoManager->GetMedium("air");
   assert(air!=NULL);
   TGeoMedium *sil = gGeoManager->GetMedium("silicon");
   assert(sil!=NULL);

   TGeoVolume *top = gGeoManager->MakeBox("TOPPER", air, 1000., 1000., 1000.);
   gGeoManager->SetTopVolume(top); // mandatory !

   double thickness(0.05);
   double distance = 1;

   for (unsigned int i=1; i<5; ++i){
     TGeoVolume *redBullCan = gGeoManager->MakeTube("redBullCan", sil, i*distance, i*distance+thickness, 20.);//, 90., 270.);
     redBullCan->SetLineColor(kRed);
     top->AddNode(redBullCan, 1, gGeoIdentity);
   }


   //--- close the geometry
   gGeoManager->CloseGeometry();

   //--- draw the ROOT box
   gGeoManager->SetVisLevel(10);
   //top->Draw("ogl");
   TFile *outfile = TFile::Open("genfitGeom.root","RECREATE");
   gGeoManager->Write();
   outfile->Close();
}
