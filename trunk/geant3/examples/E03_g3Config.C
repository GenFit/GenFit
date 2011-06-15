// $Id: E03_g3Config.C 220 2007-11-19 16:08:06Z rdm $
//
// Configuration macro for Geant3 VirtualMC for Example03

void Config()
{
  if (strstr(gEnv->GetValue("TVirtualMC",""),"TGeant3TGeo")) {
     TGeant3TGeo* geant3 = new  TGeant3TGeo("C++ Interface to Geant3 using TGeo");
     cout << "TGeant3TGeo has been created." << endl;
  } else {
     TGeant3* geant3 = new  TGeant3("C++ Interface to Geant3");
     cout << "TGeant3 has been created." << endl;
  }

  geant3->SetDRAY(1);
  geant3->SetLOSS(1);
  geant3->SetHADR(0);
}


