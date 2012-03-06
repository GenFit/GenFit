//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: DetectorConstruction.hh,v 1.1 2010-10-18 15:56:17 maire Exp $
// GEANT4 tag $Name: geant4-09-04-patch-02 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
/*
#include "core/GFAbsTrackRep.h"
#include "RKTrackRep/RKTrackRep.h"
#include "core/GFTrack.h"

#include "Geant4GM/include/Geant4GM/volumes/Factory.h"
#include "RootGM/include/RootGM/volumes/Factory.h"
#include "TGeoManager.h"*/

#include <cstdlib>

class G4Box;
class G4Tubs;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4UniformMagField;
//class DetectorMessenger;
//class MySensitiveDetector;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
    DetectorConstruction();
  //  DetectorConstruction(bool);
//    DetectorConstruction(std::vector<GFTrack>*);
  //  DetectorConstruction(std::vector<GFTrack>*,bool);
   ~DetectorConstruction();

  public:
     
    /* void SetAbsorberMaterial (G4String);
     void SetAbsorberThickness(G4double);     

     void SetGapMaterial (G4String);     
     void SetGapThickness(G4double);
     
     void SetCalorSizeYZ(G4double);          
     void SetNbOfLayers (G4int);   
      */
    // void SetMagField(G4double);
     
     G4VPhysicalVolume* Construct();

     void UpdateGeometry();
//     void setWritefile(bool);
  //   bool getWritefile();
     bool file;
  public:
  
  //   void PrintCalorParameters();
                    
     G4double GetWorldSizeX()           {return WorldSizeX;}; 
     G4double GetWorldSizeYZ()          {return WorldSizeYZ;};
     
  /*   G4double GetCalorThickness()       {return CalorThickness;};
     G4double GetCalorSizeYZ()          {return CalorSizeYZ;};
      
     G4int GetNbOfLayers()              {return NbOfLayers;}; 
     
     G4Material* GetAbsorberMaterial()  {return AbsorberMaterial;};
     G4double    GetAbsorberThickness() {return AbsorberThickness;};      
     
     G4Material* GetGapMaterial()       {return GapMaterial;};
     G4double    GetGapThickness()      {return GapThickness;};
     
*/     const G4VPhysicalVolume* GetphysiWorld() {return physiWorld;}
//	   const double GetField ();
//     const G4VPhysicalVolume* GetAbsorber()   {return physiAbsorber;};
//     const G4VPhysicalVolume* GetGap()        {return physiGap;};
                 
  private:
     
    G4Material*        AbsorberMaterial;
    /*  G4double           AbsorberThickness;
     
     G4Material*        GapMaterial;
     G4double           GapThickness;
     
     G4int              NbOfLayers;
     G4double           LayerThickness;
          
     G4double           CalorSizeYZ;
     G4double           CalorThickness;
  */
     G4Material*        defaultMaterial;
     G4double           WorldSizeYZ;
     G4double           WorldSizeX;
            
     G4Box*             solidWorld;    //pointer to the solid World 
     G4LogicalVolume*   logicWorld;    //pointer to the logical World
     G4VPhysicalVolume* physiWorld;    //pointer to the physical World

     G4Tubs*			solidBlock;    //pointer to the solid Block
     G4LogicalVolume*   logicBlock;    //pointer to the logical Block
     G4VPhysicalVolume* physiBlock;    //pointer to the physical Block

     //bool isSensitive;
//     MySensitiveDetector* pSensitivePart;
    // std::vector<GFTrack>* detTracks;

  /*
     G4Box*             solidLayer;    //pointer to the solid Layer 
     G4LogicalVolume*   logicLayer;    //pointer to the logical Layer
     G4VPhysicalVolume* physiLayer;    //pointer to the physical Layer
         
     G4Box*             solidAbsorber; //pointer to the solid Absorber
     G4LogicalVolume*   logicAbsorber; //pointer to the logical Absorber
     G4VPhysicalVolume* physiAbsorber; //pointer to the physical Absorber
     
     G4Box*             solidGap;      //pointer to the solid Gap
     G4LogicalVolume*   logicGap;      //pointer to the logical Gap
     G4VPhysicalVolume* physiGap;      //pointer to the physical Gap
  */
    // G4UniformMagField* magField;      //pointer to the magnetic field
     
//     DetectorMessenger* detectorMessenger;  //pointer to the Messenger
      
  private:
    
     void DefineMaterials();
    // void ComputeCalorParameters();
     G4VPhysicalVolume* ConstructBlock();
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
/*
inline void DetectorConstruction::ComputeCalorParameters()
{
  // Compute derived parameters of the calorimeter
     LayerThickness = AbsorberThickness + GapThickness;
     CalorThickness = NbOfLayers*LayerThickness;
     
     WorldSizeX = 1.2*CalorThickness; WorldSizeYZ = 1.2*CalorSizeYZ;
}
*/
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

