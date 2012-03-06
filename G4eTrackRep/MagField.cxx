/*
 * MagField.cc
 *
 *  Created on: Jan 26, 2012
 *      Author: poehler
 */

#include "MagField.h"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"


MagField::MagField()
: G4UniformMagField(G4ThreeVector()){
	// TODO Auto-generated constructor stub
	  GetGlobalFieldManager()->SetDetectorField(this);
	  GetGlobalFieldManager()->CreateChordFinder(this);
}

MagField::MagField(G4ThreeVector field)
: G4UniformMagField(G4ThreeVector(field)){
	// TODO Auto-generated constructor stub
	  GetGlobalFieldManager()->SetDetectorField(this);
	  GetGlobalFieldManager()->CreateChordFinder(this);
}


MagField::~MagField() {
	// TODO Auto-generated destructor stub
}

void MagField::SetFieldValue(G4double fieldValue)
{
  G4UniformMagField::SetFieldValue(G4ThreeVector(0,0,fieldValue));

}


//-------------------------------------------------------------
void MagField::SetFieldValue(G4ThreeVector fieldVector)
{
  // Find the Field Manager for the global field
  G4FieldManager* fieldMgr= GetGlobalFieldManager();

  if(fieldVector!=G4ThreeVector(0.,0.,0.))
  {
    G4UniformMagField::SetFieldValue(fieldVector);
    fieldMgr->SetDetectorField(this);
  } else {
    // If the new field's value is Zero, then it is best to
    //  insure that it is not used for propagation.
    G4MagneticField* magField = NULL;
    fieldMgr->SetDetectorField(magField);
  }
}



//-------------------------------------------------------------
G4FieldManager*  MagField::GetGlobalFieldManager()
{
  return G4TransportationManager::GetTransportationManager()->GetFieldManager();
}
