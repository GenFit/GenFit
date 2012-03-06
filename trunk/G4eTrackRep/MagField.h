/*
 * MagField.hh
 *
 *  Created on: Jan 26, 2012
 *      Author: poehler
 */

#ifndef MAGFIELD_HH_
#define MAGFIELD_HH_

#include "MagField.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"

class MagField: public G4UniformMagField {
public:
	MagField();
	MagField(G4ThreeVector);
	virtual ~MagField();

	void SetFieldValue (G4ThreeVector);
	void SetFieldValue (double);
	G4FieldManager*  GetGlobalFieldManager();
};

#endif /* MAGFIELD_HH_ */
