/* *************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//
// G3 Material Class for the G3 GUI
// Author: Andreas Morsch
// andreas.morsch@cern.ch
//

#include "G3Material.h"

ClassImp(G3Material)
G3Material::G3Material(char* name, char* title,
			       Float_t a, Float_t z, Float_t dens, Float_t radl, Float_t intl):
    TMaterial(name, title, a, z, dens, radl, intl)
{
    fId=-1;
}


void G3Material::Dump() const
{
    // Dump material information (Attn: overrides TObject::Dump()).

    printf("\n *****************************************");
    printf("\n Material Number:   %10d", fId);
    printf("\n %s", GetName());
    printf("\n Mass   Number:     %10.2f", fA);
    printf("\n Charge Number:     %10.2f", fZ);
    printf("\n Density:           %10.2f", fDensity);
    printf("\n Radiation  Length: %10.2f", fRadLength);
    printf("\n Absorption Length: %10.2f", fInterLength);        
}















