/* Copyright 2008-2009, Technische Universitaet Muenchen,
   Authors: Christian Hoeppner & Sebastian Neubert

   This file is part of GENFIT.

   GENFIT is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published
   by the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   GENFIT is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with GENFIT.  If not, see <http://www.gnu.org/licenses/>.
*/
//-----------------------------------------------------------
// File and Version Information:
// $Id$
//
// Description:
//      Abstract Interface for the Nystrom system equation
//       u''=g(u,u')
//
//
// Environment:
//      Software developed for the PANDA Detector at FAIR.
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------

#ifndef ABSNYSTROMEQM_HH
#define ABSNYSTROMEQM_HH

// Base Class Headers ----------------


// Collaborating Class Headers -------
#include <ostream> // remove if you do not need streaming op

// Collaborating Class Declarations --
#include "TVectorT.h"

class AbsNystromEQM {
public:

  // Constructors/Destructors ---------
  AbsNystromEQM(){;}
  virtual ~AbsNystromEQM(){;}

  virtual TVectorT<double> eval(const TVectorT<double>&u,
				const TVectorT<double>& uprim, 
				const TVectorT<double>& par) = 0; 

};

#endif

//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
