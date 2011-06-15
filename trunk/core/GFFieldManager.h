/* Copyright 2008-2010, Technische Universitaet Muenchen,
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
/** @addtogroup genfit
 * @{
 */


#ifndef GFFIELDMANAGER_H
#define GFFIELDMANAGER_H

#include"GFAbsBField.h"
#include<iostream>

/** @brief Singleton which provides access to magnetic field for track representations
 *
 *  @author Christian H&ouml;ppner (Technische Universit&auml;t M&uuml;nchen, original author)
 *  @author Sebastian Neubert  (Technische Universit&auml;t M&uuml;nchen, original author)
 * 
 */

class GFFieldManager{
 private:
  GFFieldManager(){}
  static GFFieldManager* fInstance;
  static GFAbsBField* fField;

 public:
  GFAbsBField* getField(){
    if(fField==NULL){
      std::cerr << "Appareantly GFFieldManager hasnt been initialized with a correct GFAbsBField pointer -> abort" << std::endl;
      throw;
    }
    return fField;
  }

  static TVector3 getFieldVal(const TVector3& x){
    if(fInstance==NULL){
      std::cerr << "Appareantly GFFieldManager hasnt been instantiated yet, call getInstance() and init() before getFieldVal() -> abort" << std::endl;
      throw;
    }
    if(fField==NULL){
      std::cerr << "Appareantly GFFieldManager hasnt been initialized with a correct GFAbsBField pointer -> abort" << std::endl;
      throw;
    }
    return fField->get(x);
  }

  //! set the magntic field here. Magnetic field classes must be derived from GFAbsBField
  void init(GFAbsBField* b) {
    fField=b;
  }

  static GFFieldManager* getInstance(){
    if(fInstance==NULL) {
      fInstance = new GFFieldManager();
    }
    return fInstance;
  }


};

/** @} */ 
#endif
