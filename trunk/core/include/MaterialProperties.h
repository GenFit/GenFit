/* Copyright 2008-2010, Technische Universitaet Muenchen,
   Authors: Christian Hoeppner & Sebastian Neubert & Johannes Rauch

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

#ifndef genfit_MaterialProperties_h
#define genfit_MaterialProperties_h

#include <TObject.h>


namespace genfit {


/**
 * @brief Material properties needed e.g. for material effects calculation.
 */
class MaterialProperties : public TObject {

 public:

  //! Compare material parameters
  friend bool operator== (const MaterialProperties& lhs, const MaterialProperties& rhs);
  friend bool operator!= (const MaterialProperties& lhs, const MaterialProperties& rhs);

  MaterialProperties();
  MaterialProperties(const double& density,
                     const double& Z,
                     const double& A,
                     const double& radiationLength,
                     const double& mEE);

  double getDensity() const {return density_;}
  double getZ() const {return Z_;}
  double getA() const {return A_;}
  double getRadLen() const {return radiationLength_;}
  double getMEE() const {return mEE_;}

  void getMaterialProperties(double& density,
                             double& Z,
                             double& A,
                             double& radiationLength,
                             double& mEE) const;

  void setMaterialProperties(const double& density,
                             const double& Z,
                             const double& A,
                             const double& radiationLength,
                             const double& mEE);

  void Print(const Option_t* = "") const;

 private:

  // material variables
  //! density of material
  double density_;
  //! Atomic number Z of material
  double Z_;
  //! Mass number A of material
  double A_;
  //! radiation length X0
  double radiationLength_;
  //! mean excitation energy [eV]
  double mEE_;


 public:
  ClassDef(MaterialProperties, 1)

};


inline MaterialProperties::MaterialProperties() :
  density_(0),
  Z_(0),
  A_(0),
  radiationLength_(0),
  mEE_(0)
{
  ;
}

inline MaterialProperties::MaterialProperties(const double& density,
                   const double& Z,
                   const double& A,
                   const double& radiationLength,
                   const double& mEE) :
  density_(density),
  Z_(Z),
  A_(A),
  radiationLength_(radiationLength),
  mEE_(mEE)
{
  ;
}

} /* End of namespace genfit */
/** @} */

#endif // genfit_MaterialProperties_h
