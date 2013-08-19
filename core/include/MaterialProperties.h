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

#include "TObject.h"

namespace genfit {

  class MaterialProperties {

 public:

  /** Compares material parameters, but not segmentLength_ */
  friend bool operator== (const MaterialProperties& lhs, const MaterialProperties& rhs);
  friend bool operator!= (const MaterialProperties& lhs, const MaterialProperties& rhs);

  MaterialProperties();
  MaterialProperties(const double& density,
                     const double& Z,
                     const double& A,
                     const double& radiationLength,
                     const double& mEE,
                     const double& segmentLength);

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
  double density_; // density of material
  double Z_; // Atomic number Z of material
  double A_; // Mass number A of material
  double radiationLength_; // radiation length
  double mEE_; // mean excitation energy [eV]

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
                   const double& mEE,
                   const double& segmentLength) :
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
