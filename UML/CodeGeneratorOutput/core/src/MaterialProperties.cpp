/* Copyright 2008-2009, Technische Universitaet Muenchen,
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

#include "MaterialProperties.h"

namespace genfit {

MaterialProperties::MaterialProperties() :
  density_(0),
  Z_(0),
  A_(0),
  radiationLength_(0),
  mEE_(0)
{
  ;
}

MaterialProperties::MaterialProperties(const double& density,
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
