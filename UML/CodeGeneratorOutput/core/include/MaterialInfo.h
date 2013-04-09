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

#ifndef genfit_MaterialInfo_h
#define genfit_MaterialInfo_h

#include <memory>

#include "SharedPlanePtr.h"

//typedef std::shared_ptr<const genfit::DetPlane> sharedPlanePtr;

namespace genfit {

class MaterialInfo {

 private:

  sharedPlanePtr sharedPlane_; // Material boundary


};

} /* End of namespace genfit */

#endif // genfit_MaterialInfo_h
