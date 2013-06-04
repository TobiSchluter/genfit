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

#include "AbsMeasurement.h"

#include <assert.h>


namespace genfit {

AbsMeasurement::AbsMeasurement(const TVectorD& rawHitCoords, const TMatrixDSym& rawHitCov, int detId, int hitId, TrackPoint* trackPoint)
  : rawHitCoords_(rawHitCoords), rawHitCov_(rawHitCov), detId_(detId), hitId_(hitId), trackPoint_(trackPoint)
{
  assert(rawHitCov_.GetNrows() == rawHitCoords_.GetNrows());
}


AbsMeasurement::~AbsMeasurement()
{
  ;
}

} /* End of namespace genfit */
