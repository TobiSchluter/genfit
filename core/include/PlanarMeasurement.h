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

#ifndef genfit_PlanarMeasurement_h
#define genfit_PlanarMeasurement_h

#include "AbsMeasurement.h"
#include "MeasurementOnPlane.h"


namespace genfit {

class AbsTrackRep;

/** @brief Measurement class implementing a planar hit geometry (1 or 2D).
 *
 *  @author Christian H&ouml;ppner (Technische Universit&auml;t M&uuml;nchen, original author)
 *  @author Sebastian Neubert  (Technische Universit&auml;t M&uuml;nchen, original author)
 *  @author Johannes Rauch  (Technische Universit&auml;t M&uuml;nchen, original author)
 *
 * The main feature of this type of hit is, that the detector plane
 * is defined by the detector hardware. 
 */
class PlanarMeasurement : public AbsMeasurement {

 public:
  PlanarMeasurement(int nDim = 1);
  PlanarMeasurement(const TVectorD& rawHitCoords, const TMatrixDSym& rawHitCov, int detId, int hitId, TrackPoint* trackPoint);

  virtual SharedPlanePtr constructPlane(const StateOnPlane* state) const override;

  virtual MeasurementOnPlane constructMeasurementOnPlane(const AbsTrackRep*, const SharedPlanePtr) const override;

  virtual const TMatrixD& getHMatrix(const AbsTrackRep*) const override;

  virtual void setPlane(SharedPlanePtr physicalPlane, int planeId = -1) {physicalPlane_ = physicalPlane; planeId_ = planeId;}

 protected:
  SharedPlanePtr physicalPlane_;
  int planeId_;
};

} /* End of namespace genfit */

#endif // genfit_PlanarMeasurement_h
