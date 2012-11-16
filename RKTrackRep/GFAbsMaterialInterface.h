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

/** @addtogroup RKTrackRep
 * @{
 */

#ifndef GFABSMATERIALINTERFACE_H
#define GFABSMATERIALINTERFACE_H

#include <TObject.h>
#include <TVector3.h>


class GFAbsMaterialInterface : public TObject {
 public:
  GFAbsMaterialInterface(){;};
  virtual ~GFAbsMaterialInterface(){;};

  /** @brief Initialize the navigator at given position and with given direction
   */
  virtual void initTrack(const TVector3& pos, const TVector3& dir) = 0;

  /** @brief Get material parameters in current material
   */
  virtual void getMaterialParameters(double& density,
                                     double& Z,
                                     double& A,
                                     double& radiationLength,
                                     double& mEE) = 0;

  /** @brief Make a step until maxStep or the next boundary is reached.
   * After making a step to a boundary, the position has to be beyond the boundary,
   * i.e. in the current material has to be that beyond the boundary.
   * The actual step made is returned.
   */
  virtual double findNextBoundaryAndStep(double maxStep) = 0;

public:
  ClassDef(GFAbsMaterialInterface, 1);

};

#endif

/** @} */
