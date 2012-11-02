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
 * @{ */



#ifndef GFWIREPOINTHIT_H
#define GFWIREPOINTHIT_H

#include "TMatrixT.h"
#include "TObject.h"

#include "GFDetPlane.h"

class GFAbsRecoHit;
class GFAbsTrackRep;

/** @brief  class for hits in wire detectors (STT and DCH)
 *  which can measure the coordinate along the wire
 *  @author Christian H&ouml;ppner (Technische Universit&auml;t M&uuml;nchen, original author)
 *  @author Lia Lavezzi (INFN Pavia, original author)
 *  @author Sebastian Neubert  (Technische Universit&auml;t M&uuml;nchen, original author)
 * This  is not valid for any kind of plane orientation
 * choice: to use it you MUST choose a plane described by u 
 * and v axes with v coincident with the wire (and u orthogonal
 * to it, obviously).
 * The hit will be described by 8 coordinates:
 * w_x1, w_y1, w_z1, w_x2, w_y2, w_z2, rdrift, zreco
 * where w_ji (with j = x, y, z and i = 1, 2) are the wire
 * extremities coordinates; rdrift = distance from the wire (u 
 * coordinate in the plane) and zreco = coordinate along the
 * wire (in the plane reference frame, v coordinate).
 *
 */
class GFWirepointHit : public GFWireHit {
public:


  GFWirepointHit();
 
  /** @brief Hit coordinates in detector plane.
   */
  TMatrixT<double> hitCoord(GFAbsRecoHit*,const GFDetPlane&);

  /** @brief Hit covariances in detector plane.
   */
  TMatrixT<double> hitCov(GFAbsRecoHit*,const GFDetPlane&);

  /** @brief Check if the detector plane is valid
   */
  void checkPlane(GFAbsRecoHit*,const GFDetPlane&);

  virtual ~GFWirepointHit(){;}

 private:
  // Private Data Members ------------
  GFDetPlane fDetPlane;

 public:
  ClassDef(GFWirepointHit,0);

};

#endif

/** @} */
