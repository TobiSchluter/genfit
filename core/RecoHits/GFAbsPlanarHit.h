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

#ifndef GFABSPLANARHIT_H
#define GFABSPLANARHIT_H

#include <iostream>
#include <string>

#include "TMatrixT.h"
#include "TObject.h"

#include "GFAbsRecoHit.h"

/** @brief Policy class implementing a planar hit geometry. 
 *
 *  @author Christian H&ouml;ppner (Technische Universit&auml;t M&uuml;nchen, original author)
 *  @author Sebastian Neubert  (Technische Universit&auml;t M&uuml;nchen, original author)
 *  @author Johannes Rauch  (Technische Universit&auml;t M&uuml;nchen, original author)
 * 
 * RecoHits for planar detectors should inherit 
 * from GFAbsPlanarHit.
 *
 * The main feature of this type of hit is, that the detector plane
 * is completely defined by the detector hardware. Derived RecoHits need only
 * to supply the physical detector plane from their geometry database.
 */

class GFAbsPlanarHit : public GFAbsRecoHit {
public:

  // Constructors/Destructors ---------
  virtual ~GFAbsPlanarHit();
  GFAbsPlanarHit(){}

  GFAbsPlanarHit(int dim) : GFAbsRecoHit(dim){;}
  

  // Accessors -----------------------
  
  /** @brief Returns the physical detector plane.
   */
  const GFDetPlane& getDetPlane(GFAbsTrackRep*) {return fPhysicalDetPlane;}
  
  virtual void getMeasurement(const GFAbsTrackRep* rep,
                              const GFDetPlane& pl,
                              const TMatrixT<double>& statePred,
                              const TMatrixT<double>& covPred,
                              TMatrixT<double>& m,
                              TMatrixT<double>& V) {
    static_cast<void>(rep);
    static_cast<void>(statePred);
    static_cast<void>(covPred);
    m.ResizeTo(fHitCoord);
    V.ResizeTo(fHitCov);
    m = fHitCoord;
    V = fHitCov;
  }

  // Modifiers -----------------------

  /** @brief Set physical detector plane. 
   * Needs to be called before hit can be used.
   *
   * For the planar detector the detector plane is fixed by the geometry of the
   * detector hardware. This method should be called in the constructor of
   * any derived RecoHit in order to setup the geometry of this hit.
   */
  void setDetPlane(const GFDetPlane& p){fPhysicalDetPlane = p;}

 protected:

  // Private Data Members ------------
  
  /** @brief Physical detector plane. Given by detector hardware.
   */
  GFDetPlane fPhysicalDetPlane;

  // Private Methods -----------------

 public:
  ClassDef(GFAbsPlanarHit,1)

};

#endif

/** @} */
