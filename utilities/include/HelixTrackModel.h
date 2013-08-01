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
/**
 *  @author Johannes Rauch (Technische Universit&auml;t M&uuml;nchen, original author)
 *
 */


/** @addtogroup utilities
 * @{
 */

#ifndef genfit_HelixTrackModel_h
#define genfit_HelixTrackModel_h

#include <TObject.h>
#include <TVector3.h>


namespace genfit {

/** @brief Detector plane genfit geometry class
 *
 * A detector plane is the principle object to define coordinate systems for
 * track fitting in genfit. Since a particle trajectory is a
 * one-dimensional object (regardless of any specific parameterization)
 * positions with respect to the track are always measured in a plane.
 *
 * Which plane is chosen depends on the type of detector. Fixed plane
 * detectors have their detector plane defined by their mechanical setup. While
 * wire chambers or time projection chambers might want to define a detector
 * plane more flexibly.
 *
 * This class parameterizes a plane in terms of an origin vector o
 * and two plane-spanning directions u and v.
 */

class HelixTrackModel : public TObject {

  /**
   * Helix track model for testing purposes
   */

 public:

  // Constructors/Destructors ---------
  HelixTrackModel(const TVector3& pos, const TVector3& mom, double charge);

  void getPosMom(double tracklength, TVector3& pos, TVector3& mom) const;
  void getPosDir(double tracklength, TVector3& pos, TVector3& dir) const {
    getPosMom(tracklength, pos, dir);
    dir.SetMag(1);
  }


 private:

  double sgn_;
  double mom_;
  double R_; // radius
  TVector3 center_;
  double alpha0_;
  double theta_;


 public:
  ClassDef(HelixTrackModel,1)

};

} /* End of namespace genfit */
/** @} */

#endif // genfit_HelixTrackModel_h
