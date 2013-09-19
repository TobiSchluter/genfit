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

#ifndef genfit_TrackCandHit_h
#define genfit_TrackCandHit_h

#include <TObject.h>


namespace genfit {

/**
 * @brief Hit object for use in TrackCand. Provides IDs and sorting parameters.
 */
class TrackCandHit : public TObject {
 public:

  virtual TrackCandHit* clone() const {return new TrackCandHit(*this);}

  // Constructors/Destructors ---------
  TrackCandHit(int detId   = -1,
               int hitId   = -1,
               int planeId = -1,
               double sortingParameter  =  0.);

  ~TrackCandHit();

  /** @brief Equality operator. Does not check sortingParameter.
   */
  friend bool operator== (const TrackCandHit& lhs, const TrackCandHit& rhs);
  friend bool operator!= (const TrackCandHit& lhs, const TrackCandHit& rhs) {
    return !(lhs == rhs);
  }

  /** @brief Compare sortingParameter, needed for sorting
   */
  friend bool operator< (const TrackCandHit& lhs, const TrackCandHit& rhs) {
    return (lhs.sortingParameter_ < rhs.sortingParameter_);
  }

  // Accessors
  int    getDetId() const {return detId_;}
  int    getHitId() const {return hitId_;}
  int    getPlaneId() const {return planeId_;}
  double getSortingParameter() const {return sortingParameter_;}

  virtual void Print(Option_t* option = "") const;

  // Modifiers
  void setSortingParameter(double sortingParameter) {sortingParameter_ = sortingParameter;}

 protected:
  // Private Data Members ------------
  int    detId_; // detId id is -1 per default
  int    hitId_; // hitId id is -1 per default
  int    planeId_; // planeId id is -1 per default
  double sortingParameter_; // sorting parameter


 public:

  ClassDef(TrackCandHit,1)

};

} /* End of namespace genfit */
/** @} */

#endif // genfit_TrackCandHit_h
