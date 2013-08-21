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

#ifndef genfit_FitStatus_h
#define genfit_FitStatus_h

#include <iostream>
#include <Rtypes.h>

namespace genfit {

class FitStatus {

 public:

  FitStatus() :
    isFitted_(false), isFitConverged_(false), trackHasChanged_(false), trackIsPruned_(false), charge_(0) {;}

  virtual ~FitStatus() {};

  virtual FitStatus* clone() const {return new FitStatus(*this);}

  bool isFitted() const {return isFitted_;}
  bool isFitConverged() const {return isFitConverged_;}
  bool hasTrackChanged() const {return trackHasChanged_;}
  bool isTrackPruned() const {return trackIsPruned_;}
  double getCharge() const {return charge_;}

  void setIsFitted(bool fitted = true) {isFitted_ = fitted;}
  void setIsFitConverged(bool fitConverged = true) {isFitConverged_ = fitConverged;}
  void setHasTrackChanged(bool trackChanged = true) {trackHasChanged_ = trackChanged;}
  void setIsTrackPruned(bool pruned = true) {trackIsPruned_ = pruned;}
  void setCharge(double charge) {charge_ = charge;}

  void Print(const Option_t* = "") const {
    std::cout << "fitStatus \n";
    if (isFitted_) {
      std::cout << " track has been fitted,";
      if (isFitConverged_)
        std::cout << " fit has converged,";
      else
        std::cout << " fit has NOT converged,";
      if (trackHasChanged_) std::cout << " track has changed since the fit,";
      if (trackIsPruned_) std::cout << " track is pruned,";
      std::cout << " fitted charge = " << charge_ << " \n";
    }
    else
      std::cout << " track has NOT been fitted,";
  }

 protected:

  bool isFitted_; // has the track been fitted
  bool isFitConverged_; // did the fit converge
  bool trackHasChanged_; // has anything in the Track been changed sinde the fit -> fit isn't valid anymore
  bool trackIsPruned_; // Information has been stripped off, no refitting possible!
  double charge_; // fitted charge

  ClassDef(FitStatus,1);

};

} /* End of namespace genfit */
/** @} */

#endif // genfit_FitStatus_h
