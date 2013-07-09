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

#include <TObject.h>

namespace genfit {

class FitStatus : public TObject {

 public:

  FitStatus() :
    isFitted_(false), isFitConverged_(false), hasTrackChanged_(false), charge_(0) {;}

  virtual ~FitStatus() {};

  virtual FitStatus* clone() const {return new FitStatus(*this);}

  bool isFitted() const {return isFitted_;}
  bool isFitConverged() const {return isFitConverged_;}
  bool hasTrackChanged() const {return hasTrackChanged_;}
  double getCharge() const {return charge_;}

  void setIsFitted(bool fitted = true) {isFitted_ = fitted;}
  void setIsFitConverged(bool fitConverged = true) {isFitConverged_ = fitConverged;}
  void setHasTrackChanged(bool trackChanged = true) {hasTrackChanged_ = trackChanged;}
  void setCharge(double charge) {charge_ = charge;}

  void Print(const Option_t* = "") const {
    std::cout << "fitStatus \n";
    if (isFitted_) {
      std::cout << " track has been fitted,";
      if (isFitConverged_) std::cout << " fit has converged,";
      if (hasTrackChanged_) std::cout << " but track has changed,";
      if (isFitConverged_) {
        std::cout << " fitted charge = " << charge_ << " \n";
      }
    }
    else
      std::cout << " track has NOT been fitted,";
  }

 protected:

  bool isFitted_; // has the track been fitted
  bool isFitConverged_; // did the fit converge
  bool hasTrackChanged_; // has anything in the Track been changed sinde the fit -> fit isn't valid anymore
  double charge_; // fitted charge

};

} /* End of namespace genfit */
/** @} */

#endif // genfit_FitStatus_h
