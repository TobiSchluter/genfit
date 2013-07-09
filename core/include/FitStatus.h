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
    isFitted(false), isFitConverged_(false), trackChanged_(false), charge_(0) {;}

  virtual ~FitStatus() {};

  bool isFitted() const {return isFitted_;}
  bool isFitConverged() const {return isFitConverged_;}
  bool hasTrackChanged() const {return hasTrackChanged_;}
  double getCharge() const {return charge_;}

  void setIsFitted(bool isFitted = true) {isFitted_ = isFitted;}
  void setIsFitConverged(bool isFitConverged = true) {isFitConverged_ = isFitConverged_;}
  void setHasTrackChanged(bool hasTrackChanged = true) {hasTrackChanged_ = hasTrackChanged_;}
  void setCharge(double charge) {charge_ = charge;}

 protected:

  bool isFitted_; // has the track been fitted
  bool isFitConverged_; // did the fit converge
  bool hasTrackChanged_; // has anything in the Track been changed sinde the fit -> fit isn't valid anymore
  double charge_; // fitted charge

  double fChi2_; // chi^2 of the forward fit
  double bChi2_; // chi^2 of the backward fit
  double fNdf_; // degrees of freedom of the forward fit
  double bNdf_; // degrees of freedom of the backward fit


} /* End of namespace genfit */
/** @} */

#endif // genfit_FitStatus_h
