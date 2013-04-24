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

/** @addtogroup RKTrackRep
 * @{
 */

#ifndef genfit_StepLimits_h
#define genfit_StepLimits_h

#include <map>

namespace genfit {

enum StepLimitType {
  // soft limits (only rough estimation, can go beyond safely)
  stp_noLimit,    // only for internal use
  stp_planeRough, // stepsize limited due to first estimation of SL distance to destination plane

  // medium limits (can go a bit further if e.g. plane or boundary will be reached)
  stp_fieldCurv,  // stepsize limited by curvature and magnetic field inhomogenities
  stp_momLoss,    // stepsize limited by stepper because maximum momLoss is reached
  stp_sMax,       // stepsize limited by SMax defined in #estimateStep()

  // hard limits (must stop there at any case!)
  stp_sMaxArg,    // stepsize limited by argument maxStepArg passed to #estimateStep()
  stp_boundary,   // stepsize limited by stepper because material boundary is encountered
  stp_plane       // stepsize limited because destination plane is reached
};

bool pairCompare( std::pair<StepLimitType, double> i, std::pair<StepLimitType, double> j) {
  return i.second < j.second;
}

class StepLimits {

 public:
  StepLimits();

  /**
   * Get limit of #type. If that limit has not yet been set, return max double value.
   */
  double getLimit(StepLimitType type) const;
  double getLimitSigned(StepLimitType type) const {
    return stepSign_*getLimit(type);
  }
  /**
   * Get the lowest limit. If medium limits are there, soft limits will be ignored.
   * If hard limits are there, medium limits can be exceeded by up to #margin
   * (default margin is 0.1, i.e. medium limits can be exceeded by up to 10%).
   * If no limit has been set yet, return std::pair<stp_noLimit, std::numeric_limits<double>::max>.
   */
  std::pair<StepLimitType, double> getLowestLimit(double margin = 0.) const;
  /**
   * Get the unsigned numerical value of the lowest limit
   */
  double getLowestLimitVal(double margin = 0.) const {
    return getLowestLimit(margin).second;
  }
  /**
   * Get the numerical value of the lowest limit, signed with #stepSign_.
   */
  double getLowestLimitSignedVal(double margin = 0.) const {
    return getLowestLimit(margin).second * stepSign_;
  }

  char getStepSign() const {return stepSign_;} // +- 1

  /**
   * absolute of #value will be taken! If limit is already lower, it will stay.
   */
  void reduceLimit(StepLimitType type, double value);
  /**
   * absolute of value will be taken! If limit is already lower, it will be set to value anyway.
   */
  void setLimit(StepLimitType type, double value);
  /**
   * sets #stepSign_ to sign of #signedVal
   */
  void setStepSign(char signedVal);
  /**
   * sets #stepSign_ to sign of #signedVal
   */
  void setStepSign(double signedVal);

  void removeLimit(StepLimitType type) {limits_.erase(type);}

  void reset();
  void Print();

 private:
  std::map<StepLimitType, double> limits_; // limits are unsigned (i.e. non-negative)
  char stepSign_;

};

} /* End of namespace genfit */

#endif // genfit_StepLimits_h
