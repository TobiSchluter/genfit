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

#include "StepLimits.h"

#include <algorithm>
#include <assert.h>
#include <limits>

namespace genfit {


bool pairCompare(std::pair<StepLimitType, double> i, std::pair<StepLimitType, double> j) {
  return i.second < j.second;
}


StepLimits::StepLimits() :
  stepSign_(1)
{
  ;
}


double StepLimits::getLimit(StepLimitType type) const {
  auto it = limits_.find(type);

  if (it == limits_.end()) {
    return std::numeric_limits<double>::max();
  }

  return it->second;
}


std::pair<StepLimitType, double> StepLimits::getLowestLimit(double margin) const {
  if (limits_.size() == 0) {
    return std::pair<StepLimitType, double>(stp_noLimit, std::numeric_limits<double>::max());
  }

  auto itMedium = limits_.upper_bound(stp_planeRough);
  auto itHard   = limits_.upper_bound(stp_sMax);

  // case 1: only soft limits
  if (itMedium == limits_.end() && itHard == limits_.end()) {
    return *min_element(limits_.begin(), limits_.end(), pairCompare );
  }

  // find minimum medium limit
  auto itMinMedium = *min_element(itMedium, itHard, pairCompare );

  // case 2: medium limits, no hard limits -> ignore soft limits
  if (itHard == limits_.end()) {
    return itMinMedium;
  }

  // find minimum hard limit
  auto itMinHard = *min_element(itHard, limits_.end(), pairCompare );

  // case 3: hard limits -> ignore soft limits, lowest hard limit may exceed lowest soft limit by up to #margin
  if (itMinHard.second <= (1+margin)*itMinMedium.second)
    return itMinHard;

  return itMinMedium;

}


void StepLimits::setLimit(StepLimitType type, double value) {
  assert (type != stp_noLimit);
  limits_[type] = value;
}


void StepLimits::reduceLimit(StepLimitType type, double value) {
  assert (type != stp_noLimit);
  value = fabs(value);

  std::map<StepLimitType, double>::iterator it;
  it = limits_.find(type);

  if (it == limits_.end()) {
    limits_[type] = value;
  }
  else {
    if (value < it->second)
      it->second = value;
  }
}


void StepLimits::setStepSign(char signedVal) {
  if (signedVal < 0)
    stepSign_ = -1;
  else
    stepSign_ = 1;
}

void StepLimits::setStepSign(double signedVal) {
  if (signedVal < 0.)
    stepSign_ = -1;
  else
    stepSign_ = 1;
}


void StepLimits::Print() {
  /*std::cout << "Stepsize has been limited due to following reasons: \n";
  if (sMax) std::cout << " sMax: stepsize only limited by SMax defined in #estimateStep() \n";
  if (sMaxArg) std::cout << " sMaxArg: stepsize limited by argument maxStepArg passed to #estimateStep()\n";
  if (planeDist) std::cout << " planeDist: stepsize limited due to first estimation of SL distance to destination plane\n";
  if (fieldCurv) std::cout << " fieldCurv: stepsize limited by curvature and magnetic field inhomogenities\n";
  if (minStep) std::cout << " minStep: stepsize set to minimum value by stepper\n";
  if (boundary) std::cout << " boundary: stepsize limited by stepper because material boundary is encountered\n";
  if (momLoss) std::cout << " boundary: stepsize limited by stepper because maximum momLoss is reached\n";
  if (atPlane) std::cout << " atPlane: stepsize limited because destination plane is reached\n";*/
}

} /* End of namespace genfit */
