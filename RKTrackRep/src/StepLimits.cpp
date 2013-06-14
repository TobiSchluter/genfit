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
#include <iostream>
#include <limits>

#include <math.h>

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

  auto itMedium = limits_.upper_bound(stp_noLimit);
  auto itHard   = limits_.upper_bound(stp_sMax);

  // find minimum medium limit
  auto itMinMedium = *min_element(itMedium, itHard, pairCompare );

  // case 2: medium limits, no hard limits
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
  limits_[type] = fabs(value);
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
  for (auto it = limits_.begin(); it != limits_.end(); ++it) {
    std::cout << "   | " << it->second << " cm due to ";
    switch (it->first) {
    case stp_noLimit:
      break;
    case stp_fieldCurv:
      std::cout << "stp_fieldCurv (medium limit): stepsize limited by curvature and magnetic field inhomogenities";
      break;
    case stp_momLoss:
      std::cout << "stp_momLoss (medium limit): stepsize limited by stepper because maximum momLoss is reached";
      break;
    case stp_sMax:
      std::cout << "stp_sMax (medium limit): stepsize limited by SMax defined in #estimateStep()";
      break;
    case stp_sMaxArg:
      std::cout << "stp_sMaxArg (hard limit): stepsize limited by argument maxStepArg passed to #estimateStep()";
      break;
    case stp_boundary:
      std::cout << "stp_boundary (hard limit): stepsize limited by stepper because material boundary is encountered";
      break;
    case stp_plane:
      std::cout << "stp_plane (hard limit):  stepsize limited because destination plane is reached";
      break;
    }
    std::cout << "\n";
  }
  std::cout << "\n";
}

} /* End of namespace genfit */
