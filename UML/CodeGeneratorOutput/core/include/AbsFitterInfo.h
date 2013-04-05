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

#ifndef genfit_AbsFitterInfo_h
#define genfit_AbsFitterInfo_h

#include "TObject.h"

namespace genfit {


  /** 
   *  This class collects all information needed and produced by a specific fitter and is specific to one #GFAbsTrackRep of the #GFTrack.
   */
class AbsFitterInfo : public TObject {

 public:

  AbsFitterInfo() {};
  virtual ~AbsFitterInfo();

  //! Deep copy ctor for polymorphic class.
  virtual AbsFitterInfo* clone() const = 0;

 private:

  AbsFitterInfo(const AbsFitterInfo&) = default;
  AbsFitterInfo(AbsFitterInfo&&) = default;
  AbsFitterInfo& operator=(const AbsFitterInfo&) = default;
  AbsFitterInfo& operator=(AbsFitterInfo&&) = default;


  ClassDef(AbsFitterInfo,1)

};

} /* End of namespace genfit */

#endif // genfit_AbsFitterInfo_h
