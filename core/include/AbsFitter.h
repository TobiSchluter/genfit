/* Copyright 2013, Ludwig-Maximilians Universität München,
   Authors: Tobias Schlüter & Johannes Rauch

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

#ifndef genfit_AbsFitter_h
#define genfit_AbsFitter_h

namespace genfit {

class Track;
class AbsTrackRep;
class TVectorD;

class AbsFitter {
 public:
  AbsFitter() {}
  virtual ~AbsFitter() {}

  virtual void processTrack(Track*, AbsTrackRep*, TVectorD&) = 0;
};

}  /* End of namespace genfit */

#endif //genfit_AbsFitter_h
