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
 *  @author Christian H&ouml;ppner (Technische Universit&auml;t M&uuml;nchen, original author)
 *
 */


/** @addtogroup genfit
 * @{
 */

#ifndef genfit_AbsFinitePlane_h
#define genfit_AbsFinitePlane_h


#include <TObject.h>


namespace genfit {

class AbsFinitePlane : public TObject {

 public:

  AbsFinitePlane() {};
  virtual ~AbsFinitePlane() {};

  //! Returns whether a u,v point is in the active plane or not. Needs to be implemented
  //! in child class.
  virtual bool isInActive(double u, double v) const = 0;

  //! Deep copy ctor for polymorphic class.
  virtual AbsFinitePlane* clone() const = 0;


 protected:

  // protect from calling copy c'tor or assignment operator from outside the class. Use #clone() if you want a copy!
  AbsFinitePlane(const AbsFinitePlane& o) : TObject(o) {;}
  AbsFinitePlane& operator=(const AbsFinitePlane&);


  ClassDef(AbsFinitePlane,1)
};

} /* End of namespace genfit */
/** @} */

#endif // genfit_AbsFinitePlane_h
