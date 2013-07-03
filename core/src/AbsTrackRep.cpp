/* Copyright 2008-2009, Technische Universitaet Muenchen,
   Authors: Christian Hoeppner & Sebastian Neubert

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

#include "AbsTrackRep.h"

#include <iostream>


namespace genfit {

AbsTrackRep::AbsTrackRep() :
  pdgCode_(0), propDir_(0)
{
  ;
}

AbsTrackRep::AbsTrackRep(int pdgCode, char propDir) :
  pdgCode_(pdgCode), propDir_(propDir)
{
  ;
}

AbsTrackRep::AbsTrackRep(const AbsTrackRep& rep) :
  TObject(rep), pdgCode_(rep.pdgCode_), propDir_(rep.propDir_)
{
  ;
}


TVectorD AbsTrackRep::get6DState(const StateOnPlane* stateInput) const {
  TVector3 pos, mom;
  getPosMom(stateInput, pos, mom);

  TVectorD stateVec(6);

  stateVec(0) = pos.X();
  stateVec(1) = pos.Y();
  stateVec(2) = pos.Z();

  stateVec(3) = mom.X();
  stateVec(4) = mom.Y();
  stateVec(5) = mom.Z();

  return stateVec;
}


void AbsTrackRep::get6DStateCov(const MeasuredStateOnPlane* stateInput, TVectorD& stateVec, TMatrixDSym& cov) const {
  TVector3 pos, mom;
  getPosMomCov(stateInput, pos, mom, cov);

  stateVec.ResizeTo(6);

  stateVec(0) = pos.X();
  stateVec(1) = pos.Y();
  stateVec(2) = pos.Z();

  stateVec(3) = mom.X();
  stateVec(4) = mom.Y();
  stateVec(5) = mom.Z();
}


void AbsTrackRep::Print(const Option_t*) const {
  std::cout << "genfit::TrackRep, pdgCode = " << pdgCode_ << ". PropDir = " << (int)propDir_ << "\n";
}


} /* End of namespace genfit */
