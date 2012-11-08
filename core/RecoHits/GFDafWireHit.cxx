/* Copyright 2011, Technische Universitaet Muenchen,
Authors: Karl Bicker

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


#include "GFDafWireHit.h"
#include <GFTools.h>
#include <GFException.h>
#include <cmath>


GFDafWireHit::GFDafWireHit(std::vector<GFAbsWireHit*> Hit)
  : GFDafHit()
{

  assert(Hit.size() == 1);

  fRawHit = Hit[0];
  fWeights.assign(2,1.);
  fHitUpd = false;

	// set initial weights according to l/r resolution (default weights are left:1 right:1)
	int lr = fRawHit->getLeftRightResolution();

	if (lr<0) { // left
	  fWeights[1] = 0.; // set right to 0
	}
	else if (lr>0){ // right
    fWeights[0] = 0.; // set left to 0
	}

	// set l/r resolution so that the plane is fixed
	fRawHit->setLeftRightResolution(1);
}


GFDafWireHit::~GFDafWireHit(){
  // set the l/r resolution of the wire hit according to what the DAF has determined
  if (fWeights[0] > fWeights[1]) fRawHit->setLeftRightResolution(-1);
  else fRawHit->setLeftRightResolution(1);
}


void GFDafWireHit::getMeasurement(const GFAbsTrackRep* rep,const GFDetPlane& pl,const TVectorT<double>& statePred,const TMatrixTSym<double>& covPred,TVectorT<double>& m, TMatrixTSym<double>& V) {

  if(fHitUpd && fDetPlane == pl) {
    m.ResizeTo(fHitCoord);
    V.ResizeTo(fHitCov);
    m = fHitCoord;
    V = fHitCov;
    return;
  }

  TMatrixTSym<double> covInv;
  TVectorT<double> coordTemp;
  TMatrixTSym<double> covTemp;
  std::vector< TVectorT<double> > coords;

  fRawHit->getMeasurement(rep, pl, statePred, covPred, coordTemp, covTemp);

  // make sure fHitCoord and fHitCov have right dimensionality and set them to 0
  fHitCoord.ResizeTo(coordTemp);
  fHitCoord.Zero();
  fHitCov.ResizeTo(covTemp);
  fHitCov.Zero();

  GFTools::invertMatrix(covTemp, covInv);

  for(unsigned int i=0; i<2; ++i) {
    coords.push_back(coordTemp);
    fHitCov += fWeights[i] * covInv;
  }
  coords[0](0) *= -1.; // invert the sign of the drift radius of the first (left) hit

  // invert fHitCov
  TMatrixTSym<double> HitCovTemp(fHitCov);
  GFTools::invertMatrix(HitCovTemp, fHitCov);

  //set the weighted-mean coord
  for(unsigned int i=0; i<2; ++i) {
    fHitCoord += fWeights[i] * covInv * coords[i];
  }
  fHitCoord = fHitCov * fHitCoord;

  //return by refernce
  m.ResizeTo(fHitCoord);
  V.ResizeTo(fHitCov);
  m = fHitCoord;
  V = fHitCov;
  fDetPlane = pl;
  fHitUpd = true;
}


ClassImp(GFDafWireHit)
