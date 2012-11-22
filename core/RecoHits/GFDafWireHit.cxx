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


GFDafWireHit::GFDafWireHit(GFAbsWireHit* hit)
  : GFDafHit()
{
  // fill fRawHits with the hit
  fRawHits.push_back(hit);

  fHitUpd = false;

  // initialise the two weigths for left and right hit
	// set initial weights according to l/r resolution (default weights are left:1 right:1)
  fWeights.assign(2,1.);
	int lr = dynamic_cast<GFAbsWireHit*>(fRawHits[0])->getLeftRightResolution();

	if (lr<0) { // left
	  fWeights[1] = 0.; // set right to 0
	}
	else if (lr>0){ // right
    fWeights[0] = 0.; // set left to 0
	}

	// set l/r resolution so that the plane is fixed
	dynamic_cast<GFAbsWireHit*>(fRawHits[0])->setLeftRightResolution(1);
}


GFDafWireHit::~GFDafWireHit(){
  // set the l/r resolution of the wire hit according to what the DAF has determined
  if (fWeights[0] > fWeights[1]) dynamic_cast<GFAbsWireHit*>(fRawHits[0])->setLeftRightResolution(-1);
  else dynamic_cast<GFAbsWireHit*>(fRawHits[0])->setLeftRightResolution(1);
}


void GFDafWireHit::getMeasurement(const GFAbsTrackRep* rep,const GFDetPlane& pl,const TVectorD& statePred,const TMatrixDSym& covPred,TVectorD& m, TMatrixDSym& V) {

  if(fHitUpd && fDetPlane == pl) {
    m.ResizeTo(fHitCoord);
    V.ResizeTo(fHitCov);
    m = fHitCoord;
    V = fHitCov;
    return;
  }

  TMatrixDSym covInv;
  TVectorD coordTemp;
  TMatrixDSym covTemp;
  std::vector<TVectorD> coords;

  dynamic_cast<GFAbsWireHit*>(fRawHits[0])->getMeasurement(rep, pl, statePred, covPred, coordTemp, covTemp);

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
  TMatrixDSym HitCovTemp(fHitCov);
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


void GFDafWireHit::getMeasurement(const GFAbsTrackRep* rep,const GFDetPlane& pl,const TVectorD& statePred,const TMatrixDSym& covPred,TVectorD& m, TMatrixDSym& V, unsigned int iHit){
  assert(iHit<2);

  fRawHits[0]->getMeasurement(rep, pl, statePred, covPred, m, V);

  if (iHit == 0){
    m[0] *= -1.; // invert the sign of the drift radius of the left hit
  }
}

ClassImp(GFDafWireHit)
