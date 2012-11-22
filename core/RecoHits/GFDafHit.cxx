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


#include "GFDafHit.h"
#include <GFTools.h>
#include <GFException.h>
#include <cmath>
#include <memory>


GFDafHit::GFDafHit(std::vector<GFAbsRecoHit*> HitsInPlane) {

	fRawHits = HitsInPlane;
	fWeights.assign(fRawHits.size(),1.);
	fHitUpd = false;

}

GFAbsRecoHit* GFDafHit::getHit(unsigned int ihit) {

	return fRawHits.at(ihit);

}

void GFDafHit::setWeights(std::vector<double> weights) {

	fWeights = weights;
	static const double minWeight = 1.e-10;
	for (unsigned int i=0; i<fWeights.size(); ++i) {
	  if (fWeights[i] < minWeight) fWeights[i] = minWeight;
	}
	fHitUpd = false;
}

const GFDetPlane& GFDafHit::getDetPlane(GFAbsTrackRep* rep) {

	return fRawHits[0]->getDetPlane(rep);

}

void GFDafHit::getMeasurement(const GFAbsTrackRep* rep,const GFDetPlane& pl,const TVectorD& statePred,const TMatrixDSym& covPred,TVectorD& m, TMatrixDSym& V) {

  if(fHitUpd && fDetPlane == pl) {
    m.ResizeTo(fHitCoord);
    V.ResizeTo(fHitCov);
    m = fHitCoord;
    V = fHitCov;
    return;
  }

  if(fRawHits.size() == 1) {
    fRawHits[0]->getMeasurement(rep, pl, statePred, covPred, fHitCoord, fHitCov);
    fHitCov = (1. / fWeights[0]) * fHitCov;
  } 
  else { // more than one hit
    TMatrixDSym covInv;
    TVectorD coordTemp;
    TMatrixDSym covTemp;
    std::vector<TVectorD> coords;
    std::vector<TMatrixDSym> covInvs;

    for(unsigned int i=0; i<fRawHits.size(); ++i) {
      fRawHits[i]->getMeasurement(rep, pl, statePred, covPred, coordTemp, covTemp);

      // make sure fHitCoord and fHitCov have right dimensionality and set them to 0
      if (i==0){
        fHitCoord.ResizeTo(coordTemp);
        fHitCoord.Zero();
        fHitCov.ResizeTo(covTemp);
        fHitCov.Zero();
      }

      coords.push_back(coordTemp);
      GFTools::invertMatrix(covTemp, covInv);
      covInvs.push_back(covInv);
      fHitCov += fWeights[i] * covInv;
    }

    // invert fHitCov
    TMatrixDSym HitCovTemp(fHitCov);
    GFTools::invertMatrix(HitCovTemp, fHitCov);

    //set the weighted-mean coord
    for(unsigned int i=0; i<coords.size(); ++i) {
      fHitCoord += fWeights[i] * covInvs[i] * coords[i];
    }
    fHitCoord = fHitCov * fHitCoord;
  }

  //return by refernce
  m.ResizeTo(fHitCoord);
  V.ResizeTo(fHitCov);
  m = fHitCoord;
  V = fHitCov;
  fDetPlane = pl;
  fHitUpd = true;
}


const TMatrixD& GFDafHit::getHMatrix(const GFAbsTrackRep* rep) {

	return fRawHits[0]->getHMatrix(rep);

}

GFDafHit* GFDafHit::clone() {

	GFDafHit* retval = new GFDafHit(fRawHits);
	retval->setWeights(fWeights);
	return retval;

}


ClassImp(GFDafHit)
