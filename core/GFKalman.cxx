/* Copyright 2008-2010, Technische Universitaet Muenchen,
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
#include "GFKalman.h"

#include "assert.h"
#include <iostream>
#include <sstream>

#include "TMath.h"

#include "GFTrack.h"
#include "GFAbsRecoHit.h"
#include "GFAbsTrackRep.h"
#include "GFBookkeeping.h"
#include "GFException.h"
#include "GFTools.h"

#define COVEXC "cov_is_zero"
//#define DEBUG

GFKalman::GFKalman():fInitialDirection(1),fNumIt(3),fBlowUpFactor(500.){;}

GFKalman::~GFKalman(){;}

void GFKalman::processTrack(GFTrack* trk){
#ifdef DEBUG
        std::cout<<"GFKalman::processTrack "<<std::endl;
#endif

  fSmooth = trk->getSmoothing();
  fSmoothFast = trk->getSmoothingFast();

  int nreps = trk->getNumReps();
  for(int i=0; i<nreps; i++) {
    GFBookkeeping* bk = trk->getBK(i);
    bk->setNhits(trk->getNumHits());
    if(fSmooth) {
      std::vector<std::string> vec_keys = bk->getVectorKeys();
      bool already_there = false;
      for(unsigned int j=0; j<vec_keys.size(); j++) {
        if(vec_keys.at(j) == "fUpSt") {
          already_there = true;
          break;
        }
      }
      if(!already_there) {
        bk->bookNumbers("fExtLen"); // extrapolated length from last hit in forward direction
        bk->bookVectors("fUpSt");
        bk->bookSymMatrices("fUpCov");
        bk->bookNumbers("bExtLen"); // extrapolated length from last hit in backward direction
        bk->bookVectors("bUpSt");
        bk->bookSymMatrices("bUpCov");
        if(fSmoothFast) {
          bk->bookVectors("fSt");
          bk->bookSymMatrices("fCov");
          bk->bookVectors("bSt");
          bk->bookSymMatrices("bCov");
        }
        bk->bookGFDetPlanes("fPl");
        bk->bookGFDetPlanes("bPl");
        if(trk->getTrackRep(i)->hasAuxInfo()) {
          bk->bookMatrices("fAuxInfo");
          bk->bookMatrices("bAuxInfo");
        }
      }
    }
  }

  int direction=fInitialDirection;
  assert(direction==1 || direction==-1);
  //  trk->clearGFBookkeeping();
  trk->clearRepAtHit();

  /*why is there a factor of two here (in the for statement)?
    Because we consider one full iteration to be one back and
    one forth fitting pass */
  for(int ipass=0; ipass<2*fNumIt; ipass++){
    if(ipass>0) blowUpCovs(trk);

    // reset X/X0 before last fitting pass
    if(ipass==(2*fNumIt)-1) {
      for(int i=0; i<nreps; ++i) {
        trk->getTrackRep(i)->resetXX0();
      }
    }

    if(direction==1){
      trk->setNextHitToFit(0);
    }
    else {
      trk->setNextHitToFit(trk->getNumHits()-1);
    }
    fittingPass(trk,direction);
    
    //save first and last plane,state&cov after the fitting pass
    if(direction==1){//forward at last hit
      for(int i=0; i<nreps; ++i){
        GFAbsTrackRep* rep = trk->getTrackRep(i);
        rep->setLastPlane( rep->getReferencePlane() );
        rep->setLastState( rep->getState() );
        rep->setLastCov( rep->getCov() );
      }
    }
    else{//backward at first hit
      for(int i=0; i<nreps; ++i){
        GFAbsTrackRep* rep = trk->getTrackRep(i);
        rep->setFirstPlane( rep->getReferencePlane() );
        rep->setFirstState( rep->getState() );
        rep->setFirstCov( rep->getCov() );
      }
    }

    //switch direction of fitting and also inside all the reps
    if(direction==1) direction=-1;
    else direction=1;
    switchDirection(trk);
  }
  
  return;
}

void
GFKalman::switchDirection(GFTrack* trk){
  int nreps=trk->getNumReps();
  for(int i=0; i<nreps; ++i){
    trk->getTrackRep(i)->switchDirection();
  }
}

void GFKalman::blowUpCovs(GFTrack* trk){
  trk->blowUpCovs(fBlowUpFactor);
}

void
GFKalman::fittingPass(GFTrack* trk, int direction){
  //loop over hits
  unsigned int nhits=trk->getNumHits();
  unsigned int starthit=trk->getNextHitToFit();

  int nreps=trk->getNumReps();
  int ihit=(int)starthit;
  
  for(int irep=0; irep<nreps; ++irep) {
    GFAbsTrackRep* arep=trk->getTrackRep(irep);
    if(arep->getStatusFlag()==0) {
      //clear chi2 sum and ndf sum in track reps
        if (direction == -1){
          arep->setChiSqu(0.);
        }
        if (direction == 1){
          arep->setForwardChiSqu(0.);
        }
      arep->setNDF(0);
      //clear failedHits and outliers
      trk->getBK(irep)->clearFailedHits();
    }
  }

  while((ihit<(int)nhits && direction==1) || (ihit>-1 && direction==-1)){
    //    GFAbsRecoHit* ahit=trk->getHit(ihit);
    // loop over reps
    for(int irep=0; irep<nreps; ++irep){
    GFAbsTrackRep* arep=trk->getTrackRep(irep);
    if(arep->getStatusFlag()==0) {
      try {
#ifdef DEBUG
        std::cout<<"++++++++++++++++++++++++++++++++++++++++\n";
        std::cout<<"GFKalman::fittingPass - process rep nr. "<<irep<<" and hit nr. "<<ihit<<std::endl;
#endif
        processHit(trk,ihit,irep,direction);
      }
      catch(GFException& e) {
        trk->addFailedHit(irep,ihit);
        std::cerr << e.what();
        e.info();
        if(e.isFatal()) {
          arep->setStatusFlag(1);
          continue; // go to next rep immediately
        }
      }
    }
    }// end loop over reps
    ihit+=direction;
  }// end loop over hits
  trk->setNextHitToFit(ihit-2*direction);
  //trk->printGFBookkeeping();
}

double GFKalman::chi2Increment(const TVectorT<double>& r,const TMatrixT<double>& H,
           const TMatrixTSym<double>& cov,const TMatrixTSym<double>& V){

  // residuals covariances:R=(V - HCH^T)
  TMatrixTSym<double> HcovHt(cov);
  HcovHt.Similarity(H);

  // instead of
  //  TMatrixTSym<double> Rinv = TMatrixTSym(V, kMinus, HcovHt)
  // because kMinus constructor doesn't work in root up to at least 5.34.
  // Bug report: <https://savannah.cern.ch/bugs/index.php?98605>

  // chisq= r^TR^(-1)r
  TMatrixTSym<double> Rinv(V - HcovHt);
  GFTools::invertMatrix(Rinv);
  double chisq = Rinv.Similarity(r);

  if(TMath::IsNaN(chisq)){
    GFException exc("chi2 is nan",__LINE__,__FILE__);
    exc.setFatal();
    std::vector< TMatrixT<double> > matrices;
    matrices.push_back(V);
    matrices.push_back(Rinv);
    matrices.push_back(cov);
    exc.setMatrices("V, R, cov",matrices);
    throw exc;
  }

  return chisq;
}


double
GFKalman::getChi2Hit(GFAbsRecoHit* hit, GFAbsTrackRep* rep)
{
  // get prototypes for matrices
  int repDim=rep->getDim();
  TVectorT<double> state(repDim);
  TMatrixTSym<double> cov(repDim);;
  GFDetPlane pl=hit->getDetPlane(rep);
  rep->extrapolate(pl,state,cov);


  const TMatrixT<double>& H = hit->getHMatrix(rep);
  TVectorT<double> m;
  TMatrixTSym<double> V;
  hit->getMeasurement(rep,pl,state,cov,m,V);

  TVectorT<double> res = m-(H*state);
  assert(res.GetNrows()>0);

  //this is where chi2 is calculated
  double chi2 = chi2Increment(res,H,cov,V);

  return chi2/res.GetNrows();
}

  
void
GFKalman::processHit(GFTrack* tr, int ihit, int irep,int direction){
  GFAbsRecoHit* hit = tr->getHit(ihit);
  GFAbsTrackRep* rep = tr->getTrackRep(irep);

  // get prototypes for matrices
  int repDim = rep->getDim();
  TVectorT<double> state(repDim);
  TMatrixTSym<double> cov(repDim);
  const GFDetPlane* ppl;

  double extLen(0.);

  /* do an extrapolation, if the trackrep irep is not given
   * at this ihit position. This will usually be the case, but
   * not if the fit turnes around
   */
  if(ihit!=tr->getRepAtHit(irep)){
    // get the (virtual) detector plane
    ppl=&hit->getDetPlane(rep);
    //do the extrapolation
    extLen = rep->extrapolate(*ppl,state,cov);
  }
  else{
    ppl = &rep->getReferencePlane();
    state = rep->getState();
    cov = rep->getCov();
    extLen = 0.;
  }
  const GFDetPlane& pl = *ppl;

  if(cov(0,0)<1.E-50){ // diagonal elements must be >=0
    GFException exc(COVEXC,__LINE__,__FILE__);
    throw exc;
  }

  GFBookkeeping* bk = tr->getBK(irep);
  if(fSmooth && fSmoothFast) {
    if(direction == 1) {
	    bk->setVector("fSt",ihit,state);
	    bk->setSymMatrix("fCov",ihit,cov);
	    if(rep->hasAuxInfo()) bk->setMatrix("fAuxInfo",ihit,*(rep->getAuxInfo(pl)));
	    bk->setDetPlane("fPl",ihit,pl);
	  } else {
	    bk->setVector("bSt",ihit,state);
	    bk->setSymMatrix("bCov",ihit,cov);
	    if(rep->hasAuxInfo()) bk->setMatrix("bAuxInfo",ihit,*(rep->getAuxInfo(pl)));
	    bk->setDetPlane("bPl",ihit,pl);
	  }
  }
  
#ifdef DEBUG
  std::cerr<<"GFKalman::processHit - state and cov prediction "<<std::endl;
  state.Print();
  cov.Print();
#endif

  const TMatrixT<double>& H(hit->getHMatrix(rep));
  TVectorT<double> m;
  TMatrixTSym<double> V;
  hit->getMeasurement(rep,pl,state,cov,m,V);
  TVectorT<double> res = m-(H*state);

  // calculate kalman gain ------------------------------
  // calculate covsum (V + HCH^T)
  TMatrixTSym<double> HcovHt(cov);
  HcovHt.Similarity(H);
  
  // invert
  //TMatrixTSym<double> covSum(V + HcovHt);
  // instead of:
  //  TMatrixTSym<double> covSum(V, TMatrixTSym<double>::kPlus, HcovHt);
  // because kPlus constructor doesn't work in root up to at least 5.34.
  // Bug report: <https://savannah.cern.ch/bugs/index.php?98605>
  TMatrixTSym<double> covSumInv(V + HcovHt);
  GFTools::invertMatrix(covSumInv);

  TMatrixT<double> CHt(cov, TMatrixT<double>::kMultTranspose, H);
  TVectorT<double> update = TMatrixT<double>(CHt, TMatrixT<double>::kMult, covSumInv) * res;
#ifdef DEBUG
  std::cout<<"residual vector"; res.Print();
  std::cout<<"update = Gain*res"; update.Print();
#endif

  state+=update; // prediction overwritten!

  // And the new covariance matrix:
  covSumInv.Similarity(CHt);
  cov-=covSumInv;  // Cnew = C - C Ht (V + H C Ht)^-1 H C

  if(fSmooth) {
    if(direction == 1) {
      bk->setNumber("fExtLen",ihit,extLen);
      bk->setVector("fUpSt",ihit,state);
      bk->setSymMatrix("fUpCov",ihit,cov);
      if(rep->hasAuxInfo()) bk->setMatrix("fAuxInfo",ihit,*(rep->getAuxInfo(pl)));
      bk->setDetPlane("fPl",ihit,pl);
    } else {
	    bk->setNumber("bExtLen",ihit,extLen);
      bk->setVector("bUpSt",ihit,state);
      bk->setSymMatrix("bUpCov",ihit,cov);
      if(rep->hasAuxInfo()) bk->setMatrix("bAuxInfo",ihit,*(rep->getAuxInfo(pl)));
      bk->setDetPlane("bPl",ihit,pl);
    }
  }

  // calculate filtered chisq
  // filtered residual
  res = m-(H*state);
  double chi2 = chi2Increment(res,H,cov,V);
  int ndf = res.GetNrows();
  if (direction == -1) {
    rep->addChiSqu( chi2 );
  }
  if (direction == 1) {
    rep->addForwardChiSqu( chi2 );
  }
  rep->addNDF( ndf );

  /*
  if(direction==1){
    bk->setNumber("fChi2",ihit,chi2/ndf);
  }
  else{
    bk->setNumber("bChi2",ihit,chi2/ndf);
  }
  */

  // if we survive until here: update TrackRep
  //rep->setState(state);
  //rep->setCov(cov);
  //rep->setReferencePlane(pl);

  rep->setData(state,pl,&cov);
  tr->setRepAtHit(irep,ihit);

#ifdef DEBUG
   std::cout<<"GFKalman::processHit - updated state and cov "<<std::endl;
   rep->getState().Print();
   rep->getCov().Print();
#endif
}







