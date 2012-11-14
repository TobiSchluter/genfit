/* Copyright 2011, Technische Universitaet Muenchen,
   Authors: Karl Bicker, Christian Hoeppner

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

#include "GFDaf.h"
#include "GFTools.h"
#include "GFException.h"
#include "RecoHits/GFDafWireHit.h"
#include "RecoHits/GFAbsWireHit.h"
#include <assert.h>
#include <cmath>

//#define DEBUG


GFDaf::GFDaf(){

  setBetas(81.,8.,4.,1.,1.,1.);
  setProbCut(0.01);
  fKalman.setNumIterations(1);

};

void GFDaf::processTrack(GFTrack* trk) {

#ifdef DEBUG
        std::cout<<"GFDaf::processTrack "<<std::endl;
#endif

  fWeights.clear();

  std::vector<GFDafHit*> eff_hits = initHitsWeights(trk);
  if(eff_hits.size() == 0) {
    for(unsigned int i=0; i<trk->getNumReps(); i++) {
      trk->getTrackRep(i)->setStatusFlag(1);
    }
    return;
  }

  // prepare track
  GFTrack* mini_trk = new GFTrack();
  mini_trk->setSmoothing();

  for(unsigned int j=0; j<eff_hits.size(); j++)
    mini_trk->addHit(eff_hits.at(j), 0, j); // using dummy det and hit id, they are never used anyway


  // fit for each trackrep separately
  for(unsigned int i=0; i<trk->getNumReps(); i++) { // loop over trackreps

    trk->getBK(i)->setNhits(trk->getNumHits());
    if(trk->getTrackRep(i)->getStatusFlag()!=0) continue;

    mini_trk->addTrackRep(trk->getTrackRep(i));

    for(unsigned int iBeta=0; iBeta<fBeta.size(); iBeta++) { // loop over betas

#ifdef DEBUG
        std::cout<<"GFDaf::processTrack, trackRep nr. " << i << ", beta = " << fBeta[iBeta] << std::endl;
#endif

      for(unsigned int j=0; j<mini_trk->getNumHits(); j++) {
        GFDafHit* hit = static_cast<GFDafHit*>(mini_trk->getHit(j));
        hit->setWeights(fWeights.at(i).at(j));
      }
      if ( iBeta != 0){
        fKalman.blowUpCovs(mini_trk);
      }
      fKalman.processTrack(mini_trk);

      if(mini_trk->getTrackRep(0)->getStatusFlag() != 0) break;

      if(iBeta != fBeta.size()-1 )
        try{
          fWeights.at(i) = calcWeights(mini_trk, fBeta.at(iBeta));
        } catch(GFException& e) {
          std::cerr<<e.what();
          e.info();
          mini_trk->getTrackRep(0)->setStatusFlag(1);
          break;
        }

    } // end loop over betas

    if(trk->getSmoothing()) copySmoothing(mini_trk, trk, i);

    mini_trk->releaseTrackReps();

  } // end loop over trackreps

  saveWeights(trk, mini_trk, fWeights);

  delete mini_trk;

};

std::vector<std::vector<double> > GFDaf::calcWeights(GFTrack* trk, double beta) {

  std::vector<std::vector<double> > ret_val;

  for(unsigned int i=0; i<trk->getNumHits(); i++) {

    GFDafHit* eff_hit = static_cast<GFDafHit*>(trk->getHit(i));
    unsigned int nEffHits = eff_hit->getNumEffHits();
    std::vector<double> phi;
    double phi_sum = 0;
    double phi_cut = 0;

    std::vector<double> weights;

    TVectorT<double> smoothedState;
    TMatrixTSym<double> smoothedCov;
    GFTools::getBiasedSmoothedData(trk, 0, i, smoothedState, smoothedCov);

    const TMatrixT<double>& H( trk->getHit(i)->getHMatrix(trk->getTrackRep(0)) );
    TVectorT<double> x_smoo(H * smoothedState);

    if(x_smoo.GetNrows() == 0) {
      weights.assign(nEffHits,0.5);
      //std::cout<<"Assumed weight 0.5!!"<<std::endl;
      ret_val.push_back(weights);
      continue;
    }

    // try to get the plane
    GFDetPlane pl;
    try{
      pl = eff_hit->getDetPlane(trk->getTrackRep(0)); // couldn't figure out how to avoid GFDetPlane copy c'tor here
    } catch(GFException& e) {
      e.what();
      e.info();
      weights.assign(nEffHits,0.);
      ret_val.push_back(weights); // can not get det plane, assign weights 0
      continue; // next eff_hit
    }

    for(unsigned int j=0; j<nEffHits; j++) {
      TVectorT<double> m;
      TMatrixTSym<double> Vorig;
      try{
        eff_hit->getMeasurement(trk->getTrackRep(0),pl,smoothedState,smoothedCov,m,Vorig, j);
      } catch(GFException& e) {
        e.what();
        e.info();
        phi.push_back(0); //m and Vorig do not contain sensible values, assign weight 0
        continue;
      }

      TMatrixTSym<double> V( beta * Vorig);
      TVectorT<double> resid(m - x_smoo);
      double detV = V.Determinant();
      TMatrixTSym<double> Vinv;
      GFTools::invertMatrix(V,Vinv);

      phi.push_back((1./(pow(2.*TMath::Pi(),V.GetNrows()/2.)*sqrt(detV)))*exp(-0.5*Vinv.Similarity(resid)));
      phi_sum += phi.at(j);

      double cutVal = fchi2Cuts[V.GetNrows()];
      assert(cutVal>1.E-6);
      phi_cut += (1./(pow(2.*TMath::Pi(),V.GetNrows()/2.)*sqrt(detV)))*exp(-0.5*cutVal/beta);

    }

    for(unsigned int j=0; j<nEffHits; j++) {

      weights.push_back(phi.at(j)/(phi_sum+phi_cut));

    }
    ret_val.push_back(weights);

  }

  return ret_val;

};

void GFDaf::setProbCut(double val){

  if(fabs(val-0.01)<1.E-10){
    fchi2Cuts[1] = 6.63;
    fchi2Cuts[2] = 9.21;
    fchi2Cuts[3] = 11.34;
    fchi2Cuts[4] = 13.23;
    fchi2Cuts[5] = 15.09;
  }
  else   if(fabs(val-0.005)<1.E-10){
    fchi2Cuts[1] = 7.88;
    fchi2Cuts[2] = 10.60;
    fchi2Cuts[3] = 12.84;
    fchi2Cuts[4] = 14.86;
    fchi2Cuts[5] = 16.75;
  }
  else   if(fabs(val-0.001)<1.E-10){
    fchi2Cuts[1] = 10.83;
    fchi2Cuts[2] = 13.82;
    fchi2Cuts[3] = 16.27;
    fchi2Cuts[4] = 18.47;
    fchi2Cuts[5] = 20.51;
  }
  else{
    GFException exc("GFDafsetProbCut() value is not supported",__LINE__,__FILE__);
    exc.setFatal();
    throw exc;
  }

}

void GFDaf::setBetas(double b1,double b2,double b3,double b4,double b5,double b6,double b7,double b8,double b9,double b10){
  fBeta.clear();
  assert(b1>0);fBeta.push_back(b1);
  if(b2>0){
    assert(b2<=b1);fBeta.push_back(b2);
    if(b3>=0.) {
      assert(b3<=b2);fBeta.push_back(b3);
      if(b4>=0.) {
        assert(b4<=b3);fBeta.push_back(b4);
        if(b5>=0.) {
          assert(b5<=b4);fBeta.push_back(b5);
          if(b6>=0.) {
            assert(b6<=b5);fBeta.push_back(b6);
            if(b7>=0.) {
              assert(b7<=b6);fBeta.push_back(b7);
              if(b8>=0.) {
                assert(b8<=b7);fBeta.push_back(b8);
                if(b9>=0.) {
                  assert(b9<=b8);fBeta.push_back(b9);
                  if(b10>=0.) {
                    assert(b10<=b9);fBeta.push_back(b10);
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

std::vector<GFDafHit*> GFDaf::initHitsWeights(GFTrack* trk) {

  std::vector<GFDafHit*> eff_hits;

  std::vector< std::vector<int>* > planes;
  if(not trk->getHitsByPlane(planes)) return eff_hits;
  int nPlanes = planes.size();

  for(int i=0; i<nPlanes; i++) {

    std::vector<GFAbsRecoHit*> hits;

    for(unsigned int j=0; j<planes.at(i)->size(); j++) {
      hits.push_back(trk->getHit(planes.at(i)->at(j)) );
    }

    GFDafHit* eff_hit;
    if (hits.size()==1 && dynamic_cast<GFAbsWireHit*>(hits[0]) != NULL){
      eff_hit = new GFDafWireHit(dynamic_cast<GFAbsWireHit*>(hits[0]));
    }
    else eff_hit = new GFDafHit(hits);
    eff_hits.push_back(eff_hit);

  }

  for(unsigned int i=0; i<trk->getNumReps(); i++) {
    std::vector<std::vector<double> > rep_weights;
    for(unsigned int j=0; j<eff_hits.size(); j++) {
      rep_weights.push_back(eff_hits.at(j)->getWeights());
    }
    fWeights.push_back(rep_weights);
  }

  return eff_hits;

}

void GFDaf::copySmoothing(const GFTrack* source, GFTrack* target, int target_irep) {

  fKalman.initBookkeeping(target);

  unsigned int hit_count = 0;
  GFBookkeeping* sourceBK0 = source->getBK(0);
  GFBookkeeping* targetBKirep = target->getBK(target_irep);

  std::vector< std::string > keys;
  TVectorT<double>& vec = *(new TVectorT<double>);
  TMatrixT<double>& mat = *(new TMatrixT<double>);
  TMatrixTSym<double>& matSym = *(new TMatrixTSym<double>);
  GFDetPlane& pl = *(new GFDetPlane);
  double numVal = 0;
  double& num = numVal;


  for(unsigned int pl_i=0; pl_i<source->getNumHits(); pl_i++) {

    GFDafHit* eff_hit = static_cast<GFDafHit*>(source->getHit(pl_i));

    for(unsigned int hit_i=0; hit_i<eff_hit->getNumRealHits(); ++hit_i) {

      // copy vectors
      keys = sourceBK0->getVectorKeys();
      for (unsigned int iKey=0; iKey<keys.size(); ++iKey){
        sourceBK0->getVector(keys[iKey], pl_i, vec);
        targetBKirep->setVector(keys[iKey], hit_count, vec);
      }

      // copy matrices
      keys = sourceBK0->getMatrixKeys();
      for (unsigned int iKey=0; iKey<keys.size(); ++iKey){
        sourceBK0->getMatrix(keys[iKey], pl_i, mat);
        targetBKirep->setMatrix(keys[iKey], hit_count, mat);
      }

      // copy symmetric matrices
      keys = sourceBK0->getSymMatrixKeys();
      for (unsigned int iKey=0; iKey<keys.size(); ++iKey){
        sourceBK0->getSymMatrix(keys[iKey], pl_i, matSym);
        targetBKirep->setSymMatrix(keys[iKey], hit_count, matSym);
      }

      // copy det planes
      keys = sourceBK0->getGFDetPlaneKeys();
      for (unsigned int iKey=0; iKey<keys.size(); ++iKey){
        sourceBK0->getDetPlane(keys[iKey], pl_i, pl);
        targetBKirep->setDetPlane(keys[iKey], hit_count, pl);
      }

      // copy numbers
      keys = sourceBK0->getNumberKeys();
      for (unsigned int iKey=0; iKey<keys.size(); ++iKey){
        sourceBK0->getNumber(keys[iKey], pl_i, num);
        targetBKirep->setNumber(keys[iKey], hit_count, num);
      }

      ++hit_count;
    }
  }

  assert(hit_count == target->getNumHits());
}

void GFDaf::saveWeights(GFTrack* trk, const GFTrack* DafTrack, const std::vector<std::vector<std::vector<double> > >& weights) const {

  assert(weights.size() == trk->getNumReps());

  // get number of real and effective hits
  unsigned int nDafHits(DafTrack->getNumHits());
  std::vector<unsigned int> nRealHits, nEffHits;
  nRealHits.reserve(nDafHits);
  nEffHits.reserve(nDafHits);
  GFDafHit* hit;
  for (unsigned int i=0; i<nDafHits; ++i){
    hit = static_cast<GFDafHit*>(DafTrack->getHit(i));
    nRealHits[i] = hit->getNumRealHits();
    nEffHits[i] = hit->getNumEffHits();
  }

  TVectorT<double> vec;

  for(unsigned int rep_i = 0; rep_i < weights.size(); rep_i++) { // loop over trackReps
    GFBookkeeping* bk = trk->getBK(rep_i);
    bk->bookVectors("dafWeight");
    unsigned int hit_i = 0;


    for (unsigned int i=0; i<nDafHits; ++i){ // loop over daf hits
      if (nRealHits[i] == nEffHits[i]) {
        for (unsigned int j=0; j<nRealHits[i]; ++j){
          vec.ResizeTo(1);
          vec(0) = weights.at(rep_i).at(i).at(j);
        }
      }
      else {
        assert (nRealHits[i] == 1);
        vec.ResizeTo(nEffHits[i]);
        for (unsigned int j=0; j<nEffHits[i]; ++j){
          vec[j] = weights.at(rep_i).at(i).at(j);
        }
      }
      bk->setVector("dafWeight", hit_i++, vec);
    }

    assert(hit_i == trk->getNumHits());
  }

}

