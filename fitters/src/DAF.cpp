/* Copyright 2008-2013, Technische Universitaet Muenchen, Ludwig-Maximilians-Universität München
   Authors: Christian Hoeppner & Sebastian Neubert & Johannes Rauch & Tobias Schlüter

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

#include "DAF.h"
#include "Exception.h"
#include "KalmanFitterInfo.h"
#include "KalmanFitterRefTrack.h"
#include "KalmanFitStatus.h"
#include "Tools.h"
#include "Track.h"
#include "TrackPoint.h"

#include <assert.h>
#include <cmath>

//root stuff
#include <TMath.h>
#include <Math/QuantFuncMathCore.h>


//#define DEBUG


namespace genfit {

DAF::DAF()
  : AbsKalmanFitter(c_maxIter), deltaWeight_(0.001)
{
  kalman_.reset(new KalmanFitterRefTrack());
  kalman_->setMultipleMeasurementHandling(weightedAverage);
  kalman_->setNumIterations(1);

  setAnnealingScheme(100, 0.1, 5);
  setProbCut(0.01);
}

DAF::DAF(AbsKalmanFitter* kalman)
  : AbsKalmanFitter(c_maxIter), deltaWeight_(0.001)
{
  kalman_.reset(kalman);
  kalman_->setNumIterations(1);

  setAnnealingScheme(100, 0.1, 5);
  setProbCut(0.01);
}


void DAF::processTrack(Track* tr, const AbsTrackRep* rep) {

#ifdef DEBUG
  std::cout<<"DAF::processTrack \n";
#endif

  weights_.clear();
  std::vector<std::vector<double> > oldWeights;

  KalmanFitStatus* status;
  bool oneLastIter = false;

  unsigned int iBeta = 0;
  for(; iBeta != c_maxIter; ++iBeta) { // loop over. If no convergence is reached after 10 iterations just stop.

#ifdef DEBUG
      std::cout<<"DAF::processTrack, trackRep  " << rep << ", beta = " << betas_[iBeta] << "\n";
#endif

    kalman_->processTrack(tr);

    status = static_cast<KalmanFitStatus*>(tr->getFitStatus(rep));
    status->setIsFittedWithDaf();

    if (! status->isFitted()){
      #ifdef DEBUG
      std::cout << "DAF::Kalman could not fit!\n";
      #endif
      status->setIsFitted(false);
      break;
    }

    if (iBeta > 0)
      oldWeights = weights_.at(rep);
    getWeights(tr, rep);

    try{
      weights_.at(rep) = calcWeights(tr, rep, betas_[iBeta]);
    } catch(Exception& e) {
      std::cerr<<e.what();
      e.info();
      //std::cerr << "calc weights failed" << std::endl;
      //mini_trk->getTrackRep(0)->setStatusFlag(1);
      status->setIsFitted(false);
      status->setIsFitConverged(false);
      break;
    }
    if( oneLastIter == true){
      #ifdef DEBUG
      std::cout << "DAF::break after one last iteration\n";
      #endif
      status->setIsFitConverged();
      break;
    }
    if (iBeta > 0) {
      if ( isConvergent(oldWeights, rep) ){
        #ifdef DEBUG
        std::cout << "DAF::convergence reached in iteration " << iBeta << "\n";
        #endif
        oneLastIter = true;
      }
    }

      if(iBeta == maxIterations_-1 ){
        status->setIsFitConverged(false);
        #ifdef DEBUG
        std::cout << "DAF::number of max iterations reached!\n";
        #endif
        break;
      }
  } // end loop over betas

  status->setNumIterations(iBeta);

}


void DAF::setProbCut(const double prob_cut){
  for ( int i = 1; i != 6; ++i){
    addProbCut(prob_cut, i);
  }
}

void DAF::addProbCut(const double prob_cut, const int measDim){
  if ( prob_cut > 1.0 || prob_cut < 0.0){
    Exception exc("DAF::addProbCut prob_cut is not between 0 and 1",__LINE__,__FILE__);
    exc.setFatal();
    throw exc;
  }
  if ( measDim < 1){
    Exception exc("DAF::addProbCut measDim must be > 0",__LINE__,__FILE__);
    exc.setFatal();
    throw exc;
  }
  chi2Cuts_[measDim] = ROOT::Math::chisquared_quantile_c( prob_cut, measDim);
}


void DAF::setBetas(double b1,double b2,double b3,double b4,double b5,double b6,double b7,double b8,double b9,double b10){
  betas_.clear();
  assert(b1>0);betas_.push_back(b1);
  if(b2>0){
    assert(b2<=b1);betas_.push_back(b2);
    if(b3>=0.) {
      assert(b3<=b2);betas_.push_back(b3);
      if(b4>=0.) {
        assert(b4<=b3);betas_.push_back(b4);
        if(b5>=0.) {
          assert(b5<=b4);betas_.push_back(b5);
          if(b6>=0.) {
            assert(b6<=b5);betas_.push_back(b6);
            if(b7>=0.) {
              assert(b7<=b6);betas_.push_back(b7);
              if(b8>=0.) {
                assert(b8<=b7);betas_.push_back(b8);
                if(b9>=0.) {
                  assert(b9<=b8);betas_.push_back(b9);
                  if(b10>=0.) {
                    assert(b10<=b9);betas_.push_back(b10);
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  betas_.resize(c_maxIter,betas_.back()); //make sure main loop has a maximum of 10 iterations and also make sure the last beta value is used for if more iterations are needed then the ones set by the user.
}


void DAF::setAnnealingScheme(double bStart, double bFinal, unsigned int nSteps) {
  assert(bStart > bFinal);
  assert(bFinal > 1.E-10);
  assert(nSteps > 1);
  assert(nSteps <= c_maxIter);

  betas_.clear();

  const double x = log(bStart/bFinal)/log(2)/(nSteps-1);

  for (unsigned int i=0; i<nSteps; ++i) {
    double exp = double(nSteps-i-1)*x;
    betas_.push_back(bFinal * pow(2., exp));
  }

  betas_.resize(c_maxIter,betas_.back()); //make sure main loop has a maximum of 10 iterations and also make sure the last beta value is used for if more iterations are needed then the ones set by the user.

  /*for (unsigned int i=0; i<betas_.size(); ++i) {
    std::cout<< betas_[i] << ", ";
  }*/
}


bool DAF::isConvergent(const std::vector<std::vector<double> >& oldWeights, const AbsTrackRep* rep) const {
  const int n = oldWeights.size();
  const std::vector<std::vector<double> >& newWeights = weights_.at(rep);
  assert(n == int(newWeights.size()));
  for( int i = 0; i != n; ++i){
    const int m = oldWeights[i].size();
    assert(m == int(newWeights[i].size()));
    for( int j = 0; j != m; ++j){
      if( fabs(oldWeights[i][j] - newWeights[i][j]) > deltaWeight_ ){
        return false;
      }
    }
  }
  return true;
}


void DAF::getWeights(const Track* tr, const AbsTrackRep* rep) {

#ifdef DEBUG
      std::cout<<"DAF::getWeights \n";
#endif

  std::vector< TrackPoint* > trackPoints = tr->getPointsWithMeasurement();
  for (std::vector< TrackPoint* >::iterator tp = trackPoints.begin(); tp != trackPoints.end(); ++tp) {
    AbsFitterInfo* fi = (*tp)->getFitterInfo(rep);
    if (dynamic_cast<KalmanFitterInfo*>(fi) == NULL){
      Exception exc("DAF::getWeights ==> can only use KalmanFitterInfos",__LINE__,__FILE__);
      throw exc;
    }
    KalmanFitterInfo* kfi = static_cast<KalmanFitterInfo*>(fi);
    unsigned int nMeas = kfi->getNumMeasurements();
    std::vector<double> weights;
    weights.reserve(nMeas);
#ifdef DEBUG
      std::cout<<"(";
#endif
    for (unsigned int i=0; i<nMeas; ++i) {
      weights.push_back(kfi->getMeasurementOnPlane(i)->getWeight());
#ifdef DEBUG
      std::cout<<kfi->getMeasurementOnPlane(i)->getWeight();
      if (i<nMeas-1)
        std::cout<<", ";
#endif
    }
#ifdef DEBUG
      std::cout<<") ";
#endif
    weights_[rep].push_back(weights);
  }
#ifdef DEBUG
      std::cout<<"\n";
#endif
}


std::vector<std::vector<double> > DAF::calcWeights(Track* tr, const AbsTrackRep* rep, double beta) {

#ifdef DEBUG
      std::cout<<"DAF::calcWeights \n";
#endif

  std::vector<std::vector<double> > ret_val;

  std::vector< TrackPoint* > trackPoints = tr->getPointsWithMeasurement();
  for (std::vector< TrackPoint* >::iterator tp = trackPoints.begin(); tp != trackPoints.end(); ++tp) {

    AbsFitterInfo* fi = (*tp)->getFitterInfo(rep);
    if (dynamic_cast<KalmanFitterInfo*>(fi) == NULL){
      Exception exc("DAF::getWeights ==> can only use KalmanFitterInfos",__LINE__,__FILE__);
      throw exc;
    }
    KalmanFitterInfo* kfi = static_cast<KalmanFitterInfo*>(fi);
    unsigned int nMeas = kfi->getNumMeasurements();
    std::vector<double> weights;
    weights.reserve(nMeas);

    /*if(kfi->getStatusFlag() != 0) { // failed hit
      weights.assign(nMeas, 0.5);
      //std::cout<<"Assumed weight 0.5!!"<<std::endl;
      ret_val.push_back(weights);

#ifdef DEBUG
      std::cout<<"(";
#endif
      for (unsigned int j=0; j<nMeas; ++j){
        kfi->getMeasurementOnPlane(j)->setWeight(0.5);
#ifdef DEBUG
        std::cout<<"0.5";
        if (j<nMeas-1)
          std::cout<<", ";
#endif
      }
#ifdef DEBUG
      std::cout<<") ";
#endif
      continue;
    }*/

    std::vector<double> phi;
    double phi_sum = 0;
    double phi_cut = 0;
    const MeasuredStateOnPlane& smoothedState = kfi->getFittedState(true);

    TVectorD x_smoo(kfi->getMeasurementOnPlane()->getHMatrix() * smoothedState.getState());

    for(unsigned int j=0; j<nMeas; j++) {
      double* detV = new double(0);

      try{
        const MeasurementOnPlane* mop = kfi->getMeasurementOnPlane(j);
        int hitDim = mop->getState().GetNoElements();
        TMatrixDSym V( beta * mop->getCov());
        TVectorD resid(mop->getState() - x_smoo);
        TMatrixDSym Vinv;
        tools::invertMatrix(V, Vinv, detV); // can throw an Exception

        phi.push_back((1./(std::pow(2.*TMath::Pi(),hitDim/2)*sqrt(*detV)))*exp(-0.5*Vinv.Similarity(resid))); // std::pow(double, int) from <cmath> is faster than pow(double, double) from <math.h> when the exponent actually _is_ an integer.
        phi_sum += phi[j];
        //std::cerr << "hitDim " << hitDim << " fchi2Cuts[hitDim] " << fchi2Cuts[hitDim] << std::endl;
        double cutVal = chi2Cuts_[hitDim];
        assert(cutVal>1.E-6);
        //the follwing assumes that in the compeating hits (real hits in one DAF hit) could have different V otherwise calculation could be simplified
        phi_cut += (1./(std::pow(2.*TMath::Pi(),hitDim/2)*sqrt(*detV)))*exp(-0.5*cutVal/beta);
      }
      catch(Exception& e) {
        delete detV;
        std::cerr<<e.what();
        e.info();
        phi.push_back(0); //m and Vorig do not contain sensible values, assign weight 0
        continue;
      }

      delete detV;

    }

#ifdef DEBUG
      std::cout<<"(";
#endif
    for(unsigned int j=0; j<nMeas; j++) {
      double weight = phi[j]/(phi_sum+phi_cut);
      weights.push_back(weight);
      kfi->getMeasurementOnPlane(j)->setWeight(weight);
#ifdef DEBUG
      std::cout<<weight;
      if (j<nMeas-1)
        std::cout<<", ";
#endif
    }
#ifdef DEBUG
      std::cout<<") ";
#endif

    ret_val.push_back(weights);


  }
#ifdef DEBUG
      std::cout<<"\n";
#endif


  return ret_val;

}


} /* End of namespace genfit */
