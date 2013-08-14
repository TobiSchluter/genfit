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
  : AbsKalmanFitter(10), deltaWeight_(0.001)
{
  kalman_.reset(new KalmanFitterRefTrack());
  kalman_->setMultipleMeasurementHandling(weightedAverage);
  kalman_->setMaxIterations(1);
  static_cast<KalmanFitterRefTrack*>(kalman_.get())->setRefitAll();

  setAnnealingScheme(100, 0.1, 5); // also sets maxIterations_
  setProbCut(0.01);
}

DAF::DAF(AbsKalmanFitter* kalman)
  : AbsKalmanFitter(10), deltaWeight_(0.001)
{
  kalman_.reset(kalman);
  kalman_->setMaxIterations(1);

  if (dynamic_cast<KalmanFitterRefTrack*>(kalman_.get()) != NULL) {
    static_cast<KalmanFitterRefTrack*>(kalman_.get())->setRefitAll();
  }

  setAnnealingScheme(100, 0.1, 5); // also sets maxIterations_
  setProbCut(0.01);
}


void DAF::processTrack(Track* tr, const AbsTrackRep* rep, bool resortHits) {

#ifdef DEBUG
  std::cout<<"DAF::processTrack //////////////////////////////////////////////////////////////// \n";
#endif

  weights_.clear();
  std::vector<std::vector<double> > oldWeights;

  KalmanFitStatus* status;
  bool oneLastIter = false;

  unsigned int iBeta = 0;
  for(; iBeta != maxIterations_; ++iBeta) { // loop over. If no convergence is reached after 10 iterations just stop.

#ifdef DEBUG
      std::cout<<"DAF::processTrack, trackRep  " << rep << ", iteration " << iBeta << ", beta = " << betas_.at(iBeta) << "\n";
#endif

    kalman_->processTrack(tr, rep, resortHits);

    status = static_cast<KalmanFitStatus*>(tr->getFitStatus(rep));
    status->setIsFittedWithDaf();


    // check break conditions

    if (! status->isFitted()){
      #ifdef DEBUG
      std::cout << "DAF::Kalman could not fit!\n";
      #endif
      status->setIsFitted(false);
      break;
    }

    if( oneLastIter == true){
      #ifdef DEBUG
      std::cout << "DAF::break after one last iteration\n";
      #endif
      status->setIsFitConverged();
      break;
    }

    if(iBeta == maxIterations_-1 ){
      status->setIsFitConverged(false);
      #ifdef DEBUG
      std::cout << "DAF::number of max iterations reached!\n";
      #endif
      break;
    }


    // get and update weights
    getWeights(tr, rep);
    if (iBeta > 0)
      oldWeights = weights_;

    try{
      weights_ = calcWeights(tr, rep, betas_.at(iBeta));
    } catch(Exception& e) {
      std::cerr<<e.what();
      e.info();
      //std::cerr << "calc weights failed" << std::endl;
      //mini_trk->getTrackRep(0)->setStatusFlag(1);
      status->setIsFitted(false);
      status->setIsFitConverged(false);
      break;
    }

    // check if converged
    if (iBeta > 0) {
      if ( isConvergent(oldWeights, rep) ){
        #ifdef DEBUG
        std::cout << "DAF::convergence reached in iteration " << iBeta << " -> Do one last iteration with updated weights.\n";
        #endif
        oneLastIter = true;
        status->setIsFitConverged();
      }
    }

    // check if fit is failing utterly
    /*if (status->getForwardNdf() < 0. && betas_.at(iBeta) < 1.) {
      #ifdef DEBUG
      std::cout << "DAF:: NDF < 0; skip track! \n";
      #endif
      status->setIsFitConverged(false);
      break;
    }*/

  } // end loop over betas

  status->setNumIterations(iBeta+1);

  if (status->getForwardPVal() == 0. &&
      status->getBackwardPVal() == 0.) {
    status->setIsFitConverged(false);
  }

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
  maxIterations_ = betas_.size() + 4;
  betas_.resize(maxIterations_,betas_.back()); //make sure main loop has a maximum of 10 iterations and also make sure the last beta value is used for if more iterations are needed then the ones set by the user.
}


void DAF::setAnnealingScheme(double bStart, double bFinal, unsigned int nSteps) {
  assert(bStart > bFinal);
  assert(bFinal > 1.E-10);
  assert(nSteps > 1);

  maxIterations_ = nSteps + 4;

  betas_.clear();

  const double x = log(bStart/bFinal)/log(2)/(nSteps-1);

  for (unsigned int i=0; i<nSteps; ++i) {
    double exp = double(nSteps-i-1)*x;
    betas_.push_back(bFinal * pow(2., exp));
  }

  betas_.resize(maxIterations_,betas_.back()); //make sure main loop has a maximum of 10 iterations and also make sure the last beta value is used for if more iterations are needed then the ones set by the user.

  /*for (unsigned int i=0; i<betas_.size(); ++i) {
    std::cout<< betas_.at(i) << ", ";
  }*/
}


bool DAF::isConvergent(const std::vector<std::vector<double> >& oldWeights, const AbsTrackRep* rep) const {
  const int n = oldWeights.size();
  assert(n == int(weights_.size()));
  for( int i = 0; i != n; ++i){
    const int m = oldWeights[i].size();
    if (m != int(weights_[i].size())) {
      std::cout << "m = " << m << ", newWeights[i].size() = " << weights_[i].size() << std::endl;
      printWeights(oldWeights);
      printWeights(weights_);
      assert(m == int(weights_[i].size()));
    }
    for( int j = 0; j != m; ++j){
      if( fabs(oldWeights[i][j] - weights_[i][j]) > deltaWeight_ ){
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

  weights_.clear();

  const std::vector< TrackPoint* >& trackPoints = tr->getPointsWithMeasurement();
  for (std::vector< TrackPoint* >::const_iterator tp = trackPoints.begin(); tp != trackPoints.end(); ++tp) {
    if (! (*tp)->hasFitterInfo(rep)) {
      continue;
    }
    AbsFitterInfo* fi = (*tp)->getFitterInfo(rep);
    if (dynamic_cast<KalmanFitterInfo*>(fi) == NULL){
      Exception exc("DAF::getWeights ==> can only use KalmanFitterInfos",__LINE__,__FILE__);
      throw exc;
    }
    KalmanFitterInfo* kfi = static_cast<KalmanFitterInfo*>(fi);
    unsigned int nMeas = kfi->getNumMeasurements();
    std::vector<double> weights;
    weights.reserve(nMeas);
    for (unsigned int i=0; i<nMeas; ++i) {
      weights.push_back(kfi->getMeasurementOnPlane(i)->getWeight());
    }
    weights_.push_back(weights);
  }

#ifdef DEBUG
  printWeights(weights_);
#endif

}


std::vector<std::vector<double> > DAF::calcWeights(Track* tr, const AbsTrackRep* rep, double beta) {

#ifdef DEBUG
      std::cout<<"DAF::calcWeights \n";
#endif

  std::vector<std::vector<double> > ret_val;

  const std::vector< TrackPoint* >& trackPoints = tr->getPointsWithMeasurement();
  for (std::vector< TrackPoint* >::const_iterator tp = trackPoints.begin(); tp != trackPoints.end(); ++tp) {
    if (! (*tp)->hasFitterInfo(rep)) {
      continue;
    }
    AbsFitterInfo* fi = (*tp)->getFitterInfo(rep);
    if (dynamic_cast<KalmanFitterInfo*>(fi) == NULL){
      Exception exc("DAF::getWeights ==> can only use KalmanFitterInfos",__LINE__,__FILE__);
      throw exc;
    }
    KalmanFitterInfo* kfi = static_cast<KalmanFitterInfo*>(fi);
    unsigned int nMeas = kfi->getNumMeasurements();

    std::vector<double> phi(nMeas, 0.);
    double phi_sum = 0;
    double phi_cut = 0;
    const MeasuredStateOnPlane& smoothedState = kfi->getFittedState();

    // This assumes that all measurements on the plane have the same dimensionality.
    TVectorD x_smoo(kfi->getMeasurementOnPlane()->getHMatrix() * smoothedState.getState());
    int hitDim = x_smoo.GetNrows();
    double twoPiN = std::pow(2.*M_PI, hitDim);

    for(unsigned int j=0; j<nMeas; j++) {

      try{
        const MeasurementOnPlane* mop = kfi->getMeasurementOnPlane(j);
	assert(mop->getState().GetNrows() == hitDim);
        TVectorD resid(mop->getState() - x_smoo);
        TMatrixDSym Vinv(mop->getCov());
	double detV;
        tools::invertMatrix(Vinv, &detV); // can throw an Exception

        double chi2 = Vinv.Similarity(resid);
#ifdef DEBUG
        std::cout<<"chi2 = " << chi2 << "\n";
#endif

	double norm = 1./sqrt(twoPiN * detV);

        phi[j] = norm*exp(-0.5*chi2/beta);
        phi_sum += phi[j];
        //std::cerr << "hitDim " << hitDim << " fchi2Cuts[hitDim] " << fchi2Cuts[hitDim] << std::endl;
        double cutVal = chi2Cuts_[hitDim];
        assert(cutVal>1.E-6);
        //the following assumes that in the competing hits (real hits in one DAF hit) could have different V otherwise calculation could be simplified
        phi_cut += norm*exp(-0.5*cutVal/beta);
      }
      catch(Exception& e) {
        std::cerr << e.what();
        e.info();
      }
    }

    std::vector<double> weights(nMeas, 0.);
    for(unsigned int j=0; j<nMeas; j++) {
      double weight = phi[j]/(phi_sum+phi_cut);
      weights[j] = weight;
      kfi->getMeasurementOnPlane(j)->setWeight(weight);
    }

    ret_val.push_back(weights);
  }

#ifdef DEBUG
  printWeights(ret_val);
#endif

  return ret_val;
}


void DAF::printWeights(const std::vector<std::vector<double> >& weights) const {
  for (unsigned int i=0; i< weights.size(); ++i){
    std::cout << "(";
      for (unsigned int j=0; j<weights[i].size(); ++j) {
        if (j > 0)
          std::cout << ", ";
        std::cout << weights[i][j];
      }
    std::cout << ") ";
  }
  std::cout << "\n";
}


// Customized from generated Streamer.
void DAF::Streamer(TBuffer &R__b)
{
   // Stream an object of class genfit::DAF.

   //This works around a msvc bug and should be harmless on other platforms
   typedef ::genfit::DAF thisClass;
   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      //This works around a msvc bug and should be harmless on other platforms
      typedef genfit::AbsKalmanFitter baseClass0;
      baseClass0::Streamer(R__b);
      R__b >> deltaWeight_;
      // weights_ are only of intermediate use -> not saved
      {
         std::vector<double> &R__stl =  betas_;
         R__stl.clear();
         int R__i, R__n;
         R__b >> R__n;
         R__stl.reserve(R__n);
         for (R__i = 0; R__i < R__n; R__i++) {
            double R__t;
            R__b >> R__t;
            R__stl.push_back(R__t);
         }
      }
      {
         std::map<int,double> &R__stl =  chi2Cuts_;
         R__stl.clear();
         int R__i, R__n;
         R__b >> R__n;
         for (R__i = 0; R__i < R__n; R__i++) {
            int R__t;
            R__b >> R__t;
            double R__t2;
            R__b >> R__t2;
            typedef int Value_t;
            std::pair<Value_t const, double > R__t3(R__t,R__t2);
            R__stl.insert(R__t3);
         }
      }
      AbsKalmanFitter *p;
      R__b >> p;
      kalman_.reset(p);
      R__b.CheckByteCount(R__s, R__c, thisClass::IsA());
   } else {
      R__c = R__b.WriteVersion(thisClass::IsA(), kTRUE);
      //This works around a msvc bug and should be harmless on other platforms
      typedef genfit::AbsKalmanFitter baseClass0;
      baseClass0::Streamer(R__b);
      R__b << deltaWeight_;
      // weights_ are only of intermediate use -> not saved
      weights_.clear();
      {
         std::vector<double> &R__stl =  betas_;
         int R__n=(&R__stl) ? int(R__stl.size()) : 0;
         R__b << R__n;
         if(R__n) {
            std::vector<double>::iterator R__k;
            for (R__k = R__stl.begin(); R__k != R__stl.end(); ++R__k) {
            R__b << (*R__k);
            }
         }
      }
      {
         std::map<int,double> &R__stl =  chi2Cuts_;
         int R__n=(&R__stl) ? int(R__stl.size()) : 0;
         R__b << R__n;
         if(R__n) {
            std::map<int,double>::iterator R__k;
            for (R__k = R__stl.begin(); R__k != R__stl.end(); ++R__k) {
            R__b << ((*R__k).first );
            R__b << ((*R__k).second);
            }
         }
      }
      R__b << kalman_.get();
      R__b.SetByteCount(R__c, kTRUE);
   }
}

} /* End of namespace genfit */
