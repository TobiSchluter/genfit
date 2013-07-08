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
#include "KalmanFitterRefTrack.h"
#include "Tools.h"
#include "Exception.h"
#include <assert.h>
#include <cmath>

//root stuff
#include <TMath.h>
#include <Math/QuantFuncMathCore.h>


namespace genfit {

DAF::DAF() {
  kalman_.reset(new KalmanFitterRefTrack());
  kalman_->setMultipleMeasurementHandling(weightedAverage);
  kalman_->setNumIterations(1);

  setBetas(81.,8.,4.,1.,1.,1.);
  setProbCut(0.01);
}

DAF::DAF(AbsKalmanFitter* kalman) {
  kalman_.reset(kalman);
  kalman_->setNumIterations(1);

  setBetas(81.,8.,4.,1.,1.,1.);
  setProbCut(0.01);
}


void DAF::processTrack(Track* tr, const AbsTrackRep* rep) {

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


bool DAF::isConvergent(const std::vector<std::vector<double> >& oldWeights, AbsTrackRep* rep) const {
  const int n = oldWeights.size();
  const std::vector<std::vector<double> >& newWeights = weights_.at(rep);
  assert(n == int(newWeights.size()));
  for( int i = 0; i != n; ++i){
    const int m = oldWeights[i].size();
    assert(m == int(newWeights[i].size()));
    for( int j = 0; j != m; ++j){
      if( fdim( oldWeights[i][j] , newWeights[i][j]) > 0.001 ){ //Moritz just made the value up. has to be tested if good
        return false;
      }
    }
  }
  return true;
}


} /* End of namespace genfit */
