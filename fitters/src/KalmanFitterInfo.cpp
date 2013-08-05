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

#include "KalmanFitterInfo.h"

#include <cassert>
#include <iostream>

#include "Exception.h"
#include "Tools.h"
#include "Track.h"
#include "TrackPoint.h"

//#define DEBUG


namespace genfit {

KalmanFitterInfo::KalmanFitterInfo() :
  AbsFitterInfo()
{
  ;
}

KalmanFitterInfo::KalmanFitterInfo(const TrackPoint* trackPoint, const AbsTrackRep* rep)  :
  AbsFitterInfo(trackPoint, rep)
{
  ;
}

KalmanFitterInfo::~KalmanFitterInfo() {
  // FIXME: need smart pointers / smart containers here
  for (size_t i = 0; i < measurementsOnPlane_.size(); ++i)
    delete measurementsOnPlane_[i];
}


KalmanFitterInfo* KalmanFitterInfo::clone() const {
  KalmanFitterInfo* retVal = new KalmanFitterInfo(this->getTrackPoint(), this->getRep());
  if (hasReferenceState())
    retVal->setReferenceState(new ReferenceStateOnPlane(*getReferenceState()));
  if (hasForwardPrediction())
    retVal->setForwardPrediction(new MeasuredStateOnPlane(*getForwardPrediction()));
  if (hasForwardUpdate())
    retVal->setForwardUpdate(new KalmanFittedStateOnPlane(*getForwardUpdate()));
  if (hasBackwardPrediction())
    retVal->setBackwardPrediction(new MeasuredStateOnPlane(*getBackwardPrediction()));
  if (hasBackwardUpdate())
    retVal->setBackwardUpdate(new KalmanFittedStateOnPlane(*getBackwardUpdate()));

  retVal->measurementsOnPlane_.reserve(getNumMeasurements());
  for (std::vector<MeasurementOnPlane*>::const_iterator it = this->measurementsOnPlane_.begin(); it != this->measurementsOnPlane_.end(); ++it) {
    retVal->addMeasurementOnPlane(new MeasurementOnPlane(**it));
  }

  return retVal;
}


std::vector< genfit::MeasurementOnPlane* > KalmanFitterInfo::getMeasurementsOnPlane() const {
  std::vector< genfit::MeasurementOnPlane* > retVal;
  retVal.reserve(measurementsOnPlane_.size());

  for (std::vector<MeasurementOnPlane*>::const_iterator it = measurementsOnPlane_.begin(); it != measurementsOnPlane_.end(); ++it) {
    retVal.push_back(*it);
  }

  return retVal;
}


MeasurementOnPlane KalmanFitterInfo::getAvgWeightedMeasurementOnPlane() const {

  MeasurementOnPlane retVal(*(measurementsOnPlane_[0]));
  retVal.setWeight(1);

  if(measurementsOnPlane_.size() == 1) {
    double weight = (measurementsOnPlane_[0])->getWeight();
    if (weight != 1.) {
      retVal.getCov() *= 1. / weight;
    }
  }
  else { // more than one hit
    retVal.getState().Zero();
    retVal.getCov().Zero();

    TMatrixDSym covInv;
    std::vector<TMatrixDSym> weightedCovInvs;

    for(unsigned int i=0; i<measurementsOnPlane_.size(); ++i) {

      if (i>0) {
        // make sure we have compatible measurement types
        // TODO: replace with Exceptions!
        assert(measurementsOnPlane_[i]->getPlane() == measurementsOnPlane_[0]->getPlane());
        assert(measurementsOnPlane_[i]->getHMatrix() == measurementsOnPlane_[0]->getHMatrix());
      }

      tools::invertMatrix(measurementsOnPlane_[i]->getCov(), covInv); // invert cov
      covInv *= measurementsOnPlane_[i]->getWeight(); // weigh cov
      weightedCovInvs.push_back(covInv); // cov is already inverted and weighted
      retVal.getCov() += covInv; // cov is already inverted and weighted
    }

    // invert fHitCov
    tools::invertMatrix(retVal.getCov());

    //set the weighted-mean coord
    for(unsigned int i=0; i<measurementsOnPlane_.size(); ++i) {
      retVal.getState() += weightedCovInvs[i] * measurementsOnPlane_[i]->getState();
    }
    retVal.getState() *= retVal.getCov();
  }

  return retVal;
}


const MeasurementOnPlane* KalmanFitterInfo::getClosestMeasurementOnPlane(const StateOnPlane* sop) const {
  if(measurementsOnPlane_.size() == 0)
    return NULL;

  if(measurementsOnPlane_.size() == 1)
    return getMeasurementOnPlane(0);

  double normMin(9.99E99);
  unsigned int iMin(0);
  for (unsigned int i=0; i<getNumMeasurements(); ++i) {
    const TMatrixD& H = measurementsOnPlane_[i]->getHMatrix();
    TVectorD res = measurementsOnPlane_[i]->getState() - (H * sop->getState());
    double norm = sqrt(res.Norm2Sqr());
    if (norm < normMin) {
      normMin = norm;
      iMin = i;
    }
  }

  return getMeasurementOnPlane(iMin);
}


std::vector<double> KalmanFitterInfo::getWeights() const {
  std::vector<double> retVal;

  for (unsigned int i=0; i<getNumMeasurements(); ++i) {
    retVal.push_back(getMeasurementOnPlane(i)->getWeight());
  }

  return retVal;
}


MeasuredStateOnPlane* KalmanFitterInfo::getFittedState(bool biased) const {
  if (biased) {
    if (fittedStateBiased_.get() != NULL)
      return fittedStateBiased_.get();

    if (this->getTrackPoint()->getTrack()->getPointWithMeasurement(-1) == this->getTrackPoint()) {// last measurement
      assert(forwardUpdate_.get() != NULL);
      #ifdef DEBUG
      std::cout << "KalmanFitterInfo::getFittedState - biased at last measurement = forwardUpdate_ \n";
      #endif
      return forwardUpdate_.get();
    }
    else if (this->getTrackPoint()->getTrack()->getPointWithMeasurement(0) == this->getTrackPoint()) { // first measurement
      assert(backwardUpdate_.get() != NULL);
      #ifdef DEBUG
      std::cout << "KalmanFitterInfo::getFittedState - biased at first measurement = backwardUpdate_ \n";
      backwardUpdate_->Print();
      #endif
      return backwardUpdate_.get();
    }

    assert(forwardUpdate_.get() != NULL);
    assert(backwardPrediction_.get() != NULL);
    #ifdef DEBUG
    std::cout << "KalmanFitterInfo::getFittedState - biased = mean(forwardUpdate_, backwardPrediction_) \n";
    #endif
    fittedStateBiased_.reset(new MeasuredStateOnPlane(calcAverageState(forwardUpdate_.get(), backwardPrediction_.get())));
    return fittedStateBiased_.get();
  }
  else { // unbiased
    if (fittedStateUnbiased_.get() != NULL)
      return fittedStateUnbiased_.get();

    if (this->getTrackPoint()->getTrack()->getPointWithMeasurement(-1) == this->getTrackPoint()) { // last measurement
      assert(forwardPrediction_.get() != NULL);
      #ifdef DEBUG
      std::cout << "KalmanFitterInfo::getFittedState - unbiased at last measurement = forwardPrediction_ \n";
      #endif
      return forwardPrediction_.get();
    }
    else if (this->getTrackPoint()->getTrack()->getPointWithMeasurement(0) == this->getTrackPoint()) { // first measurement
      assert(backwardPrediction_.get() != NULL);
      #ifdef DEBUG
      std::cout << "KalmanFitterInfo::getFittedState - unbiased at first measurement = backwardPrediction_ \n";
      #endif
      return backwardPrediction_.get();
    }

    assert(forwardPrediction_.get() != NULL);
    assert(backwardPrediction_.get() != NULL);
    #ifdef DEBUG
    std::cout << "KalmanFitterInfo::getFittedState - unbiased = mean(forwardPrediction_, backwardPrediction_) \n";
    #endif
    fittedStateUnbiased_.reset(new MeasuredStateOnPlane(calcAverageState(forwardPrediction_.get(), backwardPrediction_.get())));
    return fittedStateUnbiased_.get();
  }

}


MeasurementOnPlane KalmanFitterInfo::getResidual(bool biased, unsigned int iMeasurement) const {
  // TODO: Test

  const MeasuredStateOnPlane* smoothedState = getFittedState(biased);
  const MeasurementOnPlane* measurement = measurementsOnPlane_.at(iMeasurement);
  SharedPlanePtr plane = measurement->getPlane();

  // check equality of planes and reps
  if(*(smoothedState->getPlane()) != *plane) {
    Exception e("KalmanFitterInfo::getResidual: smoothedState and measurement are not defined in the same plane.", __LINE__,__FILE__);
    throw e;
  }
  if(smoothedState->getRep() != measurement->getRep()) {
    Exception e("KalmanFitterInfo::getResidual: smoothedState and measurement are not defined wrt the same TrackRep.", __LINE__,__FILE__);
    throw e;
  }

  const TMatrixD& H = measurement->getHMatrix();

  TVectorD res = measurement->getState() - (H * smoothedState->getState());

  TMatrixDSym cov(smoothedState->getCov());
  cov.Similarity(H);
  cov += measurement->getCov();

  return MeasurementOnPlane(res, cov, plane, smoothedState->getRep(), H, measurement->getWeight());
}


void KalmanFitterInfo::setForwardPrediction(MeasuredStateOnPlane* forwardPrediction) {
  forwardPrediction_.reset(forwardPrediction);
  fittedStateUnbiased_.reset();
  fittedStateBiased_.reset();
  if (forwardPrediction_)
    setPlane(forwardPrediction_->getPlane());
}

void KalmanFitterInfo::setBackwardPrediction(MeasuredStateOnPlane* backwardPrediction) {
  backwardPrediction_.reset(backwardPrediction);
  fittedStateUnbiased_.reset();
  fittedStateBiased_.reset();
  if (backwardPrediction_)
    setPlane(backwardPrediction_->getPlane());
}

void KalmanFitterInfo::setForwardUpdate(KalmanFittedStateOnPlane* forwardUpdate) {
  forwardUpdate_.reset(forwardUpdate);
  fittedStateUnbiased_.reset();
  fittedStateBiased_.reset();
  if (forwardUpdate_)
    setPlane(forwardUpdate_->getPlane());
}

void KalmanFitterInfo::setBackwardUpdate(KalmanFittedStateOnPlane* backwardUpdate) {
  backwardUpdate_.reset(backwardUpdate);
  fittedStateUnbiased_.reset();
  fittedStateBiased_.reset();
  if (backwardUpdate_)
    setPlane(backwardUpdate_->getPlane());
}


void KalmanFitterInfo::setMeasurementsOnPlane(const std::vector< genfit::MeasurementOnPlane* >& measurementsOnPlane) {
  measurementsOnPlane_.clear();

  for (std::vector<MeasurementOnPlane*>::const_iterator m = measurementsOnPlane.begin(), mend = measurementsOnPlane.end(); m < mend; ++m) {
    addMeasurementOnPlane(*m);
  }

  if (measurementsOnPlane.size() != 0)
    setPlane(measurementsOnPlane[0]->getPlane());
}


void KalmanFitterInfo::addMeasurementsOnPlane(const std::vector< genfit::MeasurementOnPlane* >& measurementsOnPlane) {
  for (std::vector<MeasurementOnPlane*>::const_iterator m = measurementsOnPlane.begin(), mend = measurementsOnPlane.end(); m < mend; ++m) {
    addMeasurementOnPlane(*m);
  }
}


void KalmanFitterInfo::setRep(const AbsTrackRep* rep) {
  rep_ = rep;

  if (referenceState_)
    referenceState_->setRep(rep);

  if (forwardPrediction_)
    forwardPrediction_->setRep(rep);

  if (forwardUpdate_)
    forwardUpdate_->setRep(rep);

  if (backwardPrediction_)
    backwardPrediction_->setRep(rep);

  if (backwardUpdate_)
    backwardUpdate_->setRep(rep);

  for (std::vector<MeasurementOnPlane*>::iterator it = measurementsOnPlane_.begin(); it != measurementsOnPlane_.end(); ++it) {
    (*it)->setRep(rep);
  }
}


void KalmanFitterInfo::setWeights(const std::vector<double>& weights) {

  if (weights.size() != getNumMeasurements()) {
    Exception e("KalmanFitterInfo::setWeights: weights do not have the same size as measurementsOnPlane", __LINE__,__FILE__);
    throw e;
  }

  for (unsigned int i=0; i<getNumMeasurements(); ++i) {
    getMeasurementOnPlane(i)->setWeight(weights[i]);
  }
}


MeasuredStateOnPlane KalmanFitterInfo::calcAverageState(const MeasuredStateOnPlane* forwardState, const MeasuredStateOnPlane* backwardState) const {
  if (forwardState == NULL || backwardState == NULL) {
    Exception e("KalmanFitterInfo::calcAverageState: forwardState or backwardState is NULL.", __LINE__,__FILE__);
    throw e;
  }
  // check if both states are defined in the same plane
  if (forwardState->getPlane() != backwardState->getPlane()) {
    Exception e("KalmanFitterInfo::calcAverageState: forwardState and backwardState are not defined in the same plane.", __LINE__,__FILE__);
    throw e;
  }

  TMatrixDSym fCovInv, bCovInv, smoothed_cov;

  tools::invertMatrix(forwardState->getCov(), fCovInv);
  tools::invertMatrix(backwardState->getCov(), bCovInv);

  tools::invertMatrix(fCovInv + bCovInv, smoothed_cov);

  MeasuredStateOnPlane retVal(*forwardState); // copies auxInfo
  retVal.setState(smoothed_cov * (fCovInv*forwardState->getState() + bCovInv*backwardState->getState()));
  retVal.setCov(smoothed_cov);

  return retVal;
}


void KalmanFitterInfo::Print(const Option_t*) const {
  std::cout << "genfit::KalmanFitterInfo. Belongs to TrackPoint " << trackPoint_ << "; TrackRep " <<  rep_  << "\n";

  for (unsigned int i=0; i<measurementsOnPlane_.size(); ++i) {
    std::cout << "MeasurementOnPlane Nr " << i <<": "; measurementsOnPlane_[i]->Print();
  }

  if (referenceState_) {
    std::cout << "Reference state: "; referenceState_->Print();
  }
  if (forwardPrediction_) {
    std::cout << "Forward prediction_: "; forwardPrediction_->Print();
  }
  if (forwardUpdate_) {
    std::cout << "Forward update: "; forwardUpdate_->Print();
  }
  if (backwardPrediction_) {
    std::cout << "Backward prediction_: "; backwardPrediction_->Print();
  }
  if (backwardUpdate_) {
    std::cout << "Backward update: "; backwardUpdate_->Print();
  }

}


bool KalmanFitterInfo::checkConsistency() const {

  // check if in a TrackPoint
  if (!trackPoint_) {
    std::cerr << "KalmanFitterInfo::checkConsistency(): trackPoint_ is NULL" << std::endl;
    return false;
  }

  // check if there is a reference state
  /*if (!referenceState_) {
    std::cerr << "KalmanFitterInfo::checkConsistency(): referenceState_ is NULL" << std::endl;
    return false;
  }*/

  SharedPlanePtr plane = getPlane();

  if (plane.get() == NULL) {
    if (!(referenceState_ || forwardPrediction_ || forwardUpdate_ || backwardPrediction_ || backwardUpdate_ || measurementsOnPlane_.size() > 0))
      return true;
    std::cerr << "KalmanFitterInfo::checkConsistency(): plane is NULL" << std::endl;
    return false;
  }

  TVector3 oTest = plane->getO(); // see if the plane object is there
  oTest *= 47.11;

  // if more than one measurement, check if they have the same dimensionality
  if (measurementsOnPlane_.size() > 1) {
    int dim = measurementsOnPlane_[0]->getState().GetNrows();
    for (unsigned int i = 1; i<measurementsOnPlane_.size(); ++i) {
      if(measurementsOnPlane_[i]->getPlane() != plane) {
        std::cerr << "KalmanFitterInfo::checkConsistency(): measurementsOnPlane_ are not all defined with the correct plane" << std::endl;
        return false;
      }
      if(measurementsOnPlane_[i]->getState().GetNrows() != dim) {
        std::cerr << "KalmanFitterInfo::checkConsistency(): measurementsOnPlane_ do not all have the same dimensionality" << std::endl;
        return false;
      }
    }
    if (dim == 0) {
      std::cerr << "KalmanFitterInfo::checkConsistency(): measurementsOnPlane_ have dimension 0" << std::endl;
      return false;
    }
  }

  // see if everything else is defined wrt this plane and rep_
  int dim = rep_->getDim(); // check dim
  if (referenceState_) {
    if(referenceState_->getPlane() != plane) {
      std::cerr << "KalmanFitterInfo::checkConsistency(): referenceState_ is not defined with the correct plane" << std::endl;
      return false;
    }
    if (referenceState_->getRep() != rep_) {
      std::cerr << "KalmanFitterInfo::checkConsistency(): referenceState_ is not defined with the correct TrackRep" << std::endl;
      return false;
    }
    if (referenceState_->getState().GetNrows() != dim) {
      std::cerr << "KalmanFitterInfo::checkConsistency(): referenceState_ does not have the right dimension!" << std::endl;
      return false;
    }
  }

  if (forwardPrediction_) {
    if(forwardPrediction_->getPlane() != plane) {
      std::cerr << "KalmanFitterInfo::checkConsistency(): forwardPrediction_ is not defined with the correct plane" << std::endl;
      return false;
    }
    if(forwardPrediction_->getRep() != rep_) {
      std::cerr << "KalmanFitterInfo::checkConsistency(): forwardPrediction_ is not defined with the correct TrackRep" << std::endl;
      return false;
    }
    if (forwardPrediction_->getState().GetNrows() != dim || forwardPrediction_->getCov().GetNrows() != dim) {
      std::cerr << "KalmanFitterInfo::checkConsistency(): forwardPrediction_ does not have the right dimension!" << std::endl;
      return false;
    }
  }
  if (forwardUpdate_) {
    if(forwardUpdate_->getPlane() != plane) {
      std::cerr << "KalmanFitterInfo::checkConsistency(): forwardUpdate_ is not defined with the correct plane" << std::endl;
      return false;
    }
    if(forwardUpdate_->getRep() != rep_) {
      std::cerr << "KalmanFitterInfo::checkConsistency(): forwardUpdate_ is not defined with the correct TrackRep" << std::endl;
      return false;
    }
    if (forwardUpdate_->getState().GetNrows() != dim || forwardUpdate_->getCov().GetNrows() != dim) {
      std::cerr << "KalmanFitterInfo::checkConsistency(): forwardUpdate_ does not have the right dimension!" << std::endl;
      return false;
    }
  }

  if (backwardPrediction_) {
    if(backwardPrediction_->getPlane() != plane) {
      std::cerr << "KalmanFitterInfo::checkConsistency(): backwardPrediction_ is not defined with the correct plane" << std::endl;
      return false;
    }
    if(backwardPrediction_->getRep() != rep_) {
      std::cerr << "KalmanFitterInfo::checkConsistency(): backwardPrediction_ is not defined with the correct TrackRep" << std::endl;
      return false;
    }
    if (backwardPrediction_->getState().GetNrows() != dim || backwardPrediction_->getCov().GetNrows() != dim) {
      std::cerr << "KalmanFitterInfo::checkConsistency(): backwardPrediction_ does not have the right dimension!" << std::endl;
      return false;
    }
  }
  if (backwardUpdate_) {
    if(backwardUpdate_->getPlane() != plane) {
      std::cerr << "KalmanFitterInfo::checkConsistency(): backwardUpdate_ is not defined with the correct plane" << std::endl;
      return false;
    }
    if(backwardUpdate_->getRep() != rep_) {
      std::cerr << "KalmanFitterInfo::checkConsistency(): backwardUpdate_ is not defined with the correct TrackRep" << std::endl;
      return false;
    }
    if (backwardUpdate_->getState().GetNrows() != dim || backwardUpdate_->getCov().GetNrows() != dim) {
      std::cerr << "KalmanFitterInfo::checkConsistency(): backwardUpdate_ does not have the right dimension!" << std::endl;
      return false;
    }
  }

  for (std::vector<MeasurementOnPlane*>::const_iterator it = measurementsOnPlane_.begin(); it != measurementsOnPlane_.end(); ++it) {
    if((*it)->getPlane() != plane) {
      std::cerr << "KalmanFitterInfo::checkConsistency(): measurement is not defined with the correct plane" << std::endl;
      return false;
    }
    if((*it)->getRep() != rep_) {
      std::cerr << "KalmanFitterInfo::checkConsistency(): measurement is not defined with the correct TrackRep" << std::endl;
      return false;
    }
    if ((*it)->getState().GetNrows() == 0) {
      std::cerr << "KalmanFitterInfo::checkConsistency(): measurement has dimension 0!" << std::endl;
      return false;
    }
  }

  // see if there is an update w/o prediction or measurement
  if (forwardUpdate_ && !forwardPrediction_) {
    std::cerr << "KalmanFitterInfo::checkConsistency(): forwardUpdate_ w/o forwardPrediction_" << std::endl;
    return false;
  }

  if (forwardUpdate_ && measurementsOnPlane_.size() == 0) {
    std::cerr << "KalmanFitterInfo::checkConsistency(): forwardUpdate_ w/o measurement" << std::endl;
    return false;
  }


  if (backwardUpdate_ && !backwardPrediction_) {
    std::cerr << "KalmanFitterInfo::checkConsistency(): backwardUpdate_ w/o backwardPrediction_" << std::endl;
    return false;
  }

  if (backwardUpdate_ && measurementsOnPlane_.size() == 0) {
    std::cerr << "KalmanFitterInfo::checkConsistency(): backwardUpdate_ w/o measurement" << std::endl;
    return false;
  }


  return true;
}


// Modified from auto-generated Streamer to correctly deal with smart pointers.
void KalmanFitterInfo::Streamer(TBuffer &R__b)
{
   // Stream an object of class genfit::KalmanFitterInfo.

   //This works around a msvc bug and should be harmless on other platforms
   typedef ::genfit::KalmanFitterInfo thisClass;
   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      //This works around a msvc bug and should be harmless on other platforms
      typedef genfit::AbsFitterInfo baseClass0;
      baseClass0::Streamer(R__b);
      {
	ReferenceStateOnPlane *p = 0;
	R__b >> p;
	p->setRep(this->getRep());
	referenceState_.reset(p);
      }
      {
	MeasuredStateOnPlane *p = 0;
	R__b >> p;
	p->setRep(this->getRep());
	forwardPrediction_.reset(p);
      }
      {
	KalmanFittedStateOnPlane *p = 0;
	R__b >> p;
	p->setRep(this->getRep());
	forwardUpdate_.reset(p);
      }
      {
	MeasuredStateOnPlane *p = 0;
	R__b >> p;
	p->setRep(this->getRep());
	backwardPrediction_.reset(p);
      }
      {
	KalmanFittedStateOnPlane *p = 0;
	R__b >> p;
	p->setRep(this->getRep());
	backwardUpdate_.reset(p);
      }
      {
	std::vector<genfit::MeasurementOnPlane*,std::allocator<genfit::MeasurementOnPlane*> > &R__stl =  measurementsOnPlane_;
         R__stl.clear();
         TClass *R__tcl1 = TBuffer::GetClass(typeid(genfit::MeasurementOnPlane));
         if (R__tcl1==0) {
            Error("measurementsOnPlane_ streamer","Missing the TClass object for genfit::MeasurementOnPlane!");
            return;
         }
         int R__i, R__n;
         R__b >> R__n;
         R__stl.reserve(R__n);
         for (R__i = 0; R__i < R__n; R__i++) {
            genfit::MeasurementOnPlane* R__t;
            R__t = (genfit::MeasurementOnPlane*)R__b.ReadObjectAny(R__tcl1);
            R__stl.push_back(R__t);
         }
      }
      R__b.CheckByteCount(R__s, R__c, thisClass::IsA());
   } else {
      R__c = R__b.WriteVersion(thisClass::IsA(), kTRUE);
      //This works around a msvc bug and should be harmless on other platforms
      typedef genfit::AbsFitterInfo baseClass0;
      baseClass0::Streamer(R__b);
      R__b << referenceState_.get();
      R__b << forwardPrediction_.get();
      R__b << forwardUpdate_.get();
      R__b << backwardPrediction_.get();
      R__b << backwardUpdate_.get();
      {
	std::vector<genfit::MeasurementOnPlane*,std::allocator<genfit::MeasurementOnPlane*> > &R__stl =  measurementsOnPlane_;
         int R__n=(&R__stl) ? int(R__stl.size()) : 0;
         R__b << R__n;
         if(R__n) {
	   std::vector<genfit::MeasurementOnPlane*,std::allocator<genfit::MeasurementOnPlane*> >::iterator R__k;
            for (R__k = R__stl.begin(); R__k != R__stl.end(); ++R__k) {
            R__b << (*R__k);
            }
         }
      }
      R__b.SetByteCount(R__c, kTRUE);
   }
}


} /* End of namespace genfit */
