/* Copyright 2013, Ludwig-Maximilians Universität München,
   Authors: Tobias Schlüter & Johannes Rauch

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

/* This implements the simple Kalman fitter with no reference track
   that uses the stateSeed only until it forgets about it after the
   first few hits.  */

#include "KalmanFitter.h"

#include "Exception.h"
#include "KalmanFitterInfo.h"
#include "KalmanFitStatus.h"
#include "Track.h"
#include "TrackPoint.h"
#include "Tools.h"

#include <Math/ProbFunc.h>
#include <algorithm>


//#define DEBUG

using namespace genfit;


void KalmanFitter::fitTrack(Track* tr, const AbsTrackRep* rep,
    double& chi2, double& ndf,
    int direction)
{

  if (multipleMeasurementHandling_ == unweightedClosestToReference) {
    Exception exc("KalmanFitter::fitTrack ==> cannot use unweightedClosestToReference as multiple measurement handling.",__LINE__,__FILE__);
    exc.setFatal();
    throw exc;
  }

  chi2 = 0;
  ndf = -1. * rep->getDim();

#ifdef DEBUG
  std::cout << tr->getNumPoints() << " TrackPoints in this track." << std::endl;
#endif
  for (size_t i = 0; i < tr->getNumPoints(); ++i) {
    TrackPoint *tp = 0;
    KalmanFitterInfo* fi;
    assert(direction == +1 || direction == -1);
    if (direction == +1) {
      tp = tr->getPoint(i);
    }
    else {
      tp = tr->getPoint(-i-1);
    }

    if (! tp->hasFitterInfo(rep)) {
      fi = new KalmanFitterInfo(tp, rep);
      tp->setFitterInfo(fi);
    }
    else
      fi = static_cast<KalmanFitterInfo*>(tp->getFitterInfo(rep));

    #ifdef DEBUG
    std::cout << " process TrackPoint nr. " << i << "\n";
    #endif
    processTrackPoint(tr, tp, fi, rep, chi2, ndf, direction);

  }
}


void KalmanFitter::processTrack(Track* tr, const AbsTrackRep* rep, bool resortHits)
{

  if (tr->getFitStatus(rep) != NULL && tr->getFitStatus(rep)->isTrackPruned()) {
    Exception exc("KalmanFitter::processTrack: Cannot process pruned track!", __LINE__,__FILE__);
    throw exc;
  }

  currentState_.reset(new MeasuredStateOnPlane(rep));
  TVectorD seed(tr->getStateSeed());
  TMatrixDSym cov(6);
  rep->setPosMomCov(*currentState_, seed, cov);

  currentState_->getCov().UnitMatrix();
  //currentState_->getCov() *= blowUpFactor_; // FIXME find good start values

#ifdef DEBUG
  double oldChi2FW = 1e6;
  double oldChi2BW = 1e6;
  double oldPvalFW = 0.;
#endif
  double oldPvalBW = 0.;
  double chi2FW(0), ndfFW(0);
  double chi2BW(0), ndfBW(0);


  KalmanFitStatus* status = new KalmanFitStatus();
  tr->setFitStatus(status, rep);

  unsigned int nIt = 0;
  for(;;) {
    try {
      #ifdef DEBUG
      std::cout << "\033[1;21mstate pre" << std::endl;
      currentState_->Print();
      std::cout << "\033[0mfitting" << std::endl;
      #endif
      fitTrack(tr, rep, chi2FW, ndfFW, +1);
      #ifdef DEBUG
      std::cout << "\033[1;21mstate post forward" << std::endl;
      currentState_->Print();
      std::cout << "\033[0m";
      #endif

      // Backwards iteration:
      currentState_->blowUpCov(blowUpFactor_);  // blow up cov
      fitTrack(tr, rep, chi2BW, ndfBW, -1);
      #ifdef DEBUG
      std::cout << "\033[1;21mstate post backward" << std::endl;
      currentState_->Print();
      std::cout << "\033[0m";

      std::cout << "old chi2s: " << oldChi2BW << ", " << oldChi2FW
          << " new chi2s: " << chi2BW << ", " << chi2FW << std::endl;
      #endif

      ++nIt;

      double PvalBW = ROOT::Math::chisquared_cdf_c(chi2BW, ndfBW);
      #ifdef DEBUG
      double PvalFW = ROOT::Math::chisquared_cdf_c(chi2FW, ndfFW);
      #endif
      // See if p-value only changed little.  If the initial
      // parameters are very far off, initial chi^2 and the chi^2
      // after the first iteration will be both very close to zero, so
      // we need to force at least two iterations here.  Convergence
      // doesn't make much sense before running twice anyway.
      if (nIt > 1 && fabs(oldPvalBW - PvalBW) < deltaPval_)  {
        // Finished
        status->setIsFitConverged();
        break;
      }
      else {
        oldPvalBW = PvalBW;
        #ifdef DEBUG
        oldChi2BW = chi2BW;
        oldChi2FW = chi2FW;
        oldPvalFW = PvalFW;
        #endif
        currentState_->blowUpCov(blowUpFactor_);  // blow up cov
      }

      if (nIt >= maxIterations_) {
        break;
      }
    }
    catch(Exception& e) {
      std::cerr << e.what();
      status->setIsFitted(false);
      status->setIsFitConverged(false);
      return;
    }
  }

  status->setIsFitted();
  status->setCharge(rep->getCharge(*static_cast<KalmanFitterInfo*>(tr->getPointWithMeasurement(0)->getFitterInfo(rep))->getBackwardUpdate()));
  status->setNumIterations(nIt);
  status->setForwardChiSqu(chi2FW);
  status->setBackwardChiSqu(chi2BW);
  status->setForwardNdf(std::max(0., ndfFW));
  status->setBackwardNdf(std::max(0., ndfBW));
}


void
KalmanFitter::processTrackPoint(Track* tr, TrackPoint* tp, KalmanFitterInfo* fi,
    const AbsTrackRep* rep, double& chi2, double& ndf, int direction)
{
  if (!tp->hasRawMeasurements())
    return;

  // Extrapolate to TrackPoint.
  MeasuredStateOnPlane* state = new MeasuredStateOnPlane(*currentState_);
  //state.Print();

  // construct measurementsOnPlane if it has not yet been done
  if (fi->getNumMeasurements() == 0) {
    std::vector< genfit::AbsMeasurement* > rawMeasurements =  tp->getRawMeasurements();
    // construct plane with first measurement
    SharedPlanePtr plane = rawMeasurements[0]->constructPlane(currentState_.get());
    for (std::vector< genfit::AbsMeasurement* >::iterator it = rawMeasurements.begin(); it != rawMeasurements.end(); ++it) {
      fi->setMeasurementsOnPlane((*it)->constructMeasurementsOnPlane(rep, plane));
    }
  }
  const SharedPlanePtr plane = fi->getPlane();

#ifdef DEBUG
  std::cout << "its plane is at R = " << plane->getO().Perp()
	    << " with normal pointing along (" << plane->getNormal().X() << ", " << plane->getNormal().Y() << ", " << plane->getNormal().Z() << ")" << std::endl;
#endif

  //state.Print();

#ifdef DEBUG
  double extLen = rep->extrapolateToPlane(state, plane);
  std::cout << "extrapolated by " << extLen << std::endl;
#else
  rep->extrapolateToPlane(*state, plane);
#endif
  //std::cout << "after extrap: " << std::endl;
  //state.Print();

  // unique_ptr takes care of disposing of the old prediction, takes ownership of state.
  assert(direction == -1 || direction == +1);
  fi->setPrediction(state, direction);

  TVectorD stateVector(state->getState());
  TMatrixDSym cov(state->getCov());
  const MeasurementOnPlane& mOnPlane = getMeasurement(fi, direction);
  const TVectorD& measurement(mOnPlane.getState());
  const TMatrixDSym& V(mOnPlane.getCov());
  const TMatrixD& H(mOnPlane.getHMatrix());
#ifdef DEBUG
  std::cout << "State prediction: "; stateVector.Print();
  std::cout << "Cov prediction: "; state->getCov().Print();
  //cov.Print();
  //measurement.Print();
#endif

  TVectorD res(measurement - H*stateVector);
#ifdef DEBUG
  std::cout << "Residual = (" << res(0);
  if (res.GetNrows() > 1)
    std::cout << ", " << res(1);
  std::cout << ")" << std::endl;
#endif
  // If hit, do Kalman algebra.

  // calculate kalman gain ------------------------------
  // calculate covsum (V + HCH^T)
  TMatrixDSym HcovHt(cov);
  HcovHt.Similarity(H);

  TMatrixDSym covSumInv(V + HcovHt);
  tools::invertMatrix(covSumInv);

  TMatrixD CHt(cov, TMatrixD::kMultTranspose, H);
  TVectorD update = TMatrixD(CHt, TMatrixD::kMult, covSumInv) * res;

#ifdef DEBUG
  //std::cout << "STATUS:" << std::endl;
  //stateVector.Print();
  std::cout << "Update: "; update.Print();
  //cov.Print();
#endif

  stateVector += update;
  covSumInv.Similarity(CHt); // with (C H^T)^T = H C^T = H C  (C is symmetric)
  cov -= covSumInv;
#ifdef DEBUG
  std::cout << "updated state: "; stateVector.Print();
  std::cout << "updated cov: "; cov.Print();
#endif

  TVectorD resNew(measurement - H*stateVector);
#ifdef DEBUG
  std::cout << "Residual New = (" << resNew(0);

  if (resNew.GetNrows() > 1)
    std::cout << ", " << resNew(1);
  std::cout << ")" << std::endl;
#endif

  currentState_->setStateCovPlane(stateVector, cov, plane);
  currentState_->setAuxInfo(state->getAuxInfo());

  /*TDecompChol dec(cov);
  TMatrixDSym mist(cov);
  bool status = dec.Invert(mist);
  if (!status) {
#ifdef DEBUG
    std::cout << "new cov not pos. def." << std::endl;
#endif
  }*/

  // Calculate chi²
  TMatrixDSym HCHt(cov);
  HCHt.Similarity(H);
  HCHt -= V;
  HCHt *= -1;

  tools::invertMatrix(HCHt);

  double chi2inc = HCHt.Similarity(resNew);
  chi2 += chi2inc;

  double ndfInc = measurement.GetNrows();
  ndf += ndfInc;
#ifdef DEBUG
  std::cout << "chi² increment = " << chi2inc << std::endl;
#endif

  // set update
  KalmanFittedStateOnPlane* updatedSOP = new KalmanFittedStateOnPlane(*currentState_, chi2inc, ndfInc);
  fi->setUpdate(updatedSOP, direction);
}


// Modified from auto-generated streamer to deal with scoped_ptr correctly.
void KalmanFitter::Streamer(TBuffer &R__b)
{
   // Stream an object of class genfit::KalmanFitter.

   //This works around a msvc bug and should be harmless on other platforms
   typedef ::genfit::KalmanFitter thisClass;
   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      //This works around a msvc bug and should be harmless on other platforms
      typedef genfit::AbsKalmanFitter baseClass0;
      baseClass0::Streamer(R__b);
      MeasuredStateOnPlane *p = 0;
      R__b >> p;
      currentState_.reset(p);
      R__b.CheckByteCount(R__s, R__c, thisClass::IsA());
   } else {
      R__c = R__b.WriteVersion(thisClass::IsA(), kTRUE);
      //This works around a msvc bug and should be harmless on other platforms
      typedef genfit::AbsKalmanFitter baseClass0;
      baseClass0::Streamer(R__b);
      R__b << currentState_.get();
      R__b.SetByteCount(R__c, kTRUE);
   }
}
