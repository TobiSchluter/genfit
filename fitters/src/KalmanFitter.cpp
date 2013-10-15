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
#include <TDecompChol.h>
#include <TMatrixDSymEigen.h>
#include <algorithm>

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

  if (debugLvl_ > 0) {
    std::cout << tr->getNumPoints() << " TrackPoints in this track." << std::endl;
  }
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

    if (debugLvl_ > 0) {
    std::cout << " process TrackPoint nr. " << i << " (" << tp << ")\n";
    }
    processTrackPoint(tr, tp, fi, rep, chi2, ndf, direction);

  }
}


void KalmanFitter::processTrackWithRep(Track* tr, const AbsTrackRep* rep, bool resortHits)
{

  if (tr->getFitStatus(rep) != NULL && tr->getFitStatus(rep)->isTrackPruned()) {
    Exception exc("KalmanFitter::processTrack: Cannot process pruned track!", __LINE__,__FILE__);
    throw exc;
  }

  currentState_.reset(new MeasuredStateOnPlane(rep));

  TrackPoint* trackPoint = tr->getPointWithMeasurement(0);

  if (rep != tr->getCardinalRep() &&
      trackPoint->hasFitterInfo(tr->getCardinalRep()) &&
      dynamic_cast<KalmanFitterInfo*>(trackPoint->getFitterInfo(tr->getCardinalRep())) != NULL &&
      static_cast<KalmanFitterInfo*>(trackPoint->getFitterInfo(tr->getCardinalRep()))->hasPredictionsAndUpdates() ) {
    if (debugLvl_ > 0)
      std::cout << "take smoothed state of cardinal rep fit as seed \n";
    TVector3 pos, mom;
    TMatrixDSym cov;
    const MeasuredStateOnPlane& fittedState = static_cast<KalmanFitterInfo*>(trackPoint->getFitterInfo(tr->getCardinalRep()))->getFittedState(true);
    tr->getCardinalRep()->getPosMomCov(fittedState, pos, mom, cov);
    rep->setPosMomCov(*currentState_, pos, mom, cov);
    rep->setQop(*currentState_, tr->getCardinalRep()->getQop(fittedState));
  }
  else {
    rep->setPosMomCov(*currentState_, tr->getStateSeed(), tr->getCovSeed());
  }

  // Only after we have linearly propagated the error into the TrackRep can we
  // blow up the error in good faith.
  currentState_->blowUpCov(blowUpFactor_);

  double oldChi2FW(1.e6);
  double oldChi2BW(1.e6);
  double oldPvalFW(0.);

  double oldPvalBW = 0.;
  double chi2FW(0), ndfFW(0);
  double chi2BW(0), ndfBW(0);


  KalmanFitStatus* status = new KalmanFitStatus();
  tr->setFitStatus(status, rep);

  unsigned int nIt = 0;
  for(;;) {
    try {
      if (debugLvl_ > 0) {
        std::cout << "\033[1;21mstate pre" << std::endl;
        currentState_->Print();
        std::cout << "\033[0mfitting" << std::endl;
      }

      fitTrack(tr, rep, chi2FW, ndfFW, +1);

      if (debugLvl_ > 0) {
        std::cout << "\033[1;21mstate post forward" << std::endl;
        currentState_->Print();
        std::cout << "\033[0m";
      }

      // Backwards iteration:
      currentState_->blowUpCov(blowUpFactor_);  // blow up cov

      fitTrack(tr, rep, chi2BW, ndfBW, -1);

      if (debugLvl_ > 0) {
        std::cout << "\033[1;21mstate post backward" << std::endl;
        currentState_->Print();
        std::cout << "\033[0m";

        std::cout << "old chi2s: " << oldChi2BW << ", " << oldChi2FW
		  << " new chi2s: " << chi2BW << ", " << chi2FW
		  << " oldPvals " << oldPvalFW << ", " << oldPvalBW << std::endl;
      }

      ++nIt;

      double PvalBW = ROOT::Math::chisquared_cdf_c(chi2BW, ndfBW);
      double PvalFW = (debugLvl_ > 0) ? ROOT::Math::chisquared_cdf_c(chi2FW, ndfFW) : 0; // Don't calculate if not debugging as this function potentially takes a lot of time.
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
        if (debugLvl_ > 0) {
          oldChi2BW = chi2BW;
          oldChi2FW = chi2FW;
          oldPvalFW = PvalFW;
        }
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
  status->setForwardChi2(chi2FW);
  status->setBackwardChi2(chi2BW);
  status->setForwardNdf(std::max(0., ndfFW));
  status->setBackwardNdf(std::max(0., ndfBW));
}


void
KalmanFitter::processTrackPoint(Track* tr, TrackPoint* tp, KalmanFitterInfo* fi,
    const AbsTrackRep* rep, double& chi2, double& ndf, int direction)
{
  assert(direction == -1 || direction == +1);

  if (!tp->hasRawMeasurements())
    return;

  // construct measurementsOnPlane if forward fit
  if (direction == 1) {
    // delete outdated stuff
    fi->deleteForwardInfo();
    fi->deleteBackwardInfo();
    fi->deleteMeasurementInfo();

    // construct new MeasurementsOnPlane
    std::vector<double> oldWeights = fi->getWeights();
    bool oldWeightsFixed = fi->areWeightsFixed();
    const std::vector< genfit::AbsMeasurement* >& rawMeasurements =  tp->getRawMeasurements();
    // construct plane with first measurement
    SharedPlanePtr plane = rawMeasurements[0]->constructPlane(*currentState_);
    for (std::vector< genfit::AbsMeasurement* >::const_iterator it = rawMeasurements.begin(); it != rawMeasurements.end(); ++it) {
      fi->addMeasurementsOnPlane((*it)->constructMeasurementsOnPlane(rep, plane));
    }
    if (oldWeights.size() == fi->getNumMeasurements()) {
      fi->setWeights(oldWeights);
      fi->fixWeights(oldWeightsFixed);
    }
    assert(fi->getPlane() == plane);
    assert(fi->checkConsistency());
  }
  const SharedPlanePtr plane = fi->getPlane();

  if (debugLvl_ > 0) {
    std::cout << "its plane is at R = " << plane->getO().Perp()
        << " with normal pointing along (" << plane->getNormal().X() << ", " << plane->getNormal().Y() << ", " << plane->getNormal().Z() << ")" << std::endl;
  }


  // Extrapolate
  MeasuredStateOnPlane* state = new MeasuredStateOnPlane(*currentState_);
  double extLen = rep->extrapolateToPlane(*state, plane);

  if (debugLvl_ > 0) {
    std::cout << "extrapolated by " << extLen << std::endl;
  }
  //std::cout << "after extrap: " << std::endl;
  //state.Print();

  // unique_ptr takes care of disposing of the old prediction, takes ownership of state.
  fi->setPrediction(state, direction);

  TVectorD stateVector(state->getState());
  TMatrixDSym cov(state->getCov());

  // update(s)
  double chi2inc = 0;
  double ndfInc = 0;
  const std::vector<MeasurementOnPlane *> measurements = getMeasurements(fi, direction);
  for (std::vector<MeasurementOnPlane *>::const_iterator it = measurements.begin(); it != measurements.end(); ++it) {
    const MeasurementOnPlane& mOnPlane = **it;

    if (!canIgnoreWeights() && mOnPlane.getWeight() <= 1.01E-10) {
      if (debugLvl_ > 0) {
        std::cout << "Weight of measurement is almost 0, continue ... \n";
      }
      continue;
    }

    const TVectorD& measurement(mOnPlane.getState());
    const AbsHMatrix* H(mOnPlane.getHMatrix());
    // (weighted) cov
    const TMatrixDSym& V((!canIgnoreWeights() && mOnPlane.getWeight() < 0.99999) ?
                          1./mOnPlane.getWeight() * mOnPlane.getCov() :
                          mOnPlane.getCov());
    if (debugLvl_ > 1) {
      std::cout << "State prediction: "; stateVector.Print();
      std::cout << "Cov prediction: "; state->getCov().Print();
      //cov.Print();
      //measurement.Print();
    }

    TVectorD res(measurement - H->Hv(stateVector));
    if (debugLvl_ > 1) {
      std::cout << "Residual = (" << res(0);
      if (res.GetNrows() > 1)
        std::cout << ", " << res(1);
      std::cout << ")" << std::endl;
    }
    // If hit, do Kalman algebra.

    if (squareRootFormalism_)
      {
        // Square Root formalism, Anderson & Moore, 6.5.12 gives the
        // formalism for combined update and prediction in that
        // sequence.  We need the opposite sequence.  Therefore, this is
        // my own invention.
        TMatrixD F;
        TMatrixDSym noise;
        TVectorD deltaState;
        rep->getForwardJacobianAndNoise(F, noise, deltaState);
        const TMatrixDSym& oldCov(currentState_->getCov());

        // Square Roots:
        //
        // The noise matrix is only positive semi-definite, for instance
        // no noise on q/p.  A Cholesky decomposition is therefore not
        // possible.  Hence we pull out the eigenvalues, viz noise = z d
        // z^T (with z orthogonal, d diagonal) and then construct the
        // square root according to Q^T = sqrt(d) z where sqrt is
        // applied element-wise and negative eigenvalues are forced to
        // zero.  We then reduce the matrix such that the zero
        // eigenvalues don't appear in the remainder of the calculation.
        //
        // FIXME
        // Careful, this is not tested with actual noise.  Do I need the
        // transpose???
        TMatrixDSymEigen eig(noise);
        TMatrixD Q(eig.GetEigenVectors());
        const TVectorD& evs(eig.GetEigenValues());
        // Multiplication with a diagonal matrix ... eigenvalues are
        // sorted in descending order.
        int iCol = 0;
        for (; iCol < Q.GetNcols(); ++iCol) {
          double ev = evs(iCol) > 0 ? sqrt(evs(iCol)) : 0;
          //FIXME, see below we do no resizing as of now
          if (ev == 0)
            break;
          for (int j = 0; j < Q.GetNrows(); ++j)
            Q(j,iCol) *= ev;
        }
        if (iCol < Q.GetNrows()) {
          // Hit zero eigenvalue, resize matrix ...
          Q.ResizeTo(iCol, Q.GetNcols());
        }
        // This gives the original matrix:
        // TMatrixDSym(TMatrixDSym::kAtA,TMatrixD(TMatrixD::kTransposed,Q)).Print();
        // i.e., we now have the transposed we need.

        //std::cout << "cov root" << std::endl;
        TDecompChol oldCovDecomp(oldCov);
        oldCovDecomp.Decompose();
        const TMatrixD& S(oldCovDecomp.GetU()); // this actually the transposed we want.

        TDecompChol decompR(V);
        decompR.Decompose();

        TMatrixD pre(S.GetNrows() + Q.GetNrows() + V.GetNrows(),
         S.GetNcols() + V.GetNcols());
        TMatrixD SFt(S, TMatrixD::kMultTranspose, F);
        pre.SetSub(                        0,0,decompR.GetU()); /* upper right block is zero */
        pre.SetSub(             V.GetNrows(),0,H->MHt(SFt)); pre.SetSub(V.GetNrows(),             V.GetNcols(),SFt);
        if (Q.GetNrows() > 0) { // needed to suppress warnings when inserting an empty Q
          pre.SetSub(S.GetNrows()+V.GetNrows(),0,H->MHt(Q));   pre.SetSub(S.GetNrows()+V.GetNrows(),V.GetNcols(),Q);
        }


        tools::QR(pre);
        const TMatrixD& r = pre;

        TMatrixD R(r.GetSub(0, V.GetNrows()-1, 0, V.GetNcols()-1));
        //R.Print();
        TMatrixD K(TMatrixD::kTransposed, r.GetSub(0, V.GetNrows()-1, V.GetNcols(), pre.GetNcols()-1));
        //K.Print();
        //(K*R).Print();
        TMatrixD Snew(r.GetSub(V.GetNrows(), V.GetNrows() + S.GetNrows() - 1,
             V.GetNcols(), pre.GetNcols()-1));
        //Snew.Print();
        // No need for matrix inversion.
        TVectorD update(res);
        tools::transposedForwardSubstitution(R, update);
        update *= K;
        stateVector += update;
        cov = TMatrixDSym(TMatrixDSym::kAtA, Snew); // Note that this is transposed in just the right way.
      }
    else
      {
        // calculate kalman gain ------------------------------
        // calculate covsum (V + HCH^T) and invert
        TMatrixDSym covSumInv(cov);
        H->HMHt(covSumInv);
        covSumInv += V;
        tools::invertMatrix(covSumInv);

        TMatrixD CHt(H->MHt(cov));
        TVectorD update(TMatrixD(CHt, TMatrixD::kMult, covSumInv) * res);
        //TMatrixD(CHt, TMatrixD::kMult, covSumInv).Print();

        if (debugLvl_ > 1) {
          //std::cout << "STATUS:" << std::endl;
          //stateVector.Print();
          std::cout << "Update: "; update.Print();
          //cov.Print();
        }

        stateVector += update;
        covSumInv.Similarity(CHt); // with (C H^T)^T = H C^T = H C  (C is symmetric)
        cov -= covSumInv;
      }

    if (debugLvl_ > 1) {
      std::cout << "updated state: "; stateVector.Print();
      std::cout << "updated cov: "; cov.Print();
    }

    TVectorD resNew(measurement - H->Hv(stateVector));
    if (debugLvl_ > 1) {
      std::cout << "Residual New = (" << resNew(0);

      if (resNew.GetNrows() > 1)
        std::cout << ", " << resNew(1);
      std::cout << ")" << std::endl;
    }

    /*TDecompChol dec(cov);
    TMatrixDSym mist(cov);
    bool status = dec.Invert(mist);
    if (!status) {
      if (debugLvl_ > 0) {
          std::cout << "new cov not pos. def." << std::endl;
      }
    }*/

    // Calculate chi²
    TMatrixDSym HCHt(cov);
    H->HMHt(HCHt);
    HCHt -= V;
    HCHt *= -1;

    tools::invertMatrix(HCHt);

    chi2inc += HCHt.Similarity(resNew);

    if (!canIgnoreWeights()) {
      ndfInc += mOnPlane.getWeight() * mOnPlane.getState().GetNrows();
    }
    else
      ndfInc += mOnPlane.getState().GetNrows();

    if (debugLvl_ > 0) {
      std::cout << "chi² increment = " << chi2inc << std::endl;
    }
  } // end loop over measurements

  currentState_->setStateCovPlane(stateVector, cov, plane);
  currentState_->setAuxInfo(state->getAuxInfo());

  chi2 += chi2inc;
  ndf += ndfInc;

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
