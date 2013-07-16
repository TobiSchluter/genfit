/* Copyright 2008-2010, Technische Universitaet Muenchen,
   Authors: Christian Hoeppner & Sebastian Neubert & Johannes Rauch

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

/** @addtogroup genfit
 * @{
 */

#ifndef genfit_KalmanFitStatus_h
#define genfit_KalmanFitStatus_h

#include "FitStatus.h"

namespace genfit {

class KalmanFitStatus : public FitStatus {

 public:

  KalmanFitStatus() :
    FitStatus(), numIterations_(0), fittedWithDaf_(false), fittedWithReferenceTrack_(false),
    fChi2_(0), bChi2_(0), fNdf_(0), bNdf_(0) {;}

  virtual ~KalmanFitStatus() {};

  virtual FitStatus* clone() const {return new KalmanFitStatus(*this);}

  unsigned int getNumIterations() const {return numIterations_;}
  bool isFittedWithDaf() const {return fittedWithDaf_;}
  bool isFittedWithReferenceTrack() const {return fittedWithReferenceTrack_;}
  double getForwardChiSqu() const {return fChi2_;}
  double getBackwardChiSqu() const {return bChi2_;}
  double getForwardNdf() const {return fNdf_;}
  double getBackwardNdf() const {return bNdf_;}

  void setNumIterations(unsigned int numIterations) {numIterations_ = numIterations;}
  void setIsFittedWithDaf(bool fittedWithDaf = true) {fittedWithDaf_ = fittedWithDaf;}
  void setIsFittedWithReferenceTrack(bool fittedWithReferenceTrack = true) {fittedWithReferenceTrack_ = fittedWithReferenceTrack;}
  void setForwardChiSqu(double fChi2) {fChi2_ = fChi2;}
  void setBackwardChiSqu(double bChi2) {bChi2_ = bChi2;}
  void setForwardNdf(double fNdf) {fNdf_ = fNdf;}
  void setBackwardNdf(double bNdf) {bNdf_ = bNdf;}

  void Print(const Option_t* = "") const {
    FitStatus::Print();
    if (fittedWithDaf_) std::cout << " track has been fitted with DAF,";
    if (fittedWithReferenceTrack_) std::cout << " track has been fitted with reference track,";
    if (isFitted_) {
      std::cout << " numIterations = " << numIterations_ << ", ";
      std::cout << "fChi2 = " << fChi2_ << ", ";
      std::cout << "bChi2 = " << bChi2_ << ", ";
      std::cout << "fNdf = " << fNdf_ << ", ";
      std::cout << "bNdf = " << bNdf_ << "\n";
    }
    std::cout << "\n";
  }

 protected:

  unsigned int numIterations_; // number of iterations that have been performed
  bool fittedWithDaf_;
  bool fittedWithReferenceTrack_;

  double fChi2_; // chi^2 of the forward fit
  double bChi2_; // chi^2 of the backward fit
  double fNdf_; // degrees of freedom of the forward fit
  double bNdf_; // degrees of freedom of the backward fit

  ClassDef(KalmanFitStatus, 1)
};

} /* End of namespace genfit */
/** @} */

#endif // genfit_KalmanFitStatus_h
