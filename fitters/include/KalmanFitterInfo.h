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

#ifndef genfit_KalmanFitterInfo_h
#define genfit_KalmanFitterInfo_h

#include <vector>

#ifndef __CINT__
#include "boost/scoped_ptr.hpp"
#endif

#include "AbsFitterInfo.h"
#include "KalmanFittedStateOnPlane.h"
#include "MeasuredStateOnPlane.h"
#include "MeasurementOnPlane.h"
#include "ReferenceStateOnPlane.h"
#include "StateOnPlane.h"


namespace genfit {


  /** 
   *  This class collects all information needed and produced by a Kalman filter or DAF and is specific to one #GFAbsTrackRep of the #GFTrack.
   */
class KalmanFitterInfo : public AbsFitterInfo {

 public:

  KalmanFitterInfo();
  KalmanFitterInfo(const TrackPoint* trackPoint, const AbsTrackRep* rep);
  ~KalmanFitterInfo();

  virtual KalmanFitterInfo* clone() const;

  ReferenceStateOnPlane* getReferenceState() const {return referenceState_.get();}
  MeasuredStateOnPlane* getForwardPrediction() const {return forwardPrediction_.get();}
  MeasuredStateOnPlane* getBackwardPrediction() const {return backwardPrediction_.get();}
  MeasuredStateOnPlane* getPrediction(int direction) const {if (direction >=0) return forwardPrediction_.get(); return backwardPrediction_.get();}
  KalmanFittedStateOnPlane* getForwardUpdate() const {return forwardUpdate_.get();}
  KalmanFittedStateOnPlane* getBackwardUpdate() const {return backwardUpdate_.get();}
  KalmanFittedStateOnPlane* getUpdate(int direction) const {if (direction >=0) return forwardUpdate_.get(); return backwardUpdate_.get();}
  std::vector< genfit::MeasurementOnPlane* > getMeasurementsOnPlane() const;
  MeasurementOnPlane* getMeasurementOnPlane(int i = 0) const {if (i<0) i += measurementsOnPlane_.size(); return measurementsOnPlane_.at(i);}
  /** Get weighted mean of all measurements. */
  MeasurementOnPlane getAvgWeightedMeasurementOnPlane() const;
  /** Get measurements which is closest to state. */
  const MeasurementOnPlane* getClosestMeasurementOnPlane(const StateOnPlane*) const;
  unsigned int getNumMeasurements() const {return measurementsOnPlane_.size();}
  /** Get weights of measurements */
  std::vector<double> getWeights() const;
  SharedPlanePtr getPlane() const;
  /** Get unbiased (default) or biased smoothed state */
  MeasuredStateOnPlane getFittedState(bool biased = false) const;
  /** Get unbiased (default) or biased residual from ith measurement */
  MeasurementOnPlane getResidual(bool biased = false, unsigned int iMeasurement = 0) const; // also calculates covariance of the residual

  bool hasReferenceState() const {return (referenceState_.get() != NULL);}
  bool hasForwardPrediction() const {return (forwardPrediction_.get()  != NULL);}
  bool hasBackwardPrediction() const {return (backwardPrediction_.get() != NULL);}
  bool hasPrediction(int direction) const {if (direction >=0) return hasForwardPrediction(); return hasBackwardPrediction();}
  bool hasForwardUpdate() const {return (forwardUpdate_.get() != NULL);}
  bool hasBackwardUpdate() const {return (backwardUpdate_.get() != NULL);}
  bool hasUpdate(int direction) const {if (direction >=0) return hasForwardUpdate(); return hasBackwardUpdate();}

  void setReferenceState(ReferenceStateOnPlane* referenceState) {referenceState_.reset(referenceState);}
  void setForwardPrediction(MeasuredStateOnPlane* forwardPrediction) {forwardPrediction_.reset(forwardPrediction);}
  void setBackwardPrediction(MeasuredStateOnPlane* backwardPrediction) {backwardPrediction_.reset(backwardPrediction);}
  void setPrediction(MeasuredStateOnPlane* prediction, int direction)  {if (direction >=0) setForwardPrediction(prediction); else setBackwardPrediction(prediction);}
  void setForwardUpdate(KalmanFittedStateOnPlane* forwardUpdate) {forwardUpdate_.reset(forwardUpdate);}
  void setBackwardUpdate(KalmanFittedStateOnPlane* backwardUpdate) {backwardUpdate_.reset(backwardUpdate);}
  void setUpdate(KalmanFittedStateOnPlane* update, int direction)  {if (direction >=0) setForwardUpdate(update); else setBackwardUpdate(update);}
  void setMeasurementsOnPlane(const std::vector< genfit::MeasurementOnPlane* >& measurementsOnPlane);
  void addMeasurementOnPlane(MeasurementOnPlane* measurementOnPlane) { measurementsOnPlane_.push_back(measurementOnPlane); }
  void addMeasurementsOnPlane(const std::vector< genfit::MeasurementOnPlane* >& measurementsOnPlane);
  /** Set weights of measurements */
  void setWeights(const std::vector<double>&);
  void setRep(const AbsTrackRep* rep);

  void deleteForwardInfo();
  void deleteBackwardInfo();
  void deleteReferenceInfo();
  void deleteMeasurementInfo();

  virtual void Print(const Option_t* = "") const;

  virtual bool checkConsistency() const;

 private:

  MeasuredStateOnPlane calcAverageState(const MeasuredStateOnPlane* forwardState, const MeasuredStateOnPlane* backwardState) const;

#ifndef __CINT__
  boost::scoped_ptr<ReferenceStateOnPlane> referenceState_; // Ownership
  boost::scoped_ptr<MeasuredStateOnPlane> forwardPrediction_; // Ownership
  boost::scoped_ptr<KalmanFittedStateOnPlane> forwardUpdate_; // Ownership
  boost::scoped_ptr<MeasuredStateOnPlane> backwardPrediction_; // Ownership
  boost::scoped_ptr<KalmanFittedStateOnPlane> backwardUpdate_; // Ownership
#else
  class ReferenceStateOnPlane* referenceState_;
  class MeasuredStateOnPlane* forwardPrediction_;
  class KalmanFittedStateOnPlane* forwardUpdate_;
  class MeasuredStateOnPlane* backwardPrediction_;
  class KalmanFittedStateOnPlane* backwardUpdate_;
#endif

 //> TODO ! ptr implement: to the special ownership version
  /* class owned_pointer_vector : private std::vector<MeasuredStateOnPlane*> {
   public: 
    ~owned_pointer_vector() { for (size_t i = 0; i < this->size(); ++i)
                         delete this[i]; }
    size_t size() const { return this->size(); };
    void push_back(MeasuredStateOnPlane* measuredState) { this->push_back(measuredState); };
    const  MeasuredStateOnPlane* at(size_t i)  const { return this->at(i); }; 
	//owned_pointer_vector::iterator erase(owned_pointer_vector::iterator position) ;
	//owned_pointer_vector::iterator erase(owned_pointer_vector::iterator first, owned_pointer_vector::iterator last);
};
	*/



  /** 
   *  Number of measurements must be equal to size of #fRawMeasurements in #GFTrackPoint.
   * @element-type MeasurementOnPlane
   */
  std::vector<MeasurementOnPlane*> measurementsOnPlane_; // Ownership



  ClassDef(KalmanFitterInfo,1)

};

} /* End of namespace genfit */
/** @} */

#endif // genfit_KalmanFitterInfo_h
