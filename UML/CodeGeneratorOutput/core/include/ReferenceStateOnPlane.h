#ifndef genfit_ReferenceStateOnPlane_h
#define genfit_ReferenceStateOnPlane_h

#include "StateOnPlane.h"


namespace genfit {


  /** 
   *  Transport matrices describe transport TO that plane.
   */
class ReferenceStateOnPlane : public StateOnPlane {

 public:
  double forwardSegmentLength_;
  double backwardSegmentLength_;

 protected:
  TMatrixD forwardTransportMatrix;
  TMatrixD backwardTransportMatrix;
  TMatrixDSym forwardNoiseMatrix;
  TMatrixDSym backwardNoiseMatrix;

};

} /* End of namespace genfit */

#endif // genfit_ReferenceStateOnPlane_h
