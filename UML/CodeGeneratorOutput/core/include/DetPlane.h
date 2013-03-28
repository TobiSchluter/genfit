#ifndef genfit_DetPlane_h
#define genfit_DetPlane_h


#include <TVector3.h>

#include "AbsFinitePlane.h"


namespace genfit {

class DetPlane {


 private:
  TVector3 o_;
  TVector3 u_;
  TVector3 v_;

 public:


  AbsFinitePlane* finitePlane_;

};

} /* End of namespace genfit */

#endif // genfit_DetPlane_h
