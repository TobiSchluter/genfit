#ifndef genfit_RectangularFinitePlane_h
#define genfit_RectangularFinitePlane_h

#include "AbsFinitePlane.h"


namespace genfit {

class RectangularFinitePlane : public AbsFinitePlane {

 public:

  //override inActive & Print methods
  bool inActive(const double& u, const double& v) const;
  void Print(const Option_t* = "") const;

  //! give dimensions of finite rectangle: u1,u2,v1,v2
  RectangularFinitePlane(const double&, const double&, const double&, const double&);
  RectangularFinitePlane();

  virtual ~RectangularFinitePlane();

  AbsFinitePlane* clone() const {
    return new RectangularFinitePlane(*this);
  }

 private:

  double uMin_, uMax_, vMin_, vMax_;


  ClassDef(RectangularFinitePlane,1)

};

} /* End of namespace genfit */

#endif // genfit_RectangularFinitePlane_h
