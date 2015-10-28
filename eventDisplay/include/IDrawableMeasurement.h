#ifndef genfit_IDrawableMeasurement_h
#define genfit_IDrawableMeasurement_h

class TEveElementList;
class TEveBox;
class TVector3;

namespace genfit {

class MeasuredStateOnPlane;

namespace display {

/*
 *Interface class for hit drawing in TEve.
 */
class IDrawableMeasurement {
public:
  virtual void drawMeasurement(TEveElementList* list, const MeasuredStateOnPlane& fittedState) const = 0;
  virtual void drawDetector(TEveElementList* list) const = 0;

  virtual ~IDrawableMeasurement() {}
};

/**
 * @brief Utility function for drawing detector planes.  Returns a TEveBox
 * of thickness depth centered at O with corners displaced by ud and
 * uv in the directions given by U and V.
 *
 * Implementation in EventDisplay.cc.
 */
TEveBox* boxCreator(const TVector3& O, TVector3 U, TVector3 V, float ud, float vd, float depth);

}
}

#endif
