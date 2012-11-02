#ifndef PIXHIT_HH
#define PIXHIT_HH

#include "GFAbsPlanarHit.h"

class PixHit : public GFAbsPlanarHit {
public:

  PixHit();
  PixHit(const TVector3& point,const TVector3& norm, const TVector3& u, double res, bool smear=false);

  virtual ~PixHit();

  virtual GFAbsRecoHit* clone();
  
  virtual TMatrixT<double> getHMatrix(const GFAbsTrackRep* stateVector);

private:
  static const int NparHitRep = 2;

public:
  ClassDef(PixHit,1)

};

#endif

//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
