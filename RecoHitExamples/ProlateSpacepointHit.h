#ifndef ProlateSpacepointHit_HH
#define ProlateSpacepointHit_HH

#include "GFAbsProlateSpacepointHit.h"

class ProlateSpacepointHit : public GFAbsProlateSpacepointHit {
public:

  ProlateSpacepointHit();
  ProlateSpacepointHit(const TVector3& pos, const TVector3& wireDir,
                          double resPerp, double resWire, bool smear=false);

  virtual ~ProlateSpacepointHit();

  virtual GFAbsRecoHit* clone();
  
  virtual TMatrixT<double> getHMatrix(const GFAbsTrackRep* stateVector);

 public:
  ClassDef(ProlateSpacepointHit,1)

};

#endif

//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
