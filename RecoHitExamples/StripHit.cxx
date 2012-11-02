#include "StripHit.h"

#include "RKTrackRep.h"
#include "GFDetPlane.h"
#include "TRandom.h"
#include "TMath.h"

#include"math.h"

ClassImp(StripHit)


StripHit::~StripHit()
{}

StripHit::StripHit()
  : GFAbsPlanarHit(NparHitRep)
{}

StripHit::StripHit(const TVector3& point, const TVector3& norm,
		               const TVector3& u, double res,
		               bool smear)
  : GFAbsPlanarHit(NparHitRep){

  fHitCov(0,0) = res*res;

  assert(fabs(norm*u)<1.E-5);
  TVector3 v = u.Cross(norm);
  GFDetPlane plane(point,u,v);

  if(smear) fHitCoord(0,0) = gRandom->Gaus(0,res);
  else fHitCoord(0,0) = 0;

  setDetPlane(plane);
}


GFAbsRecoHit* 
StripHit::clone(){
  return new StripHit(*this);
}


TMatrixT<double>
StripHit::getHMatrix(const GFAbsTrackRep* stateVector) {
  if ( (dynamic_cast<const RKTrackRep*>(stateVector) != NULL)){
    TMatrixT<double> HMatrix(1,5);

    HMatrix(0,0) = 0.;
    HMatrix(0,1) = 0.;
    HMatrix(0,2) = 0.;
    HMatrix(0,3) = 1.;
    HMatrix(0,4) = 0.;

    return HMatrix;
  }
  else {
    std::cerr << "StripHit can only handle state vectors of type RKtrackRep -> abort" << std::endl;
    throw;
  }
 
}


