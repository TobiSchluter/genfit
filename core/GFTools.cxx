
#include <memory>
#include "GFTools.h"
#include <typeinfo>
#include <TDecompChol.h>

TVectorT<double> GFTools::getSmoothedPos(const GFTrack* trk, unsigned int irep, unsigned int ihit) {

	TVectorT<double> smoothed_state;
	TMatrixTSym<double> smoothed_cov;
	TVectorT<double> pos;

	if(GFTools::getSmoothedData(trk, irep, ihit, smoothed_state, smoothed_cov)) {

		TMatrixT<double> H = trk->getHit(ihit)->getHMatrix(trk->getTrackRep(irep));
		TVectorT<double> pos_tmp(H * smoothed_state);
		pos.ResizeTo(pos_tmp);
		pos = pos_tmp;

	}
	return pos;

}


TVector3 GFTools::getSmoothedPosXYZ(const GFTrack* trk, unsigned int irep, unsigned int ihit){

  std::auto_ptr<GFAbsTrackRep> rep(trk->getTrackRep(irep)->clone());

  TVectorT<double> smoothed_state;
  TMatrixTSym<double> smoothed_cov;
  TMatrixT<double> auxInfo;
  GFDetPlane smoothing_plane;

  getSmoothedData(trk, irep, ihit, smoothed_state, smoothed_cov, smoothing_plane, auxInfo);

  if(rep->hasAuxInfo()) {
    rep->setData(smoothed_state, smoothing_plane, &smoothed_cov, &auxInfo);
  } else {
    rep->setData(smoothed_state, smoothing_plane, &smoothed_cov);
  }

  return rep->getPos(smoothing_plane);
}

TVector3 GFTools::getBiasedSmoothedPosXYZ(const GFTrack* trk, unsigned int irep, unsigned int ihit){

  std::auto_ptr<GFAbsTrackRep> rep(trk->getTrackRep(irep)->clone());

  TVectorT<double> smoothed_state;
  TMatrixTSym<double> smoothed_cov;
  TMatrixT<double> auxInfo;
  GFDetPlane smoothing_plane;

  getBiasedSmoothedData(trk, irep, ihit, smoothed_state, smoothed_cov, smoothing_plane, auxInfo);

  if(rep->hasAuxInfo()) {
    rep->setData(smoothed_state, smoothing_plane, &smoothed_cov, &auxInfo);
  } else {
    rep->setData(smoothed_state, smoothing_plane, &smoothed_cov);
  }

  return rep->getPos(smoothing_plane);
}


TVector3 GFTools::getSmoothedMomXYZ(const GFTrack* trk, unsigned int irep, unsigned int ihit){

  std::auto_ptr<GFAbsTrackRep> rep(trk->getTrackRep(irep)->clone());

  TVectorT<double> smoothed_state;
  TMatrixTSym<double> smoothed_cov;
  TMatrixT<double> auxInfo;
  GFDetPlane smoothing_plane;

  getSmoothedData(trk, irep, ihit, smoothed_state, smoothed_cov, smoothing_plane, auxInfo);

  if(rep->hasAuxInfo()) {
    rep->setData(smoothed_state, smoothing_plane, &smoothed_cov, &auxInfo);
  } else {
    rep->setData(smoothed_state, smoothing_plane, &smoothed_cov);
  }

  return rep->getMom(smoothing_plane);
}

TVector3 GFTools::getBiasedSmoothedMomXYZ(const GFTrack* trk, unsigned int irep, unsigned int ihit){

  std::auto_ptr<GFAbsTrackRep> rep(trk->getTrackRep(irep)->clone());

  TVectorT<double> smoothed_state;
  TMatrixTSym<double> smoothed_cov;
  TMatrixT<double> auxInfo;
  GFDetPlane smoothing_plane;

  getBiasedSmoothedData(trk, irep, ihit, smoothed_state, smoothed_cov, smoothing_plane, auxInfo);

  if(rep->hasAuxInfo()) {
    rep->setData(smoothed_state, smoothing_plane, &smoothed_cov, &auxInfo);
  } else {
    rep->setData(smoothed_state, smoothing_plane, &smoothed_cov);
  }

  return rep->getMom(smoothing_plane);
}


TVectorT<double> GFTools::getBiasedSmoothedPos(const GFTrack* trk, unsigned int irep, unsigned int ihit) {

  TVectorT<double> smoothed_state;
  TMatrixTSym<double> smoothed_cov;
  TVectorT<double> pos;

	if(GFTools::getBiasedSmoothedData(trk, irep, ihit, smoothed_state, smoothed_cov)) {

		TMatrixT<double> H = trk->getHit(ihit)->getHMatrix(trk->getTrackRep(irep));
		//H.Print();smoothed_state.Print();
		TVectorT<double> pos_tmp(H * smoothed_state);
		pos.ResizeTo(pos_tmp);
		pos = pos_tmp;

	}
	return pos;

}

TMatrixTSym<double> GFTools::getSmoothedCov(const GFTrack* trk, unsigned int irep, unsigned int ihit) {

	TVectorT<double> smoothed_state;
	TMatrixTSym<double> smoothed_cov;
	GFDetPlane pl;

	GFTools::getSmoothedData(trk, irep, ihit, smoothed_state, smoothed_cov, pl);
	TMatrixT<double> H = trk->getHit(ihit)->getHMatrix(trk->getTrackRep(irep));

	TMatrixTSym<double> cov(smoothed_cov);
	cov.Similarity(H);

	return cov;

}

TMatrixTSym<double> GFTools::getBiasedSmoothedCov(const GFTrack* trk, unsigned int irep, unsigned int ihit) {

	TVectorT<double> smoothed_state;
	TMatrixTSym<double> smoothed_cov;
	GFDetPlane pl;

	GFTools::getBiasedSmoothedData(trk, irep, ihit, smoothed_state, smoothed_cov, pl);
	TMatrixT<double> H = trk->getHit(ihit)->getHMatrix(trk->getTrackRep(irep));
	TMatrixTSym<double> cov(smoothed_cov);
	cov.Similarity(H);

	return cov;

}

bool GFTools::getSmoothedData(const GFTrack* trk, unsigned int irep, unsigned int ihit, TVectorT<double>& smoothed_state, TMatrixTSym<double>& smoothed_cov) {

	GFDetPlane pl;
	TMatrixT<double> auxInfo;
	return GFTools::getSmoothedData(trk, irep, ihit, smoothed_state, smoothed_cov, pl, auxInfo);

}

bool GFTools::getBiasedSmoothedData(const GFTrack* trk, unsigned int irep, unsigned int ihit, TVectorT<double>& smoothed_state, TMatrixTSym<double>& smoothed_cov) {

	GFDetPlane pl;
	TMatrixT<double> auxInfo;
	return GFTools::getBiasedSmoothedData(trk, irep, ihit, smoothed_state, smoothed_cov, pl, auxInfo);

}

bool GFTools::getSmoothedData(const GFTrack* trk, unsigned int irep, unsigned int ihit, TVectorT<double>& smoothed_state, TMatrixTSym<double>& smoothed_cov, GFDetPlane& smoothing_plane) {

	TMatrixT<double> auxInfo;
	return GFTools::getSmoothedData(trk, irep, ihit, smoothed_state, smoothed_cov, smoothing_plane, auxInfo);

}

bool GFTools::getBiasedSmoothedData(const GFTrack* trk, unsigned int irep, unsigned int ihit, TVectorT<double>& smoothed_state, TMatrixTSym<double>& smoothed_cov, GFDetPlane& smoothing_plane) {

	TMatrixT<double> auxInfo;
	return GFTools::getBiasedSmoothedData(trk, irep, ihit, smoothed_state, smoothed_cov, smoothing_plane, auxInfo);

}

bool GFTools::getSmoothedData(const GFTrack* trk, unsigned int irep, unsigned int ihit, TVectorT<double>& smoothed_state, TMatrixTSym<double>& smoothed_cov, GFDetPlane& smoothing_plane, TMatrixT<double>& auxInfo) {

	if(!trk->getSmoothing()) {
		std::cout << "Trying to get smoothed hit coordinates from a track without smoothing! Aborting..." << std::endl;
		TMatrixT<double> failed(1,1);
		return false;
	}

	if(ihit >= trk->getNumHits()) {
		std::cout << "Hit number out of bounds while trying to get smoothed hit coordinates! Aborting..." <<std::endl;
		TMatrixT<double> failed(1,1);
		failed(0,0) = 0;
		return false;
	}

	std::auto_ptr<GFAbsTrackRep> rep(trk->getTrackRep(irep)->clone());

	TVectorT<double> fUpSt;
	TMatrixTSym<double> fUpCov;
	TMatrixT<double> fAuxInfo;
	TMatrixT<double>* fAuxInfoP;
	TVectorT<double> bUpSt;
	TMatrixTSym<double> bUpCov;
	TMatrixT<double> bAuxInfo;
	TMatrixT<double>* bAuxInfoP;

	TVectorT<double> fSt;
	TMatrixTSym<double> fCov;
	TVectorT<double> bSt;
	TMatrixTSym<double> bCov;

  GFDetPlane fPl;
  GFDetPlane bPl;

	if(trk->getTrackRep(irep)->hasAuxInfo()) {
		trk->getBK(irep)->getMatrix("fAuxInfo",ihit,auxInfo);
	}

	if(!(trk->getSmoothingFast())) {
		if(ihit == 0) {
			trk->getBK(irep)->getVector("bUpSt",ihit+1,bUpSt);
			trk->getBK(irep)->getSymMatrix("bUpCov",ihit+1,bUpCov);
			trk->getBK(irep)->getDetPlane("fPl",ihit,smoothing_plane);
			trk->getBK(irep)->getDetPlane("bPl",ihit+1,bPl);
			if(trk->getTrackRep(irep)->hasAuxInfo()) {
				trk->getBK(irep)->getMatrix("bAuxInfo",ihit+1,bAuxInfo);
				bAuxInfoP = &bAuxInfo;
			} else {
				bAuxInfoP = NULL;
			}
			if(bUpSt.GetNrows() == 0) return false;
			rep->setData(bUpSt,bPl,&bUpCov,bAuxInfoP);
			rep->extrapolate(smoothing_plane,smoothed_state,smoothed_cov);
			return true;
		}

		if(ihit == trk->getNumHits()-1) {
			trk->getBK(irep)->getVector("fUpSt",ihit-1,fUpSt);
			trk->getBK(irep)->getSymMatrix("fUpCov",ihit-1,fUpCov);
			trk->getBK(irep)->getDetPlane("fPl",ihit-1,fPl);
			trk->getBK(irep)->getDetPlane("fPl",ihit,smoothing_plane);
			if(trk->getTrackRep(irep)->hasAuxInfo()) {
				trk->getBK(irep)->getMatrix("fAuxInfo",ihit-1,fAuxInfo);
				fAuxInfoP = &fAuxInfo;
			} else {
				fAuxInfoP = NULL;
			}
			if(fUpSt.GetNrows() == 0) return false;
			rep->setData(fUpSt,fPl,&fUpCov,fAuxInfoP);
			rep->extrapolate(smoothing_plane,smoothed_state,smoothed_cov);
			return true;
		}

		trk->getBK(irep)->getVector("fUpSt",ihit-1,fUpSt);
		trk->getBK(irep)->getSymMatrix("fUpCov",ihit-1,fUpCov);
		trk->getBK(irep)->getVector("bUpSt",ihit+1,bUpSt);
		trk->getBK(irep)->getSymMatrix("bUpCov",ihit+1,bUpCov);
		if(trk->getTrackRep(irep)->hasAuxInfo()) {
			trk->getBK(irep)->getMatrix("fAuxInfo",ihit-1,fAuxInfo);
			trk->getBK(irep)->getMatrix("bAuxInfo",ihit+1,bAuxInfo);
			fAuxInfoP = &fAuxInfo;
			bAuxInfoP = &bAuxInfo;
		} else {
			fAuxInfoP = bAuxInfoP = NULL;
		}
		trk->getBK(irep)->getDetPlane("fPl",ihit-1,fPl);
		trk->getBK(irep)->getDetPlane("fPl",ihit,smoothing_plane);
		trk->getBK(irep)->getDetPlane("bPl",ihit+1,bPl);

		if(fUpSt.GetNrows() == 0 || bUpSt.GetNrows() == 0) return false;

		rep->setData(fUpSt,fPl,&fUpCov,fAuxInfoP);
		rep->extrapolate(smoothing_plane,fSt,fCov);
		rep->setData(bUpSt,bPl,&bUpCov,bAuxInfoP);
		rep->extrapolate(smoothing_plane,bSt,bCov);

	} else {

		if(ihit == 0) {
			trk->getBK(irep)->getVector("bUpSt",ihit+1,bUpSt);
			trk->getBK(irep)->getSymMatrix("bUpCov",ihit+1,bUpCov);
			trk->getBK(irep)->getDetPlane("fPl",ihit,smoothing_plane);
			trk->getBK(irep)->getDetPlane("bPl",ihit+1,bPl);
			if(trk->getTrackRep(irep)->hasAuxInfo()) {
				trk->getBK(irep)->getMatrix("bAuxInfo",ihit+1,bAuxInfo);
				bAuxInfoP = &bAuxInfo;
			} else {
				bAuxInfoP = NULL;
			}
			if(bUpSt.GetNrows() == 0) return false;
			rep->setData(bUpSt,bPl,&bUpCov,bAuxInfoP);
			rep->extrapolate(smoothing_plane,smoothed_state,smoothed_cov);
			return true;
		}

		if(ihit == trk->getNumHits()-1) {
			trk->getBK(irep)->getDetPlane("fPl", ihit, smoothing_plane);
			trk->getBK(irep)->getVector("fSt", ihit, smoothed_state);
			trk->getBK(irep)->getSymMatrix("fCov", ihit, smoothed_cov);
			return true;
		}

		trk->getBK(irep)->getDetPlane("fPl", ihit, smoothing_plane);
		trk->getBK(irep)->getDetPlane("bPl", ihit, bPl);
		trk->getBK(irep)->getVector("fSt", ihit, fSt);
		trk->getBK(irep)->getSymMatrix("fCov", ihit, fCov);

		if(smoothing_plane == bPl) {
			trk->getBK(irep)->getVector("bSt", ihit, bSt);
			trk->getBK(irep)->getSymMatrix("bCov", ihit, bCov);
		} else {
			trk->getBK(irep)->getVector("bUpSt",ihit+1,bUpSt);
			trk->getBK(irep)->getSymMatrix("bUpCov",ihit+1,bUpCov);
			trk->getBK(irep)->getDetPlane("bPl", ihit+1, bPl);
			if(trk->getTrackRep(irep)->hasAuxInfo()) {
				trk->getBK(irep)->getMatrix("bAuxInfo", ihit+1, bAuxInfo);
				bAuxInfoP = &bAuxInfo;
			} else {
				bAuxInfoP = NULL;
			}
			rep->setData(bUpSt, bPl, &bUpCov, bAuxInfoP);
			rep->extrapolate(smoothing_plane, bSt, bCov);
		}

	}

	TMatrixTSym<double> fCovInvert;
	TMatrixTSym<double> bCovInvert;

	GFTools::invertMatrix(fCov, fCovInvert);
	GFTools::invertMatrix(bCov, bCovInvert);

	GFTools::invertMatrix(fCovInvert + bCovInvert, smoothed_cov);

	smoothed_state.ResizeTo(smoothed_cov * (fCovInvert*fSt + bCovInvert*bSt));
	smoothed_state = smoothed_cov * (fCovInvert*fSt + bCovInvert*bSt);

	return true;

}

bool GFTools::getBiasedSmoothedData(const GFTrack* trk, unsigned int irep, unsigned int ihit, TVectorT<double>& smoothed_state, TMatrixTSym<double>& smoothed_cov, GFDetPlane& smoothing_plane, TMatrixT<double>& auxInfo) {

	if(!trk->getSmoothing()) {
		std::cout << "Trying to get smoothed hit coordinates from a track without smoothing! Aborting..." << std::endl;
		TMatrixT<double> failed(1,1);
		return false;
	}

	if(ihit >= trk->getNumHits()) {
		std::cout << "Hit number out of bounds while trying to get smoothed hit coordinates! Aborting..." <<std::endl;
		TMatrixT<double> failed(1,1);
		failed(0,0) = 0;
		return false;
	}

	TVectorT<double> bUpSt;
	TMatrixTSym<double> bUpCov;
	TMatrixT<double> bAuxInfo;
	TMatrixT<double>* bAuxInfoP;
	GFDetPlane bPl;

	std::auto_ptr<GFAbsTrackRep> rep(trk->getTrackRep(irep)->clone());

	if(trk->getTrackRep(irep)->hasAuxInfo()) {
		trk->getBK(irep)->getMatrix("fAuxInfo",ihit,auxInfo);
	}

	if(ihit == 0) {
		trk->getBK(irep)->getVector("bUpSt",ihit,smoothed_state);
		trk->getBK(irep)->getSymMatrix("bUpCov",ihit,smoothed_cov);
		trk->getBK(irep)->getDetPlane("bPl",ihit,smoothing_plane);
		return true;
	}

	if(ihit == trk->getNumHits()-1) {
		trk->getBK(irep)->getVector("fUpSt",ihit,smoothed_state);
		trk->getBK(irep)->getSymMatrix("fUpCov",ihit,smoothed_cov);
		trk->getBK(irep)->getDetPlane("fPl",ihit,smoothing_plane);
		return true;
	}

	TVectorT<double> fSt;
	TMatrixTSym<double> fCov;
	TVectorT<double> bSt;
	TMatrixTSym<double> bCov;

	if(!(trk->getSmoothingFast())) {

		trk->getBK(irep)->getVector("fUpSt",ihit,fSt);
		trk->getBK(irep)->getSymMatrix("fUpCov",ihit,fCov);
		trk->getBK(irep)->getVector("bUpSt",ihit+1,bUpSt);
		trk->getBK(irep)->getSymMatrix("bUpCov",ihit+1,bUpCov);
		if(trk->getTrackRep(irep)->hasAuxInfo()) {
			trk->getBK(irep)->getMatrix("bAuxInfo",ihit+1,bAuxInfo);
			bAuxInfoP = &bAuxInfo;
		} else {
			bAuxInfoP = NULL;
		}
		trk->getBK(irep)->getDetPlane("fPl",ihit,smoothing_plane);
		trk->getBK(irep)->getDetPlane("bPl",ihit+1,bPl);

		if(bUpSt.GetNrows() == 0) {
			return false;
		}

		rep->setData(bUpSt,bPl,&bUpCov,bAuxInfoP);
		rep->extrapolate(smoothing_plane,bSt,bCov);

	} else {

		trk->getBK(irep)->getVector("fUpSt",ihit,fSt);
		trk->getBK(irep)->getSymMatrix("fUpCov",ihit,fCov);
		trk->getBK(irep)->getDetPlane("fPl",ihit,smoothing_plane);
		trk->getBK(irep)->getDetPlane("bPl",ihit,bPl);

		if(smoothing_plane == bPl) {
			trk->getBK(irep)->getVector("bSt",ihit,bSt);
			trk->getBK(irep)->getSymMatrix("bCov",ihit,bCov);
		} else {
			trk->getBK(irep)->getVector("bUpSt",ihit+1,bUpSt);
			trk->getBK(irep)->getSymMatrix("bUpCov",ihit+1,bUpCov);
			if(trk->getTrackRep(irep)->hasAuxInfo()) {
				trk->getBK(irep)->getMatrix("bAuxInfo",ihit+1,bAuxInfo);
				bAuxInfoP = &bAuxInfo;
			} else {
				bAuxInfoP = NULL;
			}
			trk->getBK(irep)->getDetPlane("fPl",ihit,smoothing_plane);
			trk->getBK(irep)->getDetPlane("bPl",ihit+1,bPl);

			if(bUpSt.GetNrows() == 0) {
				return false;
			}

			rep->setData(bUpSt,bPl,&bUpCov,bAuxInfoP);
			rep->extrapolate(smoothing_plane,bSt,bCov);
		}

	}

	TMatrixTSym<double> fCovInvert;
	TMatrixTSym<double> bCovInvert;

	GFTools::invertMatrix(fCov, fCovInvert);
	GFTools::invertMatrix(bCov, bCovInvert);

	GFTools::invertMatrix(fCovInvert + bCovInvert, smoothed_cov);

	smoothed_state.ResizeTo(smoothed_cov * (fCovInvert*fSt + bCovInvert*bSt));
	smoothed_state = smoothed_cov * (fCovInvert*fSt + bCovInvert*bSt);

	return true;

}

GFDetPlane GFTools::getSmoothingPlane(const GFTrack* trk, unsigned int irep, unsigned int ihit) {

	GFDetPlane pl;
	trk->getBK(irep)->getDetPlane("fPl",ihit,pl);
	return pl;

}

double GFTools::getTrackLength(const GFTrack* trk, unsigned int irep, unsigned int startHit, unsigned int endHit){

  if(!trk->getSmoothing()) {
    GFException exc("Trying to get tracklength from a track without smoothing!",__LINE__,__FILE__);
    throw exc;
  }

  if(startHit >= trk->getNumHits() || endHit >= trk->getNumHits()) {
    GFException exc("Hit number out of bounds while trying to get tracklength!",__LINE__,__FILE__);
    throw exc;
  }

  bool inv(false);
  if(startHit > endHit) {
    unsigned int biggerOne = startHit;
    startHit = endHit;
    endHit = biggerOne;
    inv = true;
  }

  if (startHit==0 && endHit==0) endHit = trk->getNumHits()-1;
  if (startHit == endHit) return 0.;

  double totLen(0), fLen(0), bLen(0);

  for (unsigned int i=startHit; i!=endHit; ++i){
    trk->getBK(irep)->getNumber("fExtLen", i+1, fLen);
    trk->getBK(irep)->getNumber("bExtLen", i, bLen);
    totLen += 0.5*(fLen - bLen);
  }

  if (inv) return -totLen;
  return totLen;

}

void GFTools::invertMatrix(const TMatrixTSym<double>& mat, TMatrixTSym<double>& inv){
	inv.ResizeTo(mat);

	// check if numerical limits are reached (i.e at least one entry < 1E-100 and/or at least one entry > 1E100)
	if (!(mat<1.E100) || !(mat>-1.E100)){
		GFException e("cannot invert matrix GFTools::invertMatrix(), entries too big (>1e100)",
				__LINE__,__FILE__);
		e.setFatal();
		throw e;	
	}
	// do the trivial inversions for 1x1 and 2x2 matrices manually
	if (mat.GetNrows() == 1){
	  inv(0,0) = 1/mat(0,0);
	  return;
	} if (mat.GetNrows() == 2){
	  double det = mat(0,0)*mat(1,1) - mat(1,0)*mat(1,0);
	  if(fabs(det) < 1E-50){
	    GFException e("cannot invert matrix GFTools::invertMatrix(), determinant = 0",
	        __LINE__,__FILE__);
	    e.setFatal();
	    throw e;
	  }
	  det = 1./det;
	  inv(0,0) =             det * mat(1,1);
	  inv(0,1) = inv(1,0) = -det * mat(1,0);
	  inv(1,1) =             det * mat(0,0);
	  return;
	}

	// else use TDecompChol
	bool status = 0;
	TDecompChol invertAlgo(mat);

	invertAlgo.SetTol(1E-50); //this is a hack because a tolerance of 1E-22 does not make any sense for doubles only for floats

	status = invertAlgo.Invert(inv);
	if(status == 0){
		GFException e("cannot invert matrix GFTools::invertMatrix(), status = 0",
				__LINE__,__FILE__);
		e.setFatal();
		throw e;
	}
}

void GFTools::invertMatrix(TMatrixTSym<double>& mat){
	// check if numerical limits are reached (i.e at least one entry < 1E-100 and/or at least one entry > 1E100)
	if (!(mat<1.E100) || !(mat>-1.E100)){
		GFException e("cannot invert matrix GFTools::invertMatrix(), entries too big (>1e100)",
				__LINE__,__FILE__);
		e.setFatal();
		throw e;	
	}
	// do the trivial inversions for 1x1 and 2x2 matrices manually
	if (mat.GetNrows() == 1){
	  mat(0,0) = 1/mat(0,0);
	  return;
	} else if (mat.GetNrows() == 2){
	  double *arr = mat.GetMatrixArray();
	  double det = arr[0]*arr[3] - arr[1]*arr[1];
	  if(fabs(det) < 1E-50){
	    GFException e("cannot invert matrix GFTools::invertMatrix(), determinant = 0",
	        __LINE__,__FILE__);
	    e.setFatal();
	    throw e;
	  }
	  det = 1./det;
	  double temp[3];
	  temp[0] =  det * arr[3];
	  temp[1] = -det * arr[1];
	  temp[2] =  det * arr[0];
	  //double *arr = mat.GetMatrixArray();
	  arr[0] = temp[0];
	  arr[1] = arr[2] = temp[1];
	  arr[3] = temp[2];
	  return;
	}

	// else use TDecompChol
	bool status = 0;
	TDecompChol invertAlgo(mat);

	invertAlgo.SetTol(1E-50); //this is a hack because a tolerance of 1E-22 does not make any sense for doubles only for floats

	status = invertAlgo.Invert(mat);
	if(status == 0){
		GFException e("cannot invert matrix GFTools::invertMatrix(), status = 0",
				__LINE__,__FILE__);
		e.setFatal();
		throw e;
	}
}

void GFTools::updateRepSmoothed(GFTrack* trk, unsigned int irep, unsigned int ihit) {

  TVectorT<double> smoothed_state;
  TMatrixT<double> auxInfo;
	TMatrixTSym<double> smoothed_cov;
	GFDetPlane smoothing_plane;

	getBiasedSmoothedData(trk, irep, ihit, smoothed_state, smoothed_cov, smoothing_plane, auxInfo);

	GFAbsTrackRep* rep = trk->getTrackRep(irep);
	if(rep->hasAuxInfo()) {
		rep->setData(smoothed_state, smoothing_plane, &smoothed_cov, &auxInfo);
	} else {
		rep->setData(smoothed_state, smoothing_plane, &smoothed_cov);
	}

	trk->setRepAtHit(irep, ihit);
}

double GFTools::getSmoothedChiSqu(const GFTrack* trk, unsigned int irep, unsigned int ihit){
	TVectorT<double> smoothed_state;
	TMatrixTSym<double> smoothed_cov;
	GFTools::getBiasedSmoothedData(trk, irep, ihit, smoothed_state, smoothed_cov);
	TMatrixT<double> H = trk->getHit(ihit)->getHMatrix(trk->getTrackRep(irep));
	TVectorT<double> pos = H * smoothed_state;
	TVectorT<double> m = trk->getHit(ihit)->getRawHitCoord(); //measurement of hit
	TMatrixTSym<double> V = trk->getHit(ihit)->getRawHitCov(); //covariance matrix of hit
	TVectorT<double> res = m - pos;
	TMatrixTSym<double> R = V - smoothed_cov.Similarity(H);
	TMatrixTSym<double> invR;
	invertMatrix(R,invR);
	double smoothedChiSqu = invR.Similarity(res);
	return smoothedChiSqu;
}

unsigned int GFTools::getClosestHit(const GFTrack* trk, unsigned int irep, const TVector3& pos, double& distance, bool checkEveryHit){

  if(!trk->getSmoothing()) {
    GFException exc("Trying to get closest hit from a track without smoothing!",__LINE__,__FILE__);
    throw exc;
  }

	unsigned int nHits = trk->getNumHits();
	if (nHits == 1) return 0;

	TVector3 hitPos;
	double minDist(1.E99), dist;
	unsigned int minId(0);

	if (checkEveryHit || nHits<4){
		// brute force search
		for (unsigned int i=0; i<nHits; ++i){
			hitPos = getSmoothedPosXYZ(trk, irep, i);
			dist = (pos-hitPos).Mag();
			if (dist<minDist) {
				minDist = dist;
				minId = i;
			}
		}
	}
	else { // hill climbing algorithm
		double distFirst = (pos-getSmoothedPosXYZ(trk, irep, 0)).Mag();
		double distLast = (pos-getSmoothedPosXYZ(trk, irep, nHits-1)).Mag();
		if (distFirst <= distLast){
			minDist = distFirst;
			minId = 0;
			for (unsigned int i=1; i<nHits-1; ++i){
				hitPos = getSmoothedPosXYZ(trk, irep, i);
				dist = (pos-hitPos).Mag();
				if (dist<=minDist) {
					minDist = dist;
					minId = i;
				}
				else break; // distance is increasing again!
			}
		}
		else {
			minDist = distLast;
			minId = nHits-1;
			for (unsigned int i=nHits-2; i>1; --i){
				hitPos = getSmoothedPosXYZ(trk, irep, i);
				dist = (pos-hitPos).Mag();
				if (dist<=minDist) {
					minDist = dist;
					minId = i;
				}
				else break; // distance is increasing again!
			}
		}
	}

	distance = minDist;
	return minId;
}

