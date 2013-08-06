//
// Write fit results to tree, read again, compare.
//

#include <iostream>
#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>

#include <AbsFinitePlane.h>
#include <AbsKalmanFitter.h>
#include <AbsFitterInfo.h>
#include <AbsMeasurement.h>
#include <AbsTrackRep.h>
#include <ConstField.h>
#include <DetPlane.h>
#include <Exception.h>
#include <FieldManager.h>
#include <KalmanFittedStateOnPlane.h>
#include <KalmanFitterInfo.h>
#include <KalmanFitter.h>
#include <KalmanFitterRefTrack.h>
#include <KalmanFitStatus.h>
#include <DAF.h>
#include <MaterialInfo.h>
#include <MeasuredStateOnPlane.h>
#include <MeasurementOnPlane.h>
#include <PlanarMeasurement.h>
#include <ProlateSpacepointMeasurement.h>
#include <RectangularFinitePlane.h>
#include <ReferenceStateOnPlane.h>
#include <SharedPlanePtr.h>
#include <SpacepointMeasurement.h>
#include <StateOnPlane.h>
#include <Tools.h>
#include <TrackCand.h>
#include <TrackCandHit.h>
#include <Track.h>
#include <TrackPoint.h>
#include <WireMeasurement.h>
#include <WirePointMeasurement.h>

#include <MaterialEffects.h>
#include <RKTools.h>
#include <RKTrackRep.h>
#include <StepLimits.h>
#include <TGeoMaterialInterface.h>

#include <TApplication.h>
#include <TCanvas.h>
#include <TDatabasePDG.h>
#include <TEveManager.h>
#include <TGeoManager.h>
#include <TH1D.h>
#include <TRandom.h>
#include <TStyle.h>
#include <TVector3.h>
#include <vector>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include "TDatabasePDG.h"
#include <TMath.h>
#include <TString.h>

#define FILENAME "/tmp/streamerTest.root"

void handler(int sig) {
  void *array[10];
  size_t size;

  // get void*'s for all entries on the stack
  size = backtrace(array, 10);

  // print out all the frames to stderr
  fprintf(stderr, "Error: signal %d:\n", sig);
  backtrace_symbols_fd(array, size, 2);
  exit(1);
}

int randomPdg() {
  int pdg;

  switch(int(gRandom->Uniform(8))) {
  case 1:
    pdg = -11; break;
  case 2:
    pdg = 11; break;
  case 3:
    pdg = 13; break;
  case 4:
    pdg = -13; break;
  case 5:
    pdg = 211; break;
  case 6:
    pdg = -211; break;
  case 7:
    pdg = 2212; break;
  default:
    pdg = 211;
  }

  return pdg;
}


int randomSign() {
  if (gRandom->Uniform(1) > 0.5)
    return 1;
  return -1;
}


bool emptyTrackTest()
{
  TFile *f = TFile::Open(FILENAME, "RECREATE");
  f->cd();
  genfit::Track *t = new genfit::Track();
  if (!t->checkConsistency())
    return false;
  t->Write("direct");
  f->Close();
  delete f;

  f = TFile::Open(FILENAME, "READ");
  t = (genfit::Track*)f->Get("direct");
  bool result = t->checkConsistency();
  delete t;
  delete f;
  return result;
}


int main() {
  if (!emptyTrackTest())
    {
      std::cout << "enptyTrackTest failed." << std::endl;
      return 1;
    }

  enum eFitterType { SimpleKalman = 0,
		     RefKalman,
		     DafSimple,
		     DafRef};

  const unsigned int nEvents = 1;  // events to accumulate in file.
  const double BField = 15.;       // kGauss
  const double momentum = 0.4;     // GeV
  const double theta = 100;         // degree
  const double thetaDetPlane = 90;         // degree
  const double phiDetPlane = 0;         // degree
  const double pointDist = 5;      // cm; approx. distance between measurements generated w/ RKTrackRep
  const double pointDistDeg = 5;      // degree; distance between measurements generated w/ helix model
  const double resolution = 0.2;   // cm; resolution of generated measurements

  const double resolutionWire = 5*resolution;   // cm; resolution of generated measurements
  const TVector3 wireDir(0,0,1);
  const double skewAngle(5);
  const bool useSkew = true;
  const int nSuperLayer = 5;
  const double minDrift = 0.;
  const double maxDrift = 2;
  const bool idealLRResolution = false; // resolve the l/r ambiguities of the wire measurements

  const double outlierProb = -0.;
  const double outlierRange = 5;

  const double hitSwitchProb = -0.; // probability to give hits to fit in wrong order (flip two hits)

  //const eFitterType fitterId = SimpleKalman;
  const eFitterType fitterId = RefKalman;
  //const eFitterType fitterId = DafRef;
  const genfit::eMultipleMeasurementHandling mmHandling = genfit::weightedAverage;
  //const genfit::eMultipleMeasurementHandling mmHandling = genfit::unweightedClosestToReference;
  //const genfit::eMultipleMeasurementHandling mmHandling = genfit::unweightedClosestToPrediction;
  const int nIter = 10; // max number of iterations
  const double dPVal = 1.E-3; // convergence criterion

  const bool resort = false;

  const int pdg = 13;               // particle pdg code

  const bool smearPosMom = true;     // init the Reps with smeared pos and mom
  const double posSmear = 10*resolution;     // cm
  const double momSmear = 5. /180.*TMath::Pi();     // rad
  const double momMagSmear = 0.2;   // relative
  const double zSmearFac = 2;

  const bool HelixTest = true;      // use helix for creating measurements

  const bool matFX = false;         // include material effects; can only be disabled for RKTrackRep!

  const bool debug = false;

  enum eMeasurementType { Pixel = 0,
			  Spacepoint,
			  ProlateSpacepoint,
			  StripU,
			  StripV,
			  Wire,
			  WirePoint };


  std::vector<unsigned int> measurementTypes;

  /*measurementTypes.push_back(Pixel);
    measurementTypes.push_back(Pixel);
    measurementTypes.push_back(Spacepoint);
    measurementTypes.push_back(Spacepoint);
    measurementTypes.push_back(ProlateSpacepoint);
    measurementTypes.push_back(ProlateSpacepoint);
    measurementTypes.push_back(StripU);
    measurementTypes.push_back(StripU);
    measurementTypes.push_back(StripV);
    measurementTypes.push_back(StripV);
    measurementTypes.push_back(Wire);
    measurementTypes.push_back(Wire);
    measurementTypes.push_back(WirePoint);
    measurementTypes.push_back(WirePoint);*/
  for (int i = 0; i < 12; ++i)
    measurementTypes.push_back(Pixel);


  // init fitter
  genfit::AbsKalmanFitter* fitter;
  switch (fitterId) {
  case SimpleKalman:
    fitter = new genfit::KalmanFitter(nIter, dPVal);
    fitter->setMultipleMeasurementHandling(mmHandling);
    break;

  case RefKalman:
    fitter = new genfit::KalmanFitterRefTrack(nIter, dPVal);
    fitter->setMultipleMeasurementHandling(mmHandling);
    break;

  case DafSimple:
    {
      genfit::AbsKalmanFitter* DafsKalman = new genfit::KalmanFitter();
      DafsKalman->setMultipleMeasurementHandling(mmHandling);
      fitter = new genfit::DAF(DafsKalman);
    }
    break;
  case DafRef:
    fitter = new genfit::DAF();
    break;

  }

  if (dynamic_cast<genfit::DAF*>(fitter) != NULL) {
    //static_cast<genfit::DAF*>(fitter)->setBetas(100, 50, 25, 12, 6, 3, 1, 0.5, 0.1);
    //static_cast<genfit::DAF*>(fitter)->setBetas(81, 8, 4, 0.5, 0.1);
    static_cast<genfit::DAF*>(fitter)->setAnnealingScheme(100, 0.1, 5);
    //static_cast<genfit::DAF*>(fitter)->setConvergenceDeltaWeight(0.0001);
    fitter->setMaxIterations(nIter);
  }


  gRandom->SetSeed(14);
  signal(SIGSEGV, handler);   // install our handler

  // init geometry and mag. field
  /* TGeoManager* geom = */ new TGeoManager("Geometry", "Geane geometry");
  TGeoManager::Import("genfitGeom.root");
  genfit::FieldManager::getInstance()->init(new genfit::ConstField(0.,0.,BField));
  genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());

  const double charge = TDatabasePDG::Instance()->GetParticle(pdg)->Charge()/(3.);


  // prepare output tree for Tracks 
  // std::unique_ptr<genfit::Track> fitTrack(new genfit::Track());
  genfit::Track* fitTrack = new genfit::Track();
  TVectorD stateFinal;
  TMatrixDSym covFinal;
  genfit::DetPlane planeFinal;

  TFile* fOut = new TFile(FILENAME, "RECREATE");
  fOut->cd();
  TTree* tResults = new TTree("tResults", "results from track fit");
  tResults->Branch("gfTrack", "genfit::Track", &fitTrack, 32000, -1);
  tResults->Branch("stateFinal", &stateFinal);
  tResults->Branch("covFinal", &covFinal, 32000, -1);
  tResults->Branch("planeFinal", &planeFinal, 32000, -1);

  unsigned int nTotalIter(0);

  // main loop, accumulates events in file
  for (unsigned int iEvent=0; iEvent<nEvents; ++iEvent){

    if (debug || (iEvent+1)%10==0) std::cout << "Event Nr. " << iEvent+1 << std::endl;


    // true start values
    TVector3 pos(0, 0, 0);
    TVector3 mom(1.,0,0);
    mom.SetPhi(gRandom->Uniform(0.,2*TMath::Pi()));
    //mom.SetTheta(gRandom->Uniform(0.4*TMath::Pi(),0.6*TMath::Pi()));
    mom.SetTheta(theta*TMath::Pi()/180);
    mom.SetMag(momentum);

    // calc helix parameters
    TVector3 dir2D(mom);
    dir2D.SetZ(0);
    dir2D.SetMag(1.);
    double R = 100.*mom.Perp()/(0.0299792458*BField);
    double sgn = 1;
    if (charge<0) sgn=-1.;
    TVector3 center = pos + sgn * R * dir2D.Orthogonal();
    double alpha0 = (pos-center).Phi();



    // smeared start values
    TVector3 posM(pos);
    TVector3 momM(mom);
    if (smearPosMom) {
      posM.SetX(gRandom->Gaus(posM.X(),posSmear));
      posM.SetY(gRandom->Gaus(posM.Y(),posSmear));
      posM.SetZ(gRandom->Gaus(posM.Z(),zSmearFac*posSmear));

      momM.SetPhi(gRandom->Gaus(mom.Phi(),momSmear));
      momM.SetTheta(gRandom->Gaus(mom.Theta(),momSmear));
      momM.SetMag(gRandom->Gaus(mom.Mag(), momMagSmear*mom.Mag()));
    }


    // trackrep for creating measurements
    genfit::AbsTrackRep* rep = new genfit::RKTrackRep(pdg);
    genfit::StateOnPlane stateRef(rep);
    rep->setPosMom(&stateRef, pos, mom);

    // smeared start state
    genfit::StateOnPlane stateSmeared(rep);
    rep->setPosMom(&stateSmeared, posM, momM);

    //rep->setPropDir(1);

    if (!matFX) genfit::MaterialEffects::getInstance()->setNoEffects();

    // remember original initial state
    const genfit::StateOnPlane stateRefOrig(stateRef);

    // create smeared measurements
    std::vector<genfit::AbsMeasurement*> measurements;

    // true values for left right. 0 for non wire measurements
    std::vector<int> leftRightTrue;
    std::vector<bool> outlierTrue;

    TVector3 point, dir;
    int wireCounter = 0;
    if (debug) std::cout << "Start creating measurements ... \n";
    try{
      for (unsigned int i=0; i<measurementTypes.size(); ++i){
	// get current position and momentum
	if (!HelixTest) {
	  rep->getPosMom(&stateRef, point, dir);
	  dir.SetMag(1);
	}
	else {
	  double angle = alpha0 - sgn * pointDistDeg/180.*TMath::Pi()*i;
	  TVector3 radius(R,0,0);
	  radius.SetPhi(angle);
	  point = center + radius;
	  point.SetZ(pos.Z() + ((alpha0-angle)*R * TMath::Tan(mom.Theta()-TMath::Pi()*0.5) ));
	  if (debug) {std::cout << "(true) track position"; point.Print();}

	  dir = mom;
	  dir.SetPhi(mom.Phi()+(angle-alpha0));
	  dir.SetMag(1);
	}


	// create measurement
	genfit::AbsMeasurement* measurement;

	TVector3 currentWireDir(wireDir);

	if (useSkew && (int)((double)wireCounter/(double)nSuperLayer)%2 == 1) {
	  TVector3 perp(wireDir.Cross(dir));
	  if ((int)((double)wireCounter/(double)nSuperLayer)%4 == 1){
	    currentWireDir.Rotate(skewAngle*TMath::Pi()/180, wireDir.Cross(perp));
	  }
	  else currentWireDir.Rotate(-skewAngle*TMath::Pi()/180, wireDir.Cross(perp));
	}
	currentWireDir.SetMag(1.);

	TVector3 planeNorm(dir);
	planeNorm.SetTheta(thetaDetPlane*TMath::Pi()/180);
	planeNorm.SetPhi(planeNorm.Phi()+phiDetPlane);
	static const TVector3 z(0,0,1);
	static const TVector3 x(1,0,0);

	int lr = 1;
	TVector3 wirePerp;
	if (measurementTypes[i] == Wire ||
	    measurementTypes[i] == WirePoint){
	  wirePerp = dir.Cross(currentWireDir);
	  if (gRandom->Uniform(-1,1) >= 0) {
	    wirePerp *= -1.;
	    lr = -1;
	  }
	  wirePerp.SetMag(gRandom->Uniform(minDrift, maxDrift));

	  leftRightTrue.push_back(lr);
	}
	else leftRightTrue.push_back(0);

	bool outlier(false);
	if (outlierProb > gRandom->Uniform(1.)) {
	  outlier = true;
	  if(debug)  std::cerr << "create outlier" << std::endl;
	}

	outlierTrue.push_back(outlier);

	switch(measurementTypes[i]){
	case Pixel: {
	  if (debug) std::cerr << "create PixHit" << std::endl;

	  genfit::SharedPlanePtr plane(new genfit::DetPlane(point, planeNorm.Cross(z), (planeNorm.Cross(z)).Cross(planeNorm)));

	  TVectorD hitCoords(2);
	  if (outlier) {
	    hitCoords(0) = gRandom->Uniform(-outlierRange, outlierRange);
	    hitCoords(1) = gRandom->Uniform(-outlierRange, outlierRange);
	  }
	  else {
	    hitCoords(0) = gRandom->Gaus(0,resolution);
	    hitCoords(1) = gRandom->Gaus(0,resolution);
	  }

	  TMatrixDSym hitCov(2);
	  hitCov(0,0) = resolution*resolution;
	  hitCov(1,1) = resolution*resolution;

	  measurement = new genfit::PlanarMeasurement(hitCoords, hitCov, 0, i, nullptr);
	  static_cast<genfit::PlanarMeasurement*>(measurement)->setPlane(plane, i);
	}
	  break;

	case Spacepoint: {
	  if (debug) std::cerr << "create SpacepointHit" << std::endl;

	  TVectorD hitCoords(3);
	  if (outlier) {
	    hitCoords(0) = gRandom->Uniform(point.X()-outlierRange, point.X()+outlierRange);
	    hitCoords(1) = gRandom->Uniform(point.Y()-outlierRange, point.Y()+outlierRange);
	    hitCoords(2) = gRandom->Uniform(point.Z()-outlierRange, point.Z()+outlierRange);
	  }
	  else {
	    hitCoords(0) = gRandom->Gaus(point.X(),resolution);
	    hitCoords(1) = gRandom->Gaus(point.Y(),resolution);
	    hitCoords(2) = gRandom->Gaus(point.Z(),resolution);
	  }

	  TMatrixDSym hitCov(3);
	  hitCov(0,0) = resolution*resolution;
	  hitCov(1,1) = resolution*resolution;
	  hitCov(2,2) = resolution*resolution;

	  measurement = new genfit::SpacepointMeasurement(hitCoords, hitCov, 1, i, nullptr);
	}
	  break;

	case ProlateSpacepoint: {
	  if (debug) std::cerr << "create ProlateSpacepointHit" << std::endl;

	  TVectorD hitCoords(3);
	  if (outlier) {
	    hitCoords(0) = gRandom->Uniform(point.X()-outlierRange, point.X()+outlierRange);
	    hitCoords(1) = gRandom->Uniform(point.Y()-outlierRange, point.Y()+outlierRange);
	    hitCoords(2) = gRandom->Uniform(point.Z()-outlierRange, point.Z()+outlierRange);
	  }
	  else {
	    hitCoords(0) = point.X();
	    hitCoords(1) = point.Y();
	    hitCoords(2) = point.Z();
	  }

	  TMatrixDSym hitCov(3);
	  hitCov(0,0) = resolution*resolution;
	  hitCov(1,1) = resolution*resolution;
	  hitCov(2,2) = resolutionWire*resolutionWire;

	  // rotation matrix
	  TVector3 xp = currentWireDir.Orthogonal();
	  xp.SetMag(1);
	  TVector3 yp = currentWireDir.Cross(xp);
	  yp.SetMag(1);

	  TMatrixD rot(3,3);

	  rot(0,0) = xp.X();  rot(0,1) = yp.X();  rot(0,2) = currentWireDir.X();
	  rot(1,0) = xp.Y();  rot(1,1) = yp.Y();  rot(1,2) = currentWireDir.Y();
	  rot(2,0) = xp.Z();  rot(2,1) = yp.Z();  rot(2,2) = currentWireDir.Z();

	  // smear
	  TVectorD smearVec(3);
	  smearVec(0) = resolution;
	  smearVec(1) = resolution;
	  smearVec(2) = resolutionWire;
	  smearVec *= rot;
	  if (!outlier) {
	    hitCoords(0) += gRandom->Gaus(0, smearVec(0));
	    hitCoords(1) += gRandom->Gaus(0, smearVec(1));
	    hitCoords(2) += gRandom->Gaus(0, smearVec(2));
	  }


	  // rotate cov
	  hitCov.Similarity(rot);

	  measurement = new genfit::ProlateSpacepointMeasurement(hitCoords, hitCov, 2, i, nullptr);

	  static_cast<genfit::ProlateSpacepointMeasurement*>(measurement)->setLargestErrorDirection(currentWireDir);
	}
	  break;

	case StripU: case StripV: {
	  if (debug) std::cerr << "create StripHit" << std::endl;

	  TVector3 vU, vV;
	  if (measurementTypes[i] == StripU) {
	    vU = planeNorm.Cross(z);
	    vV = (planeNorm.Cross(z)).Cross(planeNorm);
	  } else {
	    vU = (planeNorm.Cross(z)).Cross(planeNorm);
	    vV = planeNorm.Cross(z);
	  }
	  genfit::SharedPlanePtr plane(new genfit::DetPlane(point, vU, vV));

	  TVectorD hitCoords(1);
	  if (outlier)
	    hitCoords(0) = gRandom->Uniform(-outlierRange, outlierRange);
	  else
	    hitCoords(0) = gRandom->Gaus(0,resolution);

	  TMatrixDSym hitCov(1);
	  hitCov(0,0) = resolution*resolution;

	  measurement = new genfit::PlanarMeasurement(hitCoords, hitCov, 3, i, nullptr);
	  static_cast<genfit::PlanarMeasurement*>(measurement)->setPlane(plane, i);
	}
	  break;

	case Wire: {
	  if (debug) std::cerr << "create WireHit" << std::endl;

	  if (outlier) {
	    wirePerp.SetMag(gRandom->Uniform(outlierRange));
	  }

	  TVectorD hitCoords(7);
	  hitCoords(0) = (point-wirePerp-currentWireDir).X();
	  hitCoords(1) = (point-wirePerp-currentWireDir).Y();
	  hitCoords(2) = (point-wirePerp-currentWireDir).Z();

	  hitCoords(3) = (point-wirePerp+currentWireDir).X();
	  hitCoords(4) = (point-wirePerp+currentWireDir).Y();
	  hitCoords(5) = (point-wirePerp+currentWireDir).Z();

	  if (outlier)
	    hitCoords(6) = gRandom->Uniform(outlierRange);
	  else
	    hitCoords(6) = gRandom->Gaus(wirePerp.Mag(), resolution);

	  TMatrixDSym hitCov(7);
	  hitCov(6,6) = resolution*resolution;


	  measurement = new genfit::WireMeasurement(hitCoords, hitCov, 4, i, nullptr);
	  if (idealLRResolution){
	    static_cast<genfit::WireMeasurement*>(measurement)->setLeftRightResolution(lr);
	  }
	  ++wireCounter;
	}
	  break;

	case WirePoint: {
	  if (debug) std::cerr << "create WirePointHit" << std::endl;

	  if (outlier) {
	    wirePerp.SetMag(gRandom->Uniform(outlierRange));
	  }

	  TVectorD hitCoords(8);
	  hitCoords(0) = (point-wirePerp-currentWireDir).X();
	  hitCoords(1) = (point-wirePerp-currentWireDir).Y();
	  hitCoords(2) = (point-wirePerp-currentWireDir).Z();

	  hitCoords(3) = (point-wirePerp+currentWireDir).X();
	  hitCoords(4) = (point-wirePerp+currentWireDir).Y();
	  hitCoords(5) = (point-wirePerp+currentWireDir).Z();

	  if (outlier) {
	    hitCoords(6) = gRandom->Uniform(outlierRange);
	    hitCoords(7) = gRandom->Uniform(currentWireDir.Mag()-outlierRange, currentWireDir.Mag()+outlierRange);
	  }
	  else {
	    hitCoords(6) = gRandom->Gaus(wirePerp.Mag(), resolution);
	    hitCoords(7) = gRandom->Gaus(currentWireDir.Mag(), resolutionWire);
	  }


	  TMatrixDSym hitCov(8);
	  hitCov(6,6) = resolution*resolution;
	  hitCov(7,7) = resolutionWire*resolutionWire;

	  measurement = new genfit::WirePointMeasurement(hitCoords, hitCov, 5, i, nullptr);
	  if (idealLRResolution){
	    static_cast<genfit::WirePointMeasurement*>(measurement)->setLeftRightResolution(lr);
	  }
	  ++wireCounter;
	}
	  break;

	default:
	  std::cerr << "measurement type not defined!" << std::endl;
	  exit(0);
	}
	measurements.push_back(measurement);

	/*if (debug) {
	  std::cout << "(smeared) measurement coordinates"; measurement->getRawHitCoords().Print();
	  }*/

	if (!HelixTest) {
	  // stepalong (approximately)
	  dir.SetMag(pointDist);
	  genfit::SharedPlanePtr pl(new genfit::DetPlane(point+dir, dir));
	  rep->extrapolateToPlane(&stateRef, pl);
	}
      }

      assert(measurementTypes.size() == leftRightTrue.size());
    }
    catch(genfit::Exception& e){
      std::cerr<<"Exception, next track"<<std::endl;
      std::cerr << e.what();
      continue; // here is a memleak!
    }

    if (debug) std::cout << "... done creating measurements \n";



    // create track
    fitTrack = new genfit::Track(rep, rep->get6DState(&stateSmeared)); //initialized with smeared rep
    //if (debug) fitTrack->Print("C");

    //fitTrack->addTrackRep(rep->clone()); // check if everything works fine with more than one rep

    // add measurements
    for(unsigned int i=0; i<measurements.size(); ++i){
      std::vector<genfit::AbsMeasurement*> measVec;
      measVec.push_back(measurements[i]);
      if (i>0 && hitSwitchProb > gRandom->Uniform(1.))
	fitTrack->insertPoint(new genfit::TrackPoint(measVec, fitTrack), -2);
      else
	fitTrack->insertPoint(new genfit::TrackPoint(measVec, fitTrack));

      //if (debug) fitTrack->Print("C");
    }

    // print trackCand
    //if (debug) fitTrack->Print();

    assert(fitTrack->checkConsistency());

    // do the fit
    try{
      if (debug) std::cout<<"Starting the fitter"<<std::endl;
	fitter->processTrack(fitTrack, resort);
	if (debug) std::cout<<"fitter is finished!"<<std::endl;
    }
    catch(genfit::Exception& e){
      std::cerr << e.what();
      std::cerr << "Exception, next track" << std::endl;
      continue;
    }


    //if (debug) fitTrack->Print("C");
    assert(fitTrack->checkConsistency());

    nTotalIter += static_cast<genfit::KalmanFitStatus*>(fitTrack->getFitStatus(rep))->getNumIterations();

    // check if fit was successful
    //if (!fitter->isTrackFitted(fitTrack, rep)) {
    if (! fitTrack->getFitStatus(rep)->isFitted()) {
      std::cout << "Track could not be fitted successfully! \n";
      continue;
    }

    genfit::KalmanFittedStateOnPlane* kfsop = new genfit::KalmanFittedStateOnPlane(*(static_cast<genfit::KalmanFitterInfo*>(fitTrack->getPointWithMeasurement(0)->getFitterInfo(rep))->getBackwardUpdate()));
    if (debug) {
      std::cout << "state before extrapolating back to reference plane \n";
      kfsop->Print();
    }

    // extrapolate back to reference plane.
    try{
      rep->extrapolateToPlane(kfsop, stateRefOrig.getPlane());;
    }
    catch(genfit::Exception& e){
      std::cerr<<"Exception, next track"<<std::endl;
      std::cerr << e.what();
      continue; // here is a memleak!
    }

    rep->get6DStateCov(kfsop, stateFinal, covFinal);
    planeFinal = *stateRefOrig.getPlane();
    tResults->Fill();
    //fitTrack->Print();
    delete fitTrack;
    fitTrack = 0;
  }// end loop over events

  fOut->Write();
  delete fOut;

  fOut = TFile::Open(FILENAME, "READ");
  fOut->GetObject("tResults", tResults);
  TMatrixDSym *pMatrix = 0;
  tResults->SetBranchAddress("covFinal", &pMatrix);
  genfit::DetPlane *plane = 0;
  tResults->SetBranchAddress("planeFinal", &plane);
  tResults->SetBranchAddress("gfTrack", &fitTrack);

  for (Long_t nEntry = 0; nEntry < tResults->GetEntries(); ++nEntry) {
    tResults->GetEntry(nEntry);
    //fitTrack->Print();
    if (!fitTrack->checkConsistency())
      {
	std::cout << "stored track inconsistent" << std::endl;
	return 1;
      }
    genfit::AbsTrackRep *rep = fitTrack->getCardinalRep();

    genfit::KalmanFittedStateOnPlane* kfsop = new genfit::KalmanFittedStateOnPlane(*(static_cast<genfit::KalmanFitterInfo*>(fitTrack->getPointWithMeasurement(0)->getFitterInfo(rep))->getBackwardUpdate()));

    // extrapolate back to reference plane.
    try{
      rep->extrapolateToPlane(kfsop, genfit::SharedPlanePtr(plane));
    }
    catch(genfit::Exception& e){
      std::cerr<<"Exception, next track"<<std::endl;
      std::cerr << e.what();
    }

    rep->get6DStateCov(kfsop, stateFinal, covFinal);
    for (int i = 0; i < covFinal.GetNrows(); ++i)
      for (int j = 0; j < covFinal.GetNrows(); ++j)
	if ((covFinal - *pMatrix)(i,j) != 0) {
	  std::cout << "test failed.";
	  return 1;
	}
  }
  std::cout << "stored tracks are identical to fitted tracks, as far as tested." << std::endl;
  delete fitTrack;
  std::cout << "deleteing didn't segfault" << std::endl;
  return 0;
}
