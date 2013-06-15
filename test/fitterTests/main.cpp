#include <iostream>
#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>

#include <AbsFinitePlane.h>
#include <AbsFitterInfo.h>
#include <AbsMeasurement.h>
#include <AbsTrackRep.h>
#include <ConstField.h>
#include <DetPlane.h>
#include <Exception.h>
#include <FieldManager.h>
#include <KalmanFittedStateOnPlane.h>
#include <KalmanFitter.h>
#include <KalmanFitterRefTrack.h>
#include <KalmanFitterInfo.h>
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

#include <EventDisplay.h>

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

#include <memory>


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

//---------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------


//#define VALGRIND

int main() {
  std::cout<<"main"<<std::endl;

  const unsigned int nEvents = 2;
  const double BField = 15.;       // kGauss
  const double momentum = 0.1;     // GeV
  const double theta = 150;         // degree
  const double thetaDetPlane = 90;         // degree
  const double phiDetPlane = 0;         // degree
  const double pointDist = 5;      // cm; approx. distance between measurements generated w/ RKTrackRep
  const double pointDistDeg = 5;      // degree; distance between measurements generated w/ helix model
  const double resolution = 0.02;   // cm; resolution of generated measurements

  const double resolutionWire = 5*resolution;   // cm; resolution of generated measurements
  const TVector3 wireDir(0,0,1);
  const double skewAngle(5);
  const bool useSkew = false;
  const int nSuperLayer = 5;
  const double minDrift = 0;
  const double maxDrift = 2;
  const bool idealLRResolution = false; // resolve the l/r ambiguities of the wire measurements

  const int fitterId = 1; // 1 = SimpleKalman; 2 = KalmanFitterRefTrack; 3 = DAF

  const int pdg = 13;               // particle pdg code

  const bool smearPosMom = false;     // init the Reps with smeared pos and mom
  const double posSmear = 20*resolution;     // cm
  const double momSmear = 0.1*momentum;     // GeV
  const double zSmearFac = 10;

  const bool HelixTest = false;      // use helix for creating measurements

  const bool matFX = false;         // include material effects; can only be disabled for RKTrackRep!

  const bool debug = true;

  enum eMeasurementType { Pixel = 0,
			  Spacepoint,
			  ProlateSpacepoint,
			  Strip,
			  Wire,
			  WirePoint, };
  std::vector<unsigned int> measurementTypes;
  for (int i = 0; i < 4; ++i)
    measurementTypes.push_back(Pixel);



  // init fitters
  genfit::KalmanFitter simpleKalman(6);
  genfit::KalmanFitterRefTrack kalmanFitterRefTrack(1);


  gRandom->SetSeed(10);

  // init event display
#ifndef VALGRIND
  genfit::EventDisplay* display = genfit::EventDisplay::getInstance();
  display->reset();
#endif

  signal(SIGSEGV, handler);   // install our handler

  // init geometry and mag. field
  TGeoManager* geom = new TGeoManager("Geometry", "Geane geometry");
  TGeoManager::Import("genfitGeom.root");
  genfit::FieldManager::getInstance()->init(new genfit::ConstField(0.,0.,BField));
  genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());

  const double charge = TDatabasePDG::Instance()->GetParticle(pdg)->Charge()/(3.);


  // prepare output tree for Tracks
  std::unique_ptr<genfit::Track> fitTrack(new genfit::Track());
#ifndef VALGRIND
  // init rootapp (for drawing histograms)
  TApplication* rootapp = new TApplication("rootapp", 0, 0);
  TString outname = "out_Rep";
  outname += "_degPlane";
  outname += ".root";
  TFile *file = TFile::Open(outname,"RECREATE");
  TTree *tree = new TTree("t","Tracks");
  tree->Branch("fitTracks","Track",&fitTrack);


  // create histograms
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetOptFit(1111);

  TH1D *hmomRes = new TH1D("hmomRes","mom res",2000,-0.01*momentum,0.01*momentum);
  TH1D *hupRes = new TH1D("hupRes","u' res",2000,-0.05,0.05);
  TH1D *hvpRes = new TH1D("hvpRes","v' res",2000,-0.05,0.05);
  TH1D *huRes = new TH1D("huRes","u res",2000,-0.05,0.05);
  TH1D *hvRes = new TH1D("hvRes","v res",2000,-0.05,0.05);

  TH1D *hqopPu = new TH1D("hqopPu","q/p pull",200,-6.,6.);
  TH1D *pVal = new TH1D("pVal","p-value",100,0.,1.00000001);
  TH1D *hupPu = new TH1D("hupPu","u' pull",200,-6.,6.);
  TH1D *hvpPu = new TH1D("hvpPu","v' pull",200,-6.,6.);
  TH1D *huPu = new TH1D("huPu","u pull",200,-6.,6.);
  TH1D *hvPu = new TH1D("hvPu","v pull",200,-6.,6.);

  TH1D *weights = new TH1D("weights","wire Daf vs true weights",500,-1.01,1.01);
#endif

  double maxWeight(0);

  // main loop
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



      TVector3 posErr(1.,1.,1.);
      posErr *= posSmear;
      TVector3 momErr(1.,1.,1.);
      momErr *= momSmear;

      // smeared start values
      TVector3 posM(pos);
      TVector3 momM(mom);
      if (smearPosMom) {
        posM.SetX(gRandom->Gaus(posM.X(),posSmear));
        posM.SetY(gRandom->Gaus(posM.Y(),posSmear));
        posM.SetZ(gRandom->Gaus(posM.Z(),zSmearFac*posSmear));

        momM.SetX(gRandom->Gaus(momM.X(),momSmear));
        momM.SetY(gRandom->Gaus(momM.Y(),momSmear));
        momM.SetZ(gRandom->Gaus(momM.Z(),momSmear));
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
      const genfit::StateOnPlane stateRefOrig(stateSmeared);

      // create smeared measurements
      std::vector<genfit::AbsMeasurement*> measurements;

      // true values for left right. 0 for non wire measurements
      std::vector<int> leftRightTrue;

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
          if (measurementTypes[i] == Wire
	      || measurementTypes[i] == WirePoint){
            wirePerp = dir.Cross(currentWireDir);
            if (gRandom->Uniform(-1,1) >= 0) {
              wirePerp *= -1.;
              lr = -1;
            }
            wirePerp.SetMag(gRandom->Uniform(minDrift, maxDrift));

            leftRightTrue.push_back(lr);
          }
          else leftRightTrue.push_back(0);

          switch(measurementTypes[i]){
            case Pixel: {
              if (debug) std::cerr << "create PixHit" << std::endl;

              genfit::SharedPlanePtr plane(new genfit::DetPlane(point, planeNorm.Cross(z), (planeNorm.Cross(z)).Cross(planeNorm)));

              TVectorD hitCoords(2);
              hitCoords(0) = gRandom->Gaus(0,resolution);
              hitCoords(1) = gRandom->Gaus(0,resolution);

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
              hitCoords(0) = gRandom->Gaus(point.X(),resolution);
              hitCoords(1) = gRandom->Gaus(point.Y(),resolution);
              hitCoords(2) = gRandom->Gaus(point.Z(),resolution);

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
              hitCoords(0) = point.X();
              hitCoords(1) = point.Y();
              hitCoords(2) = point.Z();

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
              hitCoords(0) += gRandom->Gaus(0, smearVec(0));
              hitCoords(1) += gRandom->Gaus(0, smearVec(1));
              hitCoords(2) += gRandom->Gaus(0, smearVec(2));


              // rotate cov
              hitCov.Similarity(rot);

              measurement = new genfit::ProlateSpacepointMeasurement(hitCoords, hitCov, 2, i, nullptr);

              static_cast<genfit::ProlateSpacepointMeasurement*>(measurement)->setLargestErrorDirection(currentWireDir);
            }
            break;

            case Strip: {
              if (debug) std::cerr << "create StripHit" << std::endl;

              genfit::SharedPlanePtr plane(new genfit::DetPlane(point, planeNorm.Cross(z), (planeNorm.Cross(z)).Cross(planeNorm)));

              TVectorD hitCoords(1);
              hitCoords(0) = gRandom->Gaus(0,resolution);

              TMatrixDSym hitCov(1);
              hitCov(0,0) = resolution*resolution;

              measurement = new genfit::PlanarMeasurement(hitCoords, hitCov, 3, i, nullptr);
              static_cast<genfit::PlanarMeasurement*>(measurement)->setPlane(plane, i);
            }
            break;

            case Wire: {
              if (debug) std::cerr << "create WireHit" << std::endl;

              TVectorD hitCoords(7);
              hitCoords(0) = (point-wirePerp-currentWireDir).X();
              hitCoords(1) = (point-wirePerp-currentWireDir).Y();
              hitCoords(2) = (point-wirePerp-currentWireDir).Z();

              hitCoords(3) = (point-wirePerp+currentWireDir).X();
              hitCoords(4) = (point-wirePerp+currentWireDir).Y();
              hitCoords(5) = (point-wirePerp+currentWireDir).Z();

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

              TVectorD hitCoords(8);
              hitCoords(0) = (point-wirePerp-currentWireDir).X();
              hitCoords(1) = (point-wirePerp-currentWireDir).Y();
              hitCoords(2) = (point-wirePerp-currentWireDir).Z();

              hitCoords(3) = (point-wirePerp+currentWireDir).X();
              hitCoords(4) = (point-wirePerp+currentWireDir).Y();
              hitCoords(5) = (point-wirePerp+currentWireDir).Z();

              hitCoords(6) = gRandom->Gaus(wirePerp.Mag(), resolution);
              hitCoords(7) = gRandom->Gaus(currentWireDir.Mag(), resolutionWire);

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

          if (debug) {
            std::cout << "(smeared) measurement coordinates"; measurement->getRawHitCoords().Print();
          }

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
        e.what();
        continue; // here is a memleak!
      }

      if (debug) std::cout << "... done creating measurements \n";



      // create track
      std::cout << "fitTrack before resetting"; fitTrack->Print();
      fitTrack.reset(new genfit::Track(rep, rep->get6DState(&stateSmeared))); //initialized with smeared rep

      //fitTrack->addTrackRep(rep->clone()); // check if everything works fine with more than one rep

      // add measurements
      for(unsigned int i=0; i<measurements.size(); ++i){
        std::vector<genfit::AbsMeasurement*> measVec;
        measVec.push_back(measurements[i]);
        fitTrack->insertPoint(new genfit::TrackPoint(measVec, fitTrack.get()));
      }

      // print trackCand
      //if (debug) fitTrack->Print();



      // do the fit
      try{
        switch (fitterId) {
          case 1:
            if (debug) std::cout<<"Starting the fitter (simple Kalman)"<<std::endl;
            simpleKalman.processTrack(fitTrack.get(), rep);
            break;

          case 2:
            if (debug) std::cout<<"Starting the fitter (reference track Kalman)"<<std::endl;
            kalmanFitterRefTrack.processTrack(fitTrack.get(), rep);
            break;

          case 3:
            //if (debug) std::cout<<"Starting the fitter (DAF)"<<std::endl;
            //daf.processTrack.processTrack(fitTrack, rep);
            //break;
          default:
            std::cout<<"no fitter selected!"<<std::endl;

        }

        if (debug) std::cout<<"fitter is finished!"<<std::endl;
      }
      catch(genfit::Exception& e){
        std::cerr << e.what();
        std::cerr << "Exception, next track" << std::endl;
        continue;
      }

      if (debug) fitTrack->Print();


      //choose trackrep to check
      genfit::AbsTrackRep* repCheck = fitTrack->getTrackRep(0);

      // FIXME
      // check if fit was successfull
      /*if(repCheck->getStatusFlag() != 0 ) {
        continue;
      }*/


#ifndef VALGRIND
      // add track to event display
      std::vector<genfit::Track*> event;
      event.push_back(new genfit::Track(*fitTrack)); // FIXME memleak
      display->addEvent(event);
#endif

      genfit::KalmanFittedStateOnPlane* kfsop = new genfit::KalmanFittedStateOnPlane(*(static_cast<genfit::KalmanFitterInfo*>(fitTrack->getPointWithMeasurement(0)->getFitterInfo(repCheck))->getBackwardUpdate()));
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
        e.what();
        continue; // here is a memleak!
      }





#ifndef VALGRIND
      // calculate pulls
      const TVectorD& referenceState = stateRefOrig.getState();

      const TVectorD& state = kfsop->getState();
      const TMatrixDSym& cov = kfsop->getCov();

      double pval = kalmanFitterRefTrack.getPVal(fitTrack.get(), repCheck); // FIXME choose fitter that has been used

      hmomRes->Fill( (charge/state[0]-momentum));
      hupRes->Fill(  (state[1]-referenceState[1]));
      hvpRes->Fill(  (state[2]-referenceState[2]));
      huRes->Fill(   (state[3]-referenceState[3]));
      hvRes->Fill(   (state[4]-referenceState[4]));

      hqopPu->Fill( (state[0]-referenceState[0]) / sqrt(cov[0][0]) );
      pVal->Fill(   pval);
      hupPu->Fill(  (state[1]-referenceState[1]) / sqrt(cov[1][1]) );
      hvpPu->Fill(  (state[2]-referenceState[2]) / sqrt(cov[2][2]) );
      huPu->Fill(   (state[3]-referenceState[3]) / sqrt(cov[3][3]) );
      hvPu->Fill(   (state[4]-referenceState[4]) / sqrt(cov[4][4]) );

      // print covariance
      if (debug) cov.Print();
      else if((iEvent)%100==0)  cov.Print();



      // check l/r resolution
     /* if (fitterId == 3) {
        for (unsigned int i=0; i<leftRightTrue.size(); ++i){
          int trueSide = leftRightTrue[i];
          if (trueSide == 0) continue; // not a wire measurement
          const TVectorT<double>& dafWeightLR = fitTrack->getBK(0)->getVector(BKKey_dafWeight, i);
          assert(dafWeightLR.GetNrows()==2);

          double weightCorrectSide, weightWrongSide;

          if (trueSide < 0) {
            weightCorrectSide = dafWeightLR[0];
            weightWrongSide =  dafWeightLR[1];
          }
          else {
            weightCorrectSide = dafWeightLR[1];
            weightWrongSide =  dafWeightLR[0];
          }
          weightWrongSide -= 1.;

          weights->Fill(weightCorrectSide);
          weights->Fill(weightWrongSide);

          if (weightCorrectSide>maxWeight) maxWeight = weightCorrectSide;
        }

      }
*/

      if (debug) std::cerr<<"Fill Tree ..." << std::endl;
      tree->Fill();
#endif


  }// end loop over events

  std::cout<<"maxWeight = " << maxWeight << std::endl;

#ifndef VALGRIND
  if (debug) std::cout<<"Write Tree ...";
  tree->Write();
  if (debug) std::cout<<"... done"<<std::endl;

  if (debug) std::cout<<"Draw histograms ...";
  // fit and draw histograms
  TCanvas* c1 = new TCanvas();
  c1->Divide(2,3);

  c1->cd(1);
  hmomRes->Fit("gaus");
  hmomRes->Draw();

  c1->cd(2);
  weights->Draw();

  c1->cd(3);
  hupRes->Fit("gaus");
  hupRes->Draw();

  c1->cd(4);
  hvpRes->Fit("gaus");
  hvpRes->Draw();

  c1->cd(5);
  huRes->Fit("gaus");
  huRes->Draw();

  c1->cd(6);
  hvRes->Fit("gaus");
  hvRes->Draw();

  c1->Write();

  TCanvas* c2 = new TCanvas();
  c2->Divide(2,3);

  c2->cd(1);
  hqopPu->Fit("gaus");
  hqopPu->Draw();

  c2->cd(2);
  pVal->Fit("pol1");
  pVal->Draw();
  c2->cd(3);
  hupPu->Fit("gaus");
  hupPu->Draw();

  c2->cd(4);
  hvpPu->Fit("gaus");
  hvpPu->Draw();

  c2->cd(5);
  huPu->Fit("gaus");
  huPu->Draw();

  c2->cd(6);
  hvPu->Fit("gaus");
  hvPu->Draw();

  c2->Write();

  if (debug) std::cout<<"... done"<<std::endl;

  // open event display
  display->setOptions("THDPMA"); // G show geometry
  display->open();

  //rootapp->Run();

  //file->Close();
#endif

  if (debug) std::cout<<"... closed file"<<std::endl;

}


