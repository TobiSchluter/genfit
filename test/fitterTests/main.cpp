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
#include <KalmanFitterInfo.h>
#include <MaterialInfo.h>
#include <MeasuredStateOnPlane.h>
#include <MeasurementOnPlane.h>
#include <PlanarMeasurement.h>
#include <ProlateSpacePointMeasurement.h>
#include <RectangularFinitePlane.h>
#include <ReferenceStateOnPlane.h>
#include <SharedPlanePtr.h>
#include <SpacePointMeasurement.h>
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


#define VALGRIND

int main() {
  std::cerr<<"main"<<std::endl;

  const unsigned int nEvents = 1;
  const double BField = 15.;       // kGauss
  const double momentum = 0.2;     // GeV
  const double theta = 150;         // degree
  const double thetaDetPlane = 120;         // degree
  const double phiDetPlane = 0;         // degree
  const double pointDist = 5;      // cm; approx. distance between hits generated w/ RKTrackRep
  const double pointDistDeg = 20;      // degree; distance between hits generated w/ helix model
  const double resolution = 0.02;   // cm; resolution of generated hits

  const double resolutionWire = 5*resolution;   // cm; resolution of generated hits
  const TVector3 wireDir(0,0,1);
  const double skewAngle(5);
  const bool useSkew = false;
  const int nSuperLayer = 5;
  const double minDrift = 0;
  const double maxDrift = 2;
  const bool idealLRResolution = false; // resolve the l/r ambiguities of the wire hits

  const bool useDaf = true;

  const int pdg = 13;               // particle pdg code

  const bool smearPosMom = false;     // init the Reps with smeared pos and mom
  const double posSmear = 20*resolution;     // cm
  const double momSmear = 0.1*momentum;     // GeV
  const double zSmearFac = 100;

  const bool HelixTest = false;      // use helix for creating hits

  const bool matFX = true;         // include material effects; can only be disabled for RKTrackRep!
  const bool smoothing = true;

  const bool debug = false;

  // 0: PixMeasurement
  // 1: SpacePointMeasurement
  // 2: ProlateSpacePointMeasurement
  // 3: StripMeasurement
  // 4: WireMeasurement
  // 5: WirePointMeasurement
  std::vector<unsigned int> measurementTypes;

  measurementTypes.push_back(0);
  measurementTypes.push_back(0);
  measurementTypes.push_back(0);
  measurementTypes.push_back(0);
  measurementTypes.push_back(0);
  measurementTypes.push_back(0);




  // init fitters
  genfit::KalmanFitter kalman;


  gRandom->SetSeed(10);

  // init event display
#ifndef VALGRIND
  GenfitDisplay* display = GenfitDisplay::getInstance();
  display->reset();
#endif

  // init geometry and mag. field
  TGeoManager* geom = new TGeoManager("Geometry", "Geane geometry");
  TGeoManager::Import("genfitGeom.root");
  genfit::FieldManager::getInstance()->init(new genfit::ConstField(0.,0.,BField));
  genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());

  const double charge = TDatabasePDG::Instance()->GetParticle(pdg)->Charge()/(3.);


  // prepare output tree for Tracks
  genfit::Track* trueTrack = new genfit::Track();
  genfit::Track* fitTrack = new genfit::Track();
#ifndef VALGRIND
  // init rootapp (for drawing histograms)
  TApplication* rootapp = new TApplication("rootapp", 0, 0);
  TString outname = "out_Rep";
  outname += "_degPlane";
  outname += ".root";
  TFile *file = TFile::Open(outname,"RECREATE");
  TTree *tree = new TTree("t","Tracks");
  tree->Branch("trueTracks","Track",&trueTrack);
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
  TH1D *pVal = new TH1D("pVal","p-value",100,0.,1.);
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


      // trackrep for creating hits
      genfit::AbsTrackRep* rephits = new genfit::RKTrackRep(pos, mom, posErr, momErr, pdg);
      ((RKTrackRep*)rephits)->setPropDir(1);

      if (!matFX) MaterialEffects::getInstance()->setNoEffects();

      // remember original initial plane and state
      DetPlane referencePlane;
      TVectorT<double> referenceState(rephits->getState());

      // create smeared hits
      std::vector<AbsRecoHit*> hits;

      // true values for left right. 0 for non wire hits
      std::vector<int> leftRightTrue;

      TVector3 point, dir;
      int wireCounter = 0;
      if (debug) std::cerr << "Start creating hits ... \n";
      try{
        for (unsigned int i=0; i<measurementTypes.size(); ++i){
          // get current position and momentum
          if (!HelixTest) rephits->getPosMom(rephits->getReferencePlane(), point, dir);
          else{
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
          // create hit
          if (i==0){ // get reference state
            TVector3 planeNorm(dir);
            planeNorm.SetTheta(thetaDetPlane*TMath::Pi()/180);
            planeNorm.SetPhi(planeNorm.Phi()+phiDetPlane);
            TVector3 z(0,0,1);
            //z.SetTheta(thetaDetPlane*TMath::Pi()/180-TMath::PiOver2());
            referencePlane = DetPlane(point, planeNorm.Cross(z), (planeNorm.Cross(z)).Cross(planeNorm));
            if (!HelixTest) rephits->extrapolate(referencePlane, referenceState);
            else{
              double AtW = dir*referencePlane.getNormal();
              if (AtW<0) AtW *= -1.;
              referenceState[0] = charge/momentum;
              referenceState[1] = -1.*dir*referencePlane.getU()/AtW;
              referenceState[2] = -1.*dir*referencePlane.getV()/AtW;
              referenceState[3] = (point-referencePlane.getO())*referencePlane.getU();
              referenceState[4] = (point-referencePlane.getO())*referencePlane.getV();
              //std::cout<<"referenceState[2][0] "<<referenceState[2][0]<<"\n";
            }
          }

          AbsRecoHit* hit;

          TVector3 currentWireDir(wireDir);

          if (useSkew && (int)((double)wireCounter/(double)nSuperLayer)%2 == 1) {
            TVector3 perp(wireDir.Cross(dir));
            if ((int)((double)wireCounter/(double)nSuperLayer)%4 == 1){
              currentWireDir.Rotate(skewAngle*TMath::Pi()/180, wireDir.Cross(perp));
            }
            else currentWireDir.Rotate(-skewAngle*TMath::Pi()/180, wireDir.Cross(perp));
          }

          TVector3 planeNorm(dir);
          planeNorm.SetTheta(thetaDetPlane*TMath::Pi()/180);
          planeNorm.SetPhi(planeNorm.Phi()+phiDetPlane);
          TVector3 z(0,0,1);
          TVector3 x(1,0,0);

          int lr = 1;
          TVector3 wirePerp;
          if (measurementTypes[i] == 4 || measurementTypes[i] == 5){
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
            case 0: // 0: PixHit
              hit = new PixHit(point, planeNorm, planeNorm.Cross(z), resolution, true);
              //hit = new PixHit(point, z, x, resolution, true);
              break;

            case 1: // 1: SpacepointHit
              hit = new SpacepointHit(point, resolution, true);
              break;

            case 2: // 2: ProlateSpacepointHit
              hit = new ProlateSpacepointHit(point, currentWireDir, resolution, resolutionWire, true);
              break;

            case 3: // 3: StripHit
              hit = new StripHit(point, planeNorm, planeNorm.Cross(z), resolution, true);
              break;

            case 4: // 4: WireHit
              hit = new WireHit(point-wirePerp-currentWireDir, point-wirePerp+currentWireDir, wirePerp.Mag(), resolution, true);
              if (idealLRResolution){
                ((WireHit*)hit)->setLeftRightResolution(lr);
              }
              ++wireCounter;
              break;

            case 5: // 5: WirePointHit
              hit = new WirePointHit(point-wirePerp-currentWireDir, point-wirePerp+currentWireDir, wirePerp.Mag(), currentWireDir.Mag(), resolution, resolutionWire, true);
              if (idealLRResolution){
                ((WireHit*)hit)->setLeftRightResolution(lr);
              }
              ++wireCounter;
              break;

            default:
              std::cerr << "hit type not defined!" << std::endl;
              exit(0);
          }
          hits.push_back(hit);

          if (debug) {std::cout << "(smeared) hit coordinates"; hit->getRawHitCoord().Print();}

          if (!HelixTest) {
            // stepalong (approximately)
            dir.SetMag(pointDist);
            DetPlane pl(point+dir, dir);
            rephits->extrapolate(pl);
          }
        }

        assert(measurementTypes.size() == leftRightTrue.size());
      }
      catch(Exception& e){
        std::cerr<<"Exception, next track"<<std::endl;
        e.what();
        continue; // here is a memleak!
      }

      if (debug) std::cerr << "... done creating hits \n";



      // trackrep to be fitted and tested
      AbsTrackRep* rep;
      rep = new RKTrackRep(posM, momM, posErr, momErr, pdg);

      // create track
      if (fitTrack != NULL) delete fitTrack;
      fitTrack = new Track(rep); //initialized with smeared rep
      //fitTrack->addTrackRep(rep->clone()); // check if everything works fine with more than one rep
      fitTrack->setSmoothing(smoothing);

      // add hits
      for(unsigned int i=0; i<hits.size(); ++i){
        fitTrack->addHit(hits[i],
                         measurementTypes[i], //detector id
                         i); // hit id
      }

      // print trackCand
      if (debug) fitTrack->getCand().Print();



      // do the fit
      try{
        if (useDaf) {
          if (debug) std::cerr<<"Starting the fitter (Daf)"<<std::endl;
          daf.processTrack(fitTrack);
        }
        else {
          if (debug) std::cerr<<"Starting the fitter (Kalman)"<<std::endl;
          kalman.processTrack(fitTrack);
        }
        if (debug) std::cerr<<"fitter is finished!"<<std::endl;
      }
      catch(Exception& e){
        e.what();
        std::cerr<<"Exception, next track"<<std::endl;
        continue;
      }

      if (debug) {
        fitTrack->getBK(0)->Print();
      }


      //choose trackrep to check
      AbsTrackRep* repCheck = fitTrack->getTrackRep(0);

      // check if fit was successfull
      if(repCheck->getStatusFlag() != 0 ) {
        continue;
      }


#ifndef VALGRIND
      // add track to event display
      std::vector<Track*> event;
      event.push_back(fitTrack);
      display->addEvent(event);
#endif

      if (debug) {
        std::cout << "cov before extrapolating back to reference plane \n";
        repCheck->getCov().Print();

        std::cout << "pos before extrapolating back to reference plane \n";
        repCheck->getPos().Print();
      }

      // extrapolate back to reference plane.
      try{
        repCheck->extrapolate(referencePlane);
      }
      catch(Exception& e){
        std::cerr<<"Exception, next track"<<std::endl;
        e.what();
        continue; // here is a memleak!
      }


      // check getPos etc methods
      const double epsilon = 1E-3;
      TVector3 getPos, getMom;


      // calculate pulls
      TVectorT<double> state(repCheck->getState());
      TMatrixTSym<double> cov(repCheck->getCov());


#ifndef VALGRIND
      hmomRes->Fill( (charge/state[0]-momentum));
      hupRes->Fill(  (state[1]-referenceState[1]));
      hvpRes->Fill(  (state[2]-referenceState[2]));
      huRes->Fill(   (state[3]-referenceState[3]));
      hvRes->Fill(   (state[4]-referenceState[4]));

      if (cov[0][0]>0) {
        hqopPu->Fill( (state[0]-referenceState[0]) / sqrt(cov[0][0]) );
        pVal->Fill(   repCheck->getPVal());
        hupPu->Fill(  (state[1]-referenceState[1]) / sqrt(cov[1][1]) );
        hvpPu->Fill(  (state[2]-referenceState[2]) / sqrt(cov[2][2]) );
        huPu->Fill(   (state[3]-referenceState[3]) / sqrt(cov[3][3]) );
        hvPu->Fill(   (state[4]-referenceState[4]) / sqrt(cov[4][4]) );
      }

      // print covariance
      if (debug) cov.Print();
      else if((iEvent)%100==0)  cov.Print();



      // check l/r resolution
      if (smoothing && useDaf) {
        for (unsigned int i=0; i<leftRightTrue.size(); ++i){
          int trueSide = leftRightTrue[i];
          if (trueSide == 0) continue; // not a wire hit
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


      if (debug && smoothing){
        std::cout << "smoothed positions \n";
        DetPlane plane;
        for(unsigned int j = 0; j < measurementTypes.size(); j++) { // loop over all hits in the track
          TVectorT<double> state;
          TMatrixTSym<double> cov;
          TMatrixT<double> auxInfo;
          Tools::getBiasedSmoothedData(fitTrack, 0, j, state, cov, plane, auxInfo);
          repCheck->setData(state, plane, &cov, &auxInfo);
          repCheck->getPos(plane).Print();
        }
      }

      if (debug) std::cerr<<"Fill Tree ..." << std::endl;
      tree->Fill();
#endif


  }// end loop over events

  std::cerr<<"maxWeight = " << maxWeight << std::endl;

#ifndef VALGRIND
  if (debug) std::cerr<<"Write Tree ...";
  tree->Write();
  if (debug) std::cerr<<"... done"<<std::endl;

  if (debug) std::cerr<<"Draw histograms ...";
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

  if (debug) std::cerr<<"... done"<<std::endl;

  // open event display
  display->setOptions("THDPMAG"); // G show geometry
  display->open();

  rootapp->Run();

  file->Close();
#endif

  if (debug) std::cerr<<"... closed file"<<std::endl;




