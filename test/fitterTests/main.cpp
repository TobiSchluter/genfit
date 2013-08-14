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
#include <AbsKalmanFitter.h>
#include <KalmanFitter.h>
#include <KalmanFitterRefTrack.h>
#include <KalmanFitterInfo.h>
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

#include <EventDisplay.h>

#include <HelixTrackModel.h>
#include <MeasurementCreator.h>

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

#ifdef VALGRIND
  #include <valgrind/callgrind.h>
#else
#define CALLGRIND_START_INSTRUMENTATION
#define CALLGRIND_STOP_INSTRUMENTATION
#define CALLGRIND_DUMP_STATS
#endif

int main() {
  std::cout<<"main"<<std::endl;

  enum eFitterType { SimpleKalman = 0,
        RefKalman,
        DafSimple,
        DafRef};

  const unsigned int nEvents = 1000;
  const double BField = 15.;       // kGauss
  const double momentum = 0.4;     // GeV
  const double theta = 120;         // degree
  const double thetaDetPlane = 90;         // degree
  const double phiDetPlane = 0;         // degree
  const double pointDist = 5;      // cm; approx. distance between measurements
  const double resolution = 0.02;   // cm; resolution of generated measurements

  const double resolutionWire = 5*resolution;   // cm; resolution of generated measurements
  const TVector3 wireDir(0,0,1);
  const double skewAngle(5);
  const bool useSkew = true;
  const int nSuperLayer = 5;
  const double minDrift = 0.;
  const double maxDrift = 2;
  const bool idealLRResolution = false; // resolve the l/r ambiguities of the wire measurements

  const double outlierProb = 0.1;
  const double outlierRange = 5;

  const double hitSwitchProb = 0.1; // probability to give hits to fit in wrong order (flip two hits)

  const int splitTrack = 4; // for track merging testing:

  //const eFitterType fitterId = SimpleKalman;
  //const eFitterType fitterId = RefKalman;
  const eFitterType fitterId = DafRef;
  const genfit::eMultipleMeasurementHandling mmHandling = genfit::weightedAverage;
  //const genfit::eMultipleMeasurementHandling mmHandling = genfit::unweightedClosestToReference;
  //const genfit::eMultipleMeasurementHandling mmHandling = genfit::unweightedClosestToPrediction;
  const int nIter = 10; // max number of iterations
  const double dPVal = 1.E-3; // convergence criterion

  const bool resort = true;

  const bool twoReps = true; // test if everything works with more than one rep in the tracks

  const int pdg = 13;               // particle pdg code

  const bool smearPosMom = true;     // init the Reps with smeared pos and mom
  const double chargeSwitchProb = 0.1; // probability to seed with wrong charge sign
  const double posSmear = 10*resolution;     // cm
  const double momSmear = 5. /180.*TMath::Pi();     // rad
  const double momMagSmear = 0.2;   // relative
  const double zSmearFac = 2;


  const bool matFX = false;         // include material effects; can only be disabled for RKTrackRep!

  const bool debug = false;

  std::vector<genfit::eMeasurementType> measurementTypes;

  /*measurementTypes.push_back(genfit::Pixel);
  measurementTypes.push_back(genfit::Spacepoint);
  measurementTypes.push_back(genfit::ProlateSpacepoint);
  measurementTypes.push_back(genfit::StripV);
  measurementTypes.push_back(genfit::StripU);
  measurementTypes.push_back(genfit::Wire);
  measurementTypes.push_back(genfit::WirePoint);*/
  for (int i = 0; i < 20; ++i)
    measurementTypes.push_back(genfit::eMeasurementType(gRandom->Uniform(7)));



  gRandom->SetSeed(14);
  signal(SIGSEGV, handler);   // install our handler

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


  // init MeasurementCreator
  genfit::MeasurementCreator measurementCreator;
  measurementCreator.setResolution(resolution);
  measurementCreator.setResolutionWire(resolutionWire);
  measurementCreator.setOutlierProb(outlierProb);
  measurementCreator.setOutlierRange(outlierRange);
  measurementCreator.setThetaDetPlane(thetaDetPlane);
  measurementCreator.setPhiDetPlane(phiDetPlane);
  measurementCreator.setWireDir(wireDir);
  measurementCreator.setMinDrift(minDrift);
  measurementCreator.setMaxDrift(maxDrift);
  measurementCreator.setIdealLRResolution(idealLRResolution);
  measurementCreator.setUseSkew(useSkew);
  measurementCreator.setSkewAngle(skewAngle);
  measurementCreator.setNSuperLayer(nSuperLayer);
  measurementCreator.setDebug(debug);


  // init event display
#ifndef VALGRIND
  genfit::EventDisplay* display = genfit::EventDisplay::getInstance();
  display->reset();
#endif


  // init geometry and mag. field
  TGeoManager* geom = new TGeoManager("Geometry", "Geane geometry");
  TGeoManager::Import("genfitGeom.root");
  genfit::FieldManager::getInstance()->init(new genfit::ConstField(0.,0.,BField));
  genfit::FieldManager::getInstance()->useCache(true, 8);
  genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());

  const double charge = TDatabasePDG::Instance()->GetParticle(pdg)->Charge()/(3.);


  // prepare output tree for Tracks 
  // std::unique_ptr<genfit::Track> fitTrack(new genfit::Track());
  genfit::Track* fitTrack = new genfit::Track();
  genfit::Track* secondTrack = new genfit::Track();
#ifndef VALGRIND
  // init rootapp (for drawing histograms)
  TApplication* rootapp = new TApplication("rootapp", 0, 0);
  /*TString outname = "out_Rep";
  outname += "_degPlane";
  outname += ".root";
  TFile *file = TFile::Open(outname,"RECREATE");
  TTree *tree = new TTree("t","Tracks");
  tree->Branch("fitTracks","Track", fitTrack->get());*/


  // create histograms
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetOptFit(1111);

  double s = measurementTypes.size();
  TH1D *hmomRes = new TH1D("hmomRes","mom res",500,-2*resolution*momentum/s,2*resolution*momentum/s);
  TH1D *hupRes = new TH1D("hupRes","u' res",500,-5*resolution/s, 5*resolution/s);
  TH1D *hvpRes = new TH1D("hvpRes","v' res",500,-5*resolution/s, 5*resolution/s);
  TH1D *huRes = new TH1D("huRes","u res",500,-5*resolution, 5*resolution);
  TH1D *hvRes = new TH1D("hvRes","v res",500,-5*resolution, 5*resolution);

  TH1D *hqopPu = new TH1D("hqopPu","q/p pull",200,-6.,6.);
  TH1D *pVal = new TH1D("pVal","p-value",100,0.,1.00000001);
  pVal->SetMinimum(0);
  TH1D *hupPu = new TH1D("hupPu","u' pull",200,-6.,6.);
  TH1D *hvpPu = new TH1D("hvpPu","v' pull",200,-6.,6.);
  TH1D *huPu = new TH1D("huPu","u pull",200,-6.,6.);
  TH1D *hvPu = new TH1D("hvPu","v pull",200,-6.,6.);

  TH1D *weights = new TH1D("weights","Daf vs true weights",500,-1.01,1.01);

  TH1D *trackLenRes = new TH1D("trackLenRes","(trueLen - FittedLen) / trueLen",500,-0.01,0.01);
#endif

  double maxWeight(0);
  unsigned int nTotalIter(0);
  unsigned int nTotalIterSecond(0);
  unsigned int nSuccessfullFits(0);
  unsigned int nSuccessfullFitsSecond(0);

  CALLGRIND_START_INSTRUMENTATION;

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

      if (debug) {
        std::cout << "start values \n";
        pos.Print();
        mom.Print();
      }

      // calc helix parameters
      genfit::HelixTrackModel* helix = new genfit::HelixTrackModel(pos, mom, charge);
      measurementCreator.setTrackModel(helix);

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
      double sign(1.);
      if (chargeSwitchProb > gRandom->Uniform(1.))
        sign = -1.;
      genfit::AbsTrackRep* rep = new genfit::RKTrackRep(sign*pdg);
      sign = 1.;
      if (chargeSwitchProb > gRandom->Uniform(1.))
        sign = -1.;
      genfit::AbsTrackRep* secondRep = new genfit::RKTrackRep(sign*-211);
      genfit::StateOnPlane stateRef(rep);
      rep->setPosMom(stateRef, pos, mom);

      // smeared start state
      genfit::StateOnPlane stateSmeared(rep);
      rep->setPosMom(stateSmeared, posM, momM);

      //rep->setPropDir(1);

      if (!matFX) genfit::MaterialEffects::getInstance()->setNoEffects();

      // remember original initial state
      const genfit::StateOnPlane stateRefOrig(stateRef);

      // create smeared measurements
      std::vector<genfit::AbsMeasurement*> measurements;

      std::vector<bool> outlierTrue;
      bool outlier;
      // true values for left right. 0 for non wire measurements
      std::vector<int> leftRightTrue;
      int lr;

      double trueLen;

      try{
        for (unsigned int i=0; i<measurementTypes.size(); ++i){
          trueLen = i*pointDist;

          measurements.push_back(measurementCreator.create(measurementTypes[i], trueLen, outlier, lr));
          outlierTrue.push_back(outlier);
          leftRightTrue.push_back(lr);
        }
        assert(measurementTypes.size() == leftRightTrue.size());
        assert(measurementTypes.size() == outlierTrue.size());
      }
      catch(genfit::Exception& e){
        std::cerr<<"Exception, next track"<<std::endl;
        std::cerr << e.what();
        continue; // here is a memleak!
      }

      if (debug) std::cout << "... done creating measurements \n";



      // create track
      fitTrack = new genfit::Track(rep, rep->get6DState(stateSmeared)); //initialized with smeared rep
      secondTrack = new genfit::Track(rep, rep->get6DState(stateSmeared)); //initialized with smeared rep
      if (twoReps) {
        fitTrack->addTrackRep(secondRep);
        secondTrack->addTrackRep(secondRep);
      }
      //if (debug) fitTrack->Print("C");

      assert(fitTrack->checkConsistency());
      //fitTrack->addTrackRep(rep->clone()); // check if everything works fine with more than one rep

      // add measurements
      for(unsigned int i=0; i<measurements.size(); ++i){
        if (splitTrack > 0 && i >= splitTrack)
          break;
        std::vector<genfit::AbsMeasurement*> measVec;
        measVec.push_back(measurements[i]);
        if (i>0 && hitSwitchProb > gRandom->Uniform(1.))
          fitTrack->insertPoint(new genfit::TrackPoint(measVec, fitTrack), -2);
        else
          fitTrack->insertPoint(new genfit::TrackPoint(measVec, fitTrack));

        assert(fitTrack->checkConsistency());
        //if (debug) fitTrack->Print("C");
      }

      if (splitTrack > 0) {
        for(unsigned int i=splitTrack; i<measurements.size(); ++i){
          std::vector<genfit::AbsMeasurement*> measVec;
          measVec.push_back(measurements[i]);
          if (i>0 && hitSwitchProb > gRandom->Uniform(1.))
            secondTrack->insertPoint(new genfit::TrackPoint(measVec, fitTrack), -2);
          else
            secondTrack->insertPoint(new genfit::TrackPoint(measVec, fitTrack));

          //if (debug) fitTrack->Print("C");
        }
      }

      assert(fitTrack->checkConsistency());
      assert(secondTrack->checkConsistency());

      //if (debug) fitTrack->Print();

      // do the fit
      try{
        if (debug) std::cout<<"Starting the fitter"<<std::endl;
        fitter->processTrack(fitTrack, resort);
        if (splitTrack > 0)
          fitter->processTrack(secondTrack, resort);

        if (debug) std::cout<<"fitter is finished!"<<std::endl;
      }
      catch(genfit::Exception& e){
        std::cerr << e.what();
        std::cerr << "Exception, next track" << std::endl;
        continue;
      }

      if (splitTrack > 0) {
        if (debug) fitTrack->Print("C");
        if (debug) secondTrack->Print("C");

        fitTrack->mergeTrack(secondTrack);

        if (debug) fitTrack->Print("C");

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
      }


      if (debug)
        fitTrack->Print("C");

      assert(fitTrack->checkConsistency());
      assert(secondTrack->checkConsistency());

#ifndef VALGRIND
      // add track to event display
      std::vector<genfit::Track*> event;
      event.push_back(fitTrack);
      display->addEvent(event);
#endif


      nTotalIter += static_cast<genfit::KalmanFitStatus*>(fitTrack->getFitStatus(rep))->getNumIterations();
      if (fitTrack->getFitStatus(rep)->isFitConverged())
        nSuccessfullFits += 1;

      if (twoReps) {
        nTotalIterSecond += static_cast<genfit::KalmanFitStatus*>(fitTrack->getFitStatus(secondRep))->getNumIterations();
        if (fitTrack->getFitStatus(secondRep)->isFitConverged())
          nSuccessfullFitsSecond += 1;
      }


      // check if fit was successful
      //if (!fitter->isTrackFitted(fitTrack, rep)) {
      if (! fitTrack->getFitStatus(rep)->isFitConverged()) {
        std::cout << "Track could not be fitted successfully! \n";
        continue;
      }


      genfit::TrackPoint* tp = fitTrack->getPointWithMeasurementAndFitterInfo(0, rep);
      if (tp == NULL) {
        std::cout << "Track has no TrackPoint with fitterInfo! \n";
        continue;
      }
      genfit::KalmanFittedStateOnPlane* kfsop = new genfit::KalmanFittedStateOnPlane(*(static_cast<genfit::KalmanFitterInfo*>(tp->getFitterInfo(rep))->getBackwardUpdate()));
      if (debug) {
        std::cout << "state before extrapolating back to reference plane \n";
        kfsop->Print();
      }

      // extrapolate back to reference plane.
      try{
        rep->extrapolateToPlane(*kfsop, stateRefOrig.getPlane());;
      }
      catch(genfit::Exception& e){
        std::cerr<<"Exception, next track"<<std::endl;
        std::cerr << e.what();
        continue; // here is a memleak!
      }

#ifndef VALGRIND
      // calculate pulls
      const TVectorD& referenceState = stateRefOrig.getState();

      const TVectorD& state = kfsop->getState();
      const TMatrixDSym& cov = kfsop->getCov();

      double pval = fitter->getPVal(fitTrack, rep);
      assert( fabs(pval - static_cast<genfit::KalmanFitStatus*>(fitTrack->getFitStatus(rep))->getBackwardPVal()) < 1E-10 );

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

      try {
        trackLenRes->Fill( (trueLen - fitTrack->getTrackLen(rep)) / trueLen );

        if (debug) {
          std::cout << "true track length = " << trueLen << "; fitted length = " << fitTrack->getTrackLen(rep) << "\n";
          std::cout << "fitted tof = " << fitTrack->getTOF(rep) << " ns\n";
        }
      }
      catch (genfit::Exception& e) {
        std::cerr << e.what();
        std::cout << "could not get TraclLen or TOF! \n";
      }

      // print covariance
      if (debug) cov.Print();
      else if((iEvent)%100==0)  cov.Print();



      // check l/r resolution and outlier rejection
      if (dynamic_cast<genfit::DAF*>(fitter) != NULL) {
        for (unsigned int i=0; i<leftRightTrue.size(); ++i){

          if (! fitTrack->getPointWithMeasurement(i)->hasFitterInfo(rep))
            continue;

          if (debug) {
            std::vector<double> dafWeights = dynamic_cast<genfit::KalmanFitterInfo*>(fitTrack->getPointWithMeasurement(i)->getFitterInfo(rep))->getWeights();
            std::cout << "hit " << i;
            if (outlierTrue[i]) std::cout << " is an OUTLIER";
            std::cout << " weights: ";
            for (unsigned int j=0; j<dafWeights.size(); ++j){
              std::cout << dafWeights[j] << "  ";
            }
            std::cout << "   l/r: " << leftRightTrue[i];
            std::cout << "\n";
          }
          int trueSide = leftRightTrue[i];
          if (trueSide == 0) continue; // not a wire measurement
          if (outlierTrue[i]) continue; // an outlier
          std::vector<double> dafWeightLR = dynamic_cast<genfit::KalmanFitterInfo*>(fitTrack->getPointWithMeasurement(i)->getFitterInfo(rep))->getWeights();
          if(dafWeightLR.size() != 2)
            continue;

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

        for (unsigned int i=0; i<outlierTrue.size(); ++i){
          if (! fitTrack->getPointWithMeasurement(i)->hasFitterInfo(rep))
            continue;

          std::vector<double> dafWeights = dynamic_cast<genfit::KalmanFitterInfo*>(fitTrack->getPointWithMeasurement(i)->getFitterInfo(rep))->getWeights();

          if (outlierTrue[i]) { // an outlier
            for (unsigned int j=0; j<dafWeights.size(); ++j){
              weights->Fill(dafWeights[j]-1);
            }
          }
          else if (leftRightTrue[i] == 0) { // only for non wire hits
            for (unsigned int j=0; j<dafWeights.size(); ++j){
              weights->Fill(dafWeights[j]);
            }
          }
        }

      }


      fitTrack->prune();
      if (debug) {
        std::cout<<" pruned track: ";
        fitTrack->Print();
      }

     /* if (debug) std::cerr<<"Fill Tree ..." << std::endl;
      tree->Fill();*/
#endif


  }// end loop over events

  CALLGRIND_STOP_INSTRUMENTATION;
  CALLGRIND_DUMP_STATS;

  std::cout<<"maxWeight = " << maxWeight << std::endl;
  std::cout<<"avg nr iterations = " << (double)nTotalIter/nSuccessfullFits << std::endl;
  std::cout<<"fit efficiency = " << (double)nSuccessfullFits/nEvents << std::endl;

  std::cout<<"avg nr iterations (2nd rep) = " << (double)nTotalIterSecond/nSuccessfullFitsSecond << std::endl;
  std::cout<<"fit efficiency (2nd rep) = " << (double)nSuccessfullFitsSecond/nEvents << std::endl;


#ifndef VALGRIND
  /*if (debug) std::cout<<"Write Tree ...";
  tree->Write();
  if (debug) std::cout<<"... done"<<std::endl;*/

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



  TCanvas* c3 = new TCanvas();
  //c3->Divide(2,3);

  c3->cd(1);
  trackLenRes->Fit("gaus");
  trackLenRes->Draw();

  c3->Write();

  if (debug) std::cout<<"... done"<<std::endl;

  // open event display
  display->setOptions("ABDEFHMPT"); // G show geometry
  if (matFX) display->setOptions("ABDEFGHMPT"); // G show geometry
  display->open();

  //rootapp->Run();

  //file->Close();
#endif

  //if (debug) std::cout<<"... closed file"<<std::endl;

}


