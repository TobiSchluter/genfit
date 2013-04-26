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
#include <KalmanFitterInfo.h>
#include <MaterialInfo.h>
#include <MeasuredStateOnPlane.h>
#include <MeasurementOnPlane.h>
#include <PlanarMeasurement.h>
#include <PlanarPixelMeasurement.h>
#include <PlanarStripMeasurement.h>
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

bool compareForthBackExtrapolation() {

  double epsilonLen = 1.E-4; // 1 mu
  double epsilonMom = 1.E-5; // 10 keV

  int pdg = randomPdg();
  genfit::AbsTrackRep* rep;
  rep = new genfit::RKTrackRep(pdg);

  //TVector3 pos(0,0,0);
  TVector3 pos(0+gRandom->Gaus(0,0.1),0+gRandom->Gaus(0,0.1),0+gRandom->Gaus(0,0.1));
  TVector3 mom(0,0.5,0.);
  mom *= randomSign();


  genfit::StateOnPlane state(rep);
  rep->setPosMom(&state, pos, mom);

  genfit::SharedPlanePtr origPlane = state.getPlane();
  genfit::SharedPlanePtr plane(new genfit::DetPlane(TVector3(0,randomSign()*10,0), TVector3(0,randomSign()*1,0)));

  genfit::StateOnPlane origState(state);

  // forth
  double extrapLen(0);
  try {
    extrapLen = rep->extrapolateToPlane(&state, plane);
  }
  catch (genfit::Exception& e) {
    std::cerr << e.what();

    delete rep;
    return false;
  }

  // back
  double backExtrapLen(0);
  try {
    backExtrapLen = rep->extrapolateToPlane(&state, origPlane);
  }
  catch (genfit::Exception& e) {
    std::cerr << e.what();

    delete rep;
    return false;
  }

  // compare
  if ((rep->getPos(&origState) - rep->getPos(&state)).Mag() > epsilonLen ||
      (rep->getMom(&origState) - rep->getMom(&state)).Mag() > epsilonMom ||
      fabs(extrapLen + backExtrapLen) > epsilonLen) {

    origState.Print();
    state.Print();

    std::cerr << "pos difference = " << (rep->getPos(&origState) - rep->getPos(&state)).Mag() << "\n";
    std::cerr << "mom difference = " << (rep->getMom(&origState) - rep->getMom(&state)).Mag() << "\n";
    std::cerr << "len difference = " << extrapLen + backExtrapLen << "\n";

    std::cerr << std::endl;

    delete rep;
    return false;
  }

  delete rep;
  return true;

}

//=====================================================================================================================
//=====================================================================================================================
//=====================================================================================================================
//=====================================================================================================================


int main() {

  const double BField = 15.;       // kGauss
  const bool debug = true;

  gRandom->SetSeed(10);
  signal(SIGSEGV, handler);   // install our handler

  // init geometry and mag. field
  TGeoManager* geom = new TGeoManager("Geometry", "Geane geometry");
  TGeoManager::Import("genfitGeom.root");
  genfit::FieldManager::getInstance()->init(new genfit::ConstField(0.,0.,BField));
  genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());

  /*genfit::MaterialEffects::getInstance()->setEnergyLossBetheBloch(false);
  genfit::MaterialEffects::getInstance()->setNoiseBetheBloch(false);
  genfit::MaterialEffects::getInstance()->setNoiseCoulomb(false);
  genfit::MaterialEffects::getInstance()->setEnergyLossBrems(false);
  genfit::MaterialEffects::getInstance()->setNoiseBrems(false);*/

  /*genfit::Track* testTrack = new genfit::Track();


  TString outname = "test.root";
  TFile *file = TFile::Open(outname,"RECREATE");
  TTree *tree = new TTree("t","Tracks");
  tree->Branch("testTracks","genfit::Track",&testTrack);





  tree->Fill();
  if (debug) std::cerr<<"Write Tree ...";
  tree->Write();
  file->Close();*/

  for (unsigned int i=0; i<100; ++i) {
    if (!compareForthBackExtrapolation()) {
      std::cout << "failed in " << i << "\n";
      break;
    }
  }



  return 0;
}






