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


int main() {

  const double BField = 15.;       // kGauss
  const bool debug = true;


  signal(SIGSEGV, handler);   // install our handler

  // init geometry and mag. field
  TGeoManager* geom = new TGeoManager("Geometry", "Geane geometry");
  TGeoManager::Import("genfitGeom.root");
  genfit::FieldManager::getInstance()->init(new genfit::ConstField(0.,0.,BField));
  genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());

  /*genfit::Track* testTrack = new genfit::Track();


  TString outname = "test.root";
  TFile *file = TFile::Open(outname,"RECREATE");
  TTree *tree = new TTree("t","Tracks");
  tree->Branch("testTracks","genfit::Track",&testTrack);





  tree->Fill();
  if (debug) std::cerr<<"Write Tree ...";
  tree->Write();
  file->Close();*/


  genfit::AbsTrackRep* rep;
  rep = new genfit::RKTrackRep(211);

  TVector3 pos(0,0,0);
  TVector3 mom(0,0.5,0.);

  genfit::StateOnPlane* state = new genfit::StateOnPlane(rep);
  rep->setPosMom(state, pos, mom);
  state->Print();

  genfit::SharedPlanePtr origPlane = state->getPlane();

  genfit::SharedPlanePtr plane(new genfit::DetPlane(TVector3(0,10,0), TVector3(0,-1,0)));

  double extrapLen(0);
  try {
    extrapLen = rep->extrapolateToPlane(state, plane);
  }
  catch (genfit::Exception& e) {
    std::cerr << e.what();
  }

  std::cout << "extrapLen = " << extrapLen << "\n";

  state->Print();



  try {
    extrapLen = rep->extrapolateToPlane(state, origPlane);
  }
  catch (genfit::Exception& e) {
    std::cerr << e.what();
  }

  state->Print();



  return 0;
}
