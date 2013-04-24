#include <iostream>
#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>

#include <AbsFinitePlane.h>
#include <AbsFitterInfo.h>
#include <AbsMeasurement.h>
#include <AbsTrackRep.h>
#include <DetPlane.h>
#include <Exception.h>
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
  signal(SIGSEGV, handler);   // install our handler
  std::cerr<<"main"<<std::endl;

  const bool debug = true;

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
  TVector3 mom(1,1,1);

  genfit::StateOnPlane* state = new genfit::StateOnPlane(rep);
  rep->setPosMom(state, pos, mom);

  rep->getPos(state).Print();
  rep->getMom(state).Print();




  return 0;
}
