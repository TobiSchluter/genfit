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
#include <GFGbl.h>
#include <MeasuredStateOnPlane.h>
#include <MeasurementOnPlane.h>
#include <FullMeasurement.h>
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
#include <RKTrackRepEnergy.h>
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
#include <getopt.h>

void setup() {
  // gStyle->SetCanvasPreferGL(true);
  TGeoManager *geom = new TGeoManager("simple1", "Simple geometry");

  //--- define some materials
  TGeoMaterial *matVacuum = new TGeoMaterial("Vacuum", 0,0,0);
  TGeoMaterial *matAl = new TGeoMaterial("Al", 26.98,13,2.7);
  //   //--- define some media
  TGeoMedium *Vacuum = new TGeoMedium("Vacuum",1, matVacuum);
  TGeoMedium *Al = new TGeoMedium("Root Material",2, matAl);

  //--- define the transformations
  TGeoTranslation *tr1 = new TGeoTranslation(0., 0, 0.);

  //--- make the top container volume
  Double_t worldx = 100.;
  Double_t worldy = 100.;
  Double_t worldz = 100.;
  TGeoVolume *top = geom->MakeBox("TOP", Vacuum, worldx, worldy, worldz);
  geom->SetTopVolume(top);
  TGeoVolume *box = geom->MakeBox("BOX", Al, worldx, .1, worldz);
  top->AddNode(box, 1, tr1);
  geom->CloseGeometry();
   
  genfit::FieldManager::getInstance()->init(new genfit::ConstField(0.,0.,15));
  genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());
}

using namespace genfit;

int main()
{
  setup();

  RKTrackRep* rk1 = new RKTrackRep(211);
  RKTrackRepEnergy* rk2 = new RKTrackRepEnergy(211);

  StateOnPlane sop1(rk1);
  StateOnPlane sop2(rk2);

  TVector3 pos(0, -1, 0);
  TVector3 mom(0, 1, 0);
  sop1.setPosMom(pos, mom);
  sop2.setPosMom(pos, mom);

  SharedPlanePtr target(new DetPlane(TVector3(0, 1, 0), TVector3(0, 1, 0)));
  sop1.extrapolateToPlane(target);
  sop1.Print();

  //rk2->setDebugLvl(1);
  sop2.extrapolateToPlane(target);
  sop2.Print();

  return 0;
}
