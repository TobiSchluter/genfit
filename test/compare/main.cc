#include <iostream>

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

#include <TGeoManager.h>

void setup() {
  // gStyle->SetCanvasPreferGL(true);
  TGeoManager *geom = new TGeoManager("simple1", "Simple geometry");

  //--- define some materials
  TGeoMaterial *matVacuum = new TGeoMaterial("Vacuum", 0,0,0);
  TGeoMaterial *matAl = new TGeoMaterial("Al", 26.98,13,2.7);
  TGeoMaterial *matAlLowRho = new TGeoMaterial("Al", 26.98,13,0.27);
  //   //--- define some media
  TGeoMedium *Vacuum = new TGeoMedium("Vacuum",1, matVacuum);
  TGeoMedium *Al = new TGeoMedium("Al",2, matAl);
  TGeoMedium *AlLowRho = new TGeoMedium("AlLowRho",3, matAl);

  //--- define the transformations
  TGeoTranslation *tr1 = new TGeoTranslation(0., 0, 0.);
  TGeoTranslation *tr2 = new TGeoTranslation(0., .5, 0.);
  TGeoTranslation *tr3 = new TGeoTranslation(0., -1, 0.);
  TGeoTranslation *tr4 = new TGeoTranslation(0., -1.002, 0.);
  TGeoTranslation *tr5 = new TGeoTranslation(0., -1.004, 0.);

  //--- make the top container volume
  Double_t worldx = 100.;
  Double_t worldy = 100.;
  Double_t worldz = 100.;
  TGeoVolume *top = geom->MakeBox("TOP", Vacuum, worldx, worldy, worldz);
  geom->SetTopVolume(top);
  TGeoVolume *box = geom->MakeBox("AlBOX", Al, worldx, .01, worldz);
  top->AddNode(box, 1, tr1);
  top->AddNode(box, 2, tr2);
  TGeoVolume *smallBox = geom->MakeBox("smallBox", Al, worldx, .001, worldz);
  top->AddNode(smallBox, 1, tr3);
  TGeoVolume *smallBox2 = geom->MakeBox("smallBoxLow", AlLowRho, worldx, .001, worldz);
  top->AddNode(smallBox2, 1, tr4);
  top->AddNode(smallBox, 2, tr5);

  geom->CloseGeometry();
   
  genfit::FieldManager::getInstance()->init(new genfit::ConstField(0.,0.,15.));
  genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());
}

using namespace genfit;

int main()
{
  setup();

  RKTrackRep* rk1 = new RKTrackRep(-211);
  RKTrackRepEnergy* rk2 = new RKTrackRepEnergy(-211);

  MeasuredStateOnPlane sop1(rk1);
  MeasuredStateOnPlane sop2(rk2);

  TVector3 pos(0, -2, 0);
  TVector3 mom(0.001, .1, 0.01);
  TMatrixDSym cov(6);
  for (int i = 0; i < 6; ++i)
    cov(i, i) = 1.e-2;
  SharedPlanePtr start(new DetPlane(TVector3(0, -2, 0), TVector3(0, 1, 0)));
  sop1.setPosMomCov(pos, mom, cov);
  sop2.setPosMomCov(pos, mom, cov);

  sop1.Print();
  sop2.Print();

  MeasuredStateOnPlane sop3(sop1);
  MeasuredStateOnPlane sop4(sop2);

  SharedPlanePtr target(new DetPlane(TVector3(0, 1, 0), TVector3(0, 1, 0)));
  SharedPlanePtr middle(new DetPlane(TVector3(0, .2, 0), TVector3(0, 1, 0)));

  std::cout << "forth ->->->->->->->->->->->->->->->->->" << std::endl;

  sop1.extrapolateToPlane(target);
  sop1.Print();
  TMatrixD jac(5,5);
  TMatrixDSym noise(5);
  TVectorD d(5);
  rk1->getForwardJacobianAndNoise(jac, noise, d);
  jac.Print();
  noise.Print();

  //rk2->setDebugLvl(1);
  sop2.extrapolateToPlane(target);
  sop2.Print();
  rk2->getForwardJacobianAndNoise(jac, noise, d);
  jac.Print();
  noise.Print();

  std::cout << "back <-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-" << std::endl;

  sop1.extrapolateToPlane(start);
  sop1.Print();

  //rk2->setDebugLvl(1);
  sop2.extrapolateToPlane(start);
  sop2.Print();

  std::cout << "forth in two steps -->-->-->-->-->-->-->-->" << std::endl;

  sop3.extrapolateToPlane(middle);
  //sop3.Print();
  sop3.extrapolateToPlane(target);
  sop3.Print();

  sop4.extrapolateToPlane(middle);
  //sop4.Print();
  sop4.extrapolateToPlane(target);
  sop4.Print();

  std::cout << "back in two steps <--<--<--<--<--<--<--<--" << std::endl;

  sop3.extrapolateToPlane(middle);
  //sop3.Print();
  sop3.extrapolateToPlane(start);
  sop3.Print();

  sop4.extrapolateToPlane(middle);
  //sop4.Print();
  sop4.extrapolateToPlane(start);
  sop4.Print();

  return 0;
}
