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
#include <RKTrackRepTime.h>
#include <StepLimits.h>
#include <TGeoMaterialInterface.h>

#include <TGeoManager.h>

void setup() {
  // gStyle->SetCanvasPreferGL(true);
  TGeoManager *geom = new TGeoManager("simple1", "Simple geometry");

  //--- define some materials
  TGeoMaterial *matVacuum = new TGeoMaterial("Vacuum", 0,0,0);
  TGeoMaterial *matAl = new TGeoMaterial("Al", 26.98,13,2.7);
  TGeoMaterial *matAir = new TGeoMaterial("Air", 14.61, 7.3, 0.1205e-2);
  //   //--- define some media
  TGeoMedium *Vacuum = new TGeoMedium("Vacuum",1, matVacuum);
  TGeoMedium *Al = new TGeoMedium("Al",2, matAl);
  TGeoMedium *Air = new TGeoMedium("Air", 1, matAir);
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
  TGeoVolume *top = geom->MakeBox("TOP", Air, worldx, worldy, worldz);
  geom->SetTopVolume(top);
  TGeoVolume *box = geom->MakeBox("AlBOX", Al, worldx, .01, worldz);
  top->AddNode(box, 1, tr1);
  top->AddNode(box, 2, tr2);
  TGeoVolume *smallBox = geom->MakeBox("smallBox", Al, worldx, .001, worldz);
  top->AddNode(smallBox, 1, tr3);
  TGeoVolume *smallBox2 = geom->MakeBox("smallBoxLow", Air, worldx, .001, worldz);
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

  AbsTrackRep* rk1 = new RKTrackRepTime(-211);
  //AbsTrackRep* rk1 = new RKTrackRepEnergy(-211);

  MeasuredStateOnPlane mop1(rk1);

  TVector3 pos(0, -2, 0);
  TVector3 mom(0, .1, .0);
  TMatrixDSym cov(6);
  for (int i = 0; i < 6; ++i)
    cov(i, i) = 1.e-2;
  SharedPlanePtr start(new DetPlane(TVector3(0, -2, 0), TVector3(0, 1, 0)));
  mop1.setPosMomCov(pos, mom, cov);

  MeasuredStateOnPlane mop3(mop1);

  SharedPlanePtr target(new DetPlane(TVector3(0, 1, 0), TVector3(0, 1, 0)));
  SharedPlanePtr middle(new DetPlane(TVector3(0, .02, 0), TVector3(0, 1, 0)));

  std::cout << "forth ->->->->->->->->->->->->->->->->->" << std::endl;

  mop1.extrapolateToPlane(target);
  mop1.Print();
  TMatrixD jac(5,5);
  TMatrixDSym noise(5);
  TVectorD d(5);
  std::cout << "analytically RKTrackRep" << std::endl;
  rk1->getForwardJacobianAndNoise(jac, noise, d);
  jac.Print();
  jac.Zero();
  std::cout << "numerically RKTrackRep:" << std::endl;
  rk1->calcJacobianNumerically(mop3, target, jac);
  jac.Print();

  return 0;
#if 0
  std::cout << "next" << std::endl;

  //rk2->setDebugLvl(1);
  mop2.extrapolateToPlane(target);
  mop2.Print();
  std::cout << "analytically RKTrackRepEnergy" << std::endl;
  TMatrixD jacf(5,5);
  rk2->getForwardJacobianAndNoise(jacf, noise, d);
  jacf.Print();
  TMatrixD jacfn(5,5);
  jac.Zero();
  std::cout << "numerically RKTrackRepEnergy:" << std::endl;
  rk2->calcJacobianNumerically(mop4, target, jacfn);
  std::cout << " jac " << std::endl;
  jacfn.Print();

  //return 0;

  std::cout << "back <-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-" << std::endl;

  MeasuredStateOnPlane mop1b(mop1);
  mop1b.extrapolateToPlane(start);
  mop1b.Print();

  std::cout << "analytically RKTrackRep" << std::endl;
  rk1->getForwardJacobianAndNoise(jac, noise, d);
  jac.Print();
  jac.Zero();
  std::cout << "numerically RKTrackRep:" << std::endl;
  rk1->calcJacobianNumerically(mop1, start, jac);
  jac.Print();

  MeasuredStateOnPlane mop2b(mop2);
  mop2b.extrapolateToPlane(start);
  mop2b.Print();
  TMatrixD jacb(5,5);
  std::cout << "analytically RKTrackRepE" << std::endl;
  rk2->getForwardJacobianAndNoise(jacb, noise, d);
  jacb.Print();
  TMatrixD jacbn(5,5);
  std::cout << "numerically RKTrackRepE:" << std::endl;
  rk2->calcJacobianNumerically(mop2, start, jacbn);
  jacbn.Print();

  std::cout << "product Jfw Jbw analytically" << std::endl;
  TMatrixD(jacb, TMatrixD::kMult, jacf).Print();
  std::cout << "product Jfw Jbw numerically" << std::endl;
  TMatrixD(jacbn, TMatrixD::kMult, jacfn).Print();


  return(0);
  std::cout << "forth in two steps -->-->-->-->-->-->-->-->" << std::endl;

  mop3.extrapolateToPlane(middle);
  //mop3.Print();
  mop3.extrapolateToPlane(target);
  mop3.Print();

  mop4.extrapolateToPlane(middle);
  //mop4.Print();
  mop4.extrapolateToPlane(target);
  mop4.Print();

  std::cout << "back in two steps <--<--<--<--<--<--<--<--" << std::endl;

  mop3.extrapolateToPlane(middle);
  //mop3.Print();
  mop3.extrapolateToPlane(start);
  mop3.Print();

  mop4.extrapolateToPlane(middle);
  //mop4.Print();
  mop4.extrapolateToPlane(start);
  mop4.Print();
#endif
  return 0;
}
