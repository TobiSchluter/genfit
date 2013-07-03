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


bool compareMatrices(TMatrixTBase<double>& A, TMatrixTBase<double>& B, double maxAbsErr, double maxRelErr) {
  bool retVal = true;
  for (int i=0; i<A.GetNrows(); ++i) {
    for (int j=0; j<A.GetNcols(); ++j) {
      double absErr = A(i,j) - B(i,j);
      if ( fabs(absErr) > maxAbsErr ) {
        double relErr = A(i,j)/B(i,j) - 1;
        if ( fabs(relErr) > maxRelErr ) {
          std::cout << "compareMatrices: A("<<i<<","<<j<<") = " << A(i,j) << "  B("<<i<<","<<j<<") = " << B(i,j) << "     absErr = " << absErr << "    relErr = " << relErr << "\n";
          retVal = false;
        }
      }
    }
  }
  return retVal;
}

bool isCovMatrix(TMatrixTBase<double>& cov) {

  if (!(cov.IsSymmetric())) {
    std::cout << "isCovMatrix: not symmetric\n";
    return false;
  }

  for (int i=0; i<cov.GetNrows(); ++i) {
    for (int j=0; j<cov.GetNcols(); ++j) {
       if (isnan(cov(i,j))) {
         std::cout << "isCovMatrix: element isnan\n";
         return false;
       }
       if (i==j && cov(i,j) < 0) {
         std::cout << "isCovMatrix: negative diagonal element\n";
         return false;
       }
    }
  }

  return true;
}



bool checkSetGetPosMom() {

  double epsilonLen = 1.E-10;
  double epsilonMom = 1.E-10;

  int pdg = randomPdg();
  genfit::AbsTrackRep* rep;
  rep = new genfit::RKTrackRep(pdg);

  //TVector3 pos(0,0,0);
  TVector3 pos(gRandom->Gaus(0,0.1),gRandom->Gaus(0,0.1),gRandom->Gaus(0,0.1));
  TVector3 mom(0,0.5,gRandom->Gaus(0,0.3));
  mom.SetMag(0.5);
  mom *= randomSign();


  genfit::StateOnPlane state(rep);
  rep->setPosMom(&state, pos, mom);

  // check if we can set another position in the same plane
  if (randomSign() == 1) {
    genfit::SharedPlanePtr plane = state.getPlane();
    const TVector3& u = plane->getU();
    const TVector3& v = plane->getV();

    // random position on plane
    pos += gRandom->Gaus() * u;
    pos += gRandom->Gaus() * v;

    // new random momentum
    mom.SetXYZ(0,0.5,gRandom->Gaus(0,0.3));
    mom.SetMag(0.5);
    mom *= randomSign();

    rep->setPosMom(&state, pos, mom);

    // check if plane has changed
    if (state.getPlane() != plane) {
      std::cout << "plane has changed unexpectedly! \n";
      delete rep;
      return false;
    }
  }


  // compare
  if ((pos - rep->getPos(&state)).Mag() > epsilonLen ||
      (mom - rep->getMom(&state)).Mag() > epsilonMom) {

    state.Print();

    std::cout << "pos difference = " << (pos - rep->getPos(&state)).Mag() << "\n";
    std::cout << "mom difference = " << (mom - rep->getMom(&state)).Mag() << "\n";

    std::cout << std::endl;

    delete rep;
    return false;
  }

  delete rep;
  return true;

}


bool compareForthBackExtrapolation() {

  double epsilonLen = 5.E-5; // 0.5 mu
  double epsilonMom = 1.E-4; // 100 keV

  int pdg = randomPdg();
  genfit::AbsTrackRep* rep;
  rep = new genfit::RKTrackRep(pdg);

  //TVector3 pos(0,0,0);
  TVector3 pos(gRandom->Gaus(0,0.1),gRandom->Gaus(0,0.1),gRandom->Gaus(0,0.1));
  TVector3 mom(0,0.5,gRandom->Gaus(0,0.3));
  mom.SetMag(0.5);
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

    std::cout << "pos difference = " << (rep->getPos(&origState) - rep->getPos(&state)).Mag() << "\n";
    std::cout << "mom difference = " << (rep->getMom(&origState) - rep->getMom(&state)).Mag() << "\n";
    std::cout << "len difference = " << extrapLen + backExtrapLen << "\n";

    std::cout << std::endl;

    delete rep;
    return false;
  }

  delete rep;
  return true;

}


bool compareForthBackJacNoise() {

  double epsilonJac = 0.004; // absolute  // best reached: 5.E-5
  double deltaJac = 0.11; // relative     // best reached: 0.01
  double epsilonNoise = 2.E-3;
  double deltaNoise = 0.02;

  int pdg = randomPdg();
  genfit::AbsTrackRep* rep;
  rep = new genfit::RKTrackRep(pdg);

  //TVector3 pos(0,0,0);
  //TVector3 mom(0,1,2);
  TVector3 pos(gRandom->Gaus(0,0.1),gRandom->Gaus(0,0.1),gRandom->Gaus(0,0.1));
  TVector3 mom(0, 0.5, gRandom->Gaus(0, 1));
  mom *= randomSign();
  mom.SetMag(gRandom->Uniform(2)+0.3);
  //mom.SetMag(3);

  TMatrixD jac_f, jac_fi, jac_b, jac_bi;
  TMatrixDSym noise_f, noise_fi, noise_b, noise_bi;


  // original state and plane
  genfit::MeasuredStateOnPlane state(rep);
  rep->setPosMom(&state, pos, mom);

  static const double smear = 0.2;
  TVector3 normal(mom);
  normal.SetMag(1);
  normal.SetXYZ(gRandom->Gaus(normal.X(), smear),
      gRandom->Gaus(normal.Y(), smear),
      gRandom->Gaus(normal.Z(), smear));
  genfit::DetPlane* origPlanePtr = new genfit::DetPlane (pos, normal);
  //genfit::DetPlane* origPlanePtr = new genfit::DetPlane (pos, TVector3(1,0,0), TVector3(0,0,1));
  double rotAngleOrig = gRandom->Uniform(2.*TMath::Pi());
  origPlanePtr->rotate(rotAngleOrig);
  genfit::SharedPlanePtr origPlane(origPlanePtr);
  //genfit::SharedPlanePtr origPlane = state.getPlane();
  rep->extrapolateToPlane(&state, origPlane);

  const genfit::StateOnPlane origState(state);


  // dest plane
  normal = mom;
  normal.SetMag(1);
  normal.SetXYZ(gRandom->Gaus(normal.X(), smear),
      gRandom->Gaus(normal.Y(), smear),
      gRandom->Gaus(normal.Z(), smear));
  TVector3 dest(mom);
  dest.SetMag(10);
  genfit::DetPlane* planePtr = new genfit::DetPlane (dest, normal);
  //genfit::DetPlane* planePtr = new genfit::DetPlane (dest, TVector3(1,0,0), TVector3(0,0,1));
  double rotAngle = gRandom->Uniform(2.*TMath::Pi());
  planePtr->rotate(rotAngle);
  genfit::SharedPlanePtr plane(planePtr);

 /* genfit::DetPlane* planePtr = new genfit::DetPlane (*origPlane);
  planePtr->setO(TVector3(0,randomSign()*10,0));
  //planePtr->rotate(rotAngle);
  genfit::SharedPlanePtr plane(planePtr);
*/

  // manual calculation
  TMatrixD jac_f_man;
  rep->calcJacobianNumerically(&origState, plane, jac_f_man);

  double extrapLen(0);

  // forth
  try {
    //std::cout << "DO FORTH EXTRAPOLATION \n";
    extrapLen = rep->extrapolateToPlane(&state, plane);
    //std::cout << "GET INFO FOR FORTH EXTRAPOLATION \n";
    rep->getForwardJacobianAndNoise(jac_f, noise_f);
    rep->getBackwardJacobianAndNoise(jac_fi, noise_fi);
  }
  catch (genfit::Exception& e) {
    std::cerr << e.what();

    delete rep;
    return false;
  }

  // back
  try {
    //std::cout << "DO BACK EXTRAPOLATION \n";
    rep->extrapolateToPlane(&state, origPlane);
    //std::cout << "GET INFO FOR BACK EXTRAPOLATION \n";
    rep->getForwardJacobianAndNoise(jac_b, noise_b);
    rep->getBackwardJacobianAndNoise(jac_bi, noise_bi);
  }
  catch (genfit::Exception& e) {
    std::cerr << e.what();

    delete rep;
    return false;
  }


  // compare
  if (!isCovMatrix(state.getCov()) ||
      !compareMatrices(jac_f, jac_f_man, 2*epsilonJac, 2*deltaJac) ||
      !compareMatrices(jac_f, jac_bi, epsilonJac, deltaJac) ||
      !compareMatrices(jac_b, jac_fi, epsilonJac, deltaJac) ||
      !compareMatrices(noise_f, noise_bi, epsilonNoise, deltaNoise) ||
      !compareMatrices(noise_b, noise_fi, epsilonNoise, deltaNoise)) {

    double rotAngleDeg = rotAngle / TMath::Pi() * 180;
    if (rotAngleDeg > 180) rotAngleDeg -= 360;
    std::cout << "rotAngle = " << rotAngleDeg << " ° \n";
    std::cout << "extrapLen = " << extrapLen << " cm \n";

    origState.Print();
    state.Print();

    std::cout << "origPlane "; origPlane->Print();
    std::cout << "plane "; plane->Print();

    std::cout << "jac_f_man = "; jac_f_man.Print();
    std::cout << "jac_f = "; jac_f.Print();
    std::cout << "jac_bi = "; jac_bi.Print();
    std::cout << "jac_b = "; jac_b.Print();
    std::cout << "jac_fi = "; jac_fi.Print();

    std::cout << "noise_f = "; noise_f.Print();
    std::cout << "noise_bi = "; noise_bi.Print();
    std::cout << "noise_b = "; noise_b.Print();
    std::cout << "noise_fi = "; noise_fi.Print();

    /*std::cout << "jac difference (jac_f - jac_bi) = "; (jac_f - jac_bi).Print();
    std::cout << "jac difference (jac_b - jac_fi) = "; (jac_b - jac_fi).Print();
    std::cout << "noise difference (noise_f - noise_bi) = "; (noise_f - noise_bi).Print();
    std::cout << "noise difference (noise_b - noise_fi) = "; (noise_b - noise_fi).Print();*/

    std::cout << std::endl;

    delete rep;
    return false;
  }

  delete rep;
  return true;
}


bool checkStopAtBoundary() {

  double epsilonLen = 1.E-4; // 1 mu
  //double epsilonMom = 1.E-5; // 10 keV

  double matRadius(1.);

  int pdg = randomPdg();
  genfit::AbsTrackRep* rep;
  rep = new genfit::RKTrackRep(pdg);

  //TVector3 pos(0,0,0);
  TVector3 pos(gRandom->Gaus(0,0.1),gRandom->Gaus(0,0.1),gRandom->Gaus(0,0.1));
  TVector3 mom(0,0.5,gRandom->Gaus(0,0.1));
  mom *= randomSign();


  genfit::StateOnPlane state(rep);
  rep->setPosMom(&state, pos, mom);

  genfit::SharedPlanePtr origPlane = state.getPlane();
  genfit::SharedPlanePtr plane(new genfit::DetPlane(TVector3(0,randomSign()*10,0), TVector3(0,randomSign()*1,0)));

  genfit::StateOnPlane origState(state);

  // forth
  try {
    rep->extrapolateToPlane(&state, plane, true);
  }
  catch (genfit::Exception& e) {
    std::cerr << e.what();

    delete rep;
    return false;
  }


  // compare
  if (fabs(rep->getPos(&state).Perp() - matRadius) > epsilonLen) {

      origState.Print();
      state.Print();

      std::cerr << "radius difference = " << rep->getPos(&state).Perp() - matRadius << "\n";

      std::cerr << std::endl;

      delete rep;
      return false;
    }

    delete rep;
    return true;

}


bool checkErrorPropagation() {

  //double epsilonLen = 1.E-4; // 1 mu
  //double epsilonMom = 1.E-5; // 10 keV

  int pdg = randomPdg();
  genfit::AbsTrackRep* rep;
  rep = new genfit::RKTrackRep(pdg);

  //TVector3 pos(0,0,0);
  TVector3 pos(gRandom->Gaus(0,0.1),gRandom->Gaus(0,0.1),gRandom->Gaus(0,0.1));
  TVector3 mom(0,0.5,gRandom->Gaus(0,0.1));
  mom *= randomSign();


  genfit::MeasuredStateOnPlane state(rep);
  rep->setPosMom(&state, pos, mom);

  genfit::SharedPlanePtr origPlane = state.getPlane();
  genfit::SharedPlanePtr plane(new genfit::DetPlane(TVector3(0,randomSign()*50,0), TVector3(0,randomSign()*1,0)));

  genfit::MeasuredStateOnPlane origState(state);

  // forth
  try {
    rep->extrapolateToPlane(&state, plane);
  }
  catch (genfit::Exception& e) {
    std::cerr << e.what();

    delete rep;
    return false;
  }


  // check
  if (!isCovMatrix(state.getCov())) {

    origState.Print();
    state.Print();

    delete rep;
    return false;
  }

  delete rep;
  return true;

}


bool checkExtrapolateToLine() {

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
  genfit::StateOnPlane origState(state);

  TVector3 linePoint(gRandom->Gaus(),randomSign()*10+gRandom->Gaus(),gRandom->Gaus());
  TVector3 lineDirection(gRandom->Gaus(),gRandom->Gaus(),randomSign()*10+gRandom->Gaus());

  // forth
  try {
    rep->extrapolateToLine(&state, linePoint, lineDirection, false);
  }
  catch (genfit::Exception& e) {
    std::cerr << e.what();

    delete rep;
    return false;
  }


  // compare
  if (fabs(state.getPlane()->distance(linePoint)) > epsilonLen ||
      fabs(state.getPlane()->distance(linePoint+lineDirection)) > epsilonLen ||
      (rep->getMom(&state).Unit() * state.getPlane()->getNormal()) > epsilonMom) {

      origState.Print();
      state.Print();

      std::cout << "distance of linePoint to plane = " << state.getPlane()->distance(linePoint) << "\n";
      std::cout << "distance of linePoint+lineDirection to plane = " << state.getPlane()->distance(linePoint+lineDirection) << "\n";
      std::cout << "direction * plane normal = " << rep->getMom(&state).Unit() * state.getPlane()->getNormal() << "\n";

      delete rep;
      return false;
    }

    delete rep;
    return true;

}


bool checkExtrapolateToPoint() {

  double epsilonLen = 1.E-4; // 1 mu
  double epsilonMom = 1.E-5; // 10 keV

  int pdg = randomPdg();
  genfit::AbsTrackRep* rep;
  rep = new genfit::RKTrackRep(pdg);

  //TVector3 pos(0,0,0);
  TVector3 pos(gRandom->Gaus(0,0.1),gRandom->Gaus(0,0.1),gRandom->Gaus(0,0.1));
  TVector3 mom(0,0.5,gRandom->Gaus(0,0.1));
  mom *= randomSign();


  genfit::StateOnPlane state(rep);
  rep->setPosMom(&state, pos, mom);

  genfit::SharedPlanePtr origPlane = state.getPlane();
  genfit::StateOnPlane origState(state);

  TVector3 point(gRandom->Gaus(),randomSign()*10+gRandom->Gaus(),gRandom->Gaus());

  // forth
  try {
    rep->extrapolateToPoint(&state, point, false);
  }
  catch (genfit::Exception& e) {
    std::cerr << e.what();

    delete rep;
    return false;
  }


  // compare
  if (fabs(state.getPlane()->distance(point)) > epsilonLen ||
      fabs((rep->getMom(&state).Unit() * state.getPlane()->getNormal())) - 1 > epsilonMom) {

      origState.Print();
      state.Print();

      std::cout << "distance of point to plane = " << state.getPlane()->distance(point) << "\n";
      std::cout << "direction * plane normal = " << rep->getMom(&state).Unit() * state.getPlane()->getNormal() << "\n";

      delete rep;
      return false;
    }

    delete rep;
    return true;

}


bool checkExtrapolateToCylinder() {

  double epsilonLen = 1.E-4; // 1 mu
  //double epsilonMom = 1.E-5; // 10 keV

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
  genfit::StateOnPlane origState(state);

  const TVector3 linePoint(gRandom->Gaus(0,5), gRandom->Gaus(0,5), gRandom->Gaus(0,5));
  const TVector3 lineDirection(gRandom->Gaus(),gRandom->Gaus(),2+gRandom->Gaus());
  const double radius = gRandom->Uniform(10);

  // forth
  try {
    rep->extrapolateToCylinder(&state, radius, linePoint, lineDirection, false);
  }
  catch (genfit::Exception& e) {
    std::cerr << e.what();

    delete rep;

    static const char* bla = "cannot solve";
    const char* what = e.what();
    if (strstr(what, bla))
      return true;
    return false;
  }

  TVector3 pocaOnLine(lineDirection);
  double t = 1./(pocaOnLine.Mag2()) * ((rep->getPos(&state)*pocaOnLine) - (linePoint*pocaOnLine));
  pocaOnLine *= t;
  pocaOnLine += linePoint;

  TVector3 radiusVec = rep->getPos(&state) - pocaOnLine;

  // compare
  if (fabs(lineDirection*radiusVec) > epsilonLen ||
      fabs(radiusVec.Mag()-radius) > epsilonLen) {

      origState.Print();
      state.Print();

      std::cout << "lineDirection*radiusVec = " << lineDirection*radiusVec << "\n";
      std::cout << "radiusVec.Mag()-radius = " << radiusVec.Mag()-radius << "\n";

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
  //const bool debug = true;

  gRandom->SetSeed(10);
  signal(SIGSEGV, handler);   // install our handler

  // init geometry and mag. field
  TGeoManager* geom = new TGeoManager("Geometry", "Geane geometry");
  TGeoManager::Import("genfitGeom.root");
  genfit::FieldManager::getInstance()->init(new genfit::ConstField(0.,0.,BField));
  genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());

  TDatabasePDG::Instance()->GetParticle(211);

  /*genfit::MaterialEffects::getInstance()->setEnergyLossBetheBloch(false);
  genfit::MaterialEffects::getInstance()->setNoiseBetheBloch(false);
  genfit::MaterialEffects::getInstance()->setNoiseCoulomb(false);
  genfit::MaterialEffects::getInstance()->setEnergyLossBrems(false);
  genfit::MaterialEffects::getInstance()->setNoiseBrems(false);*/
  //genfit::MaterialEffects::getInstance()->setNoEffects();
  /*genfit::Track* testTrack = new genfit::Track();


  TString outname = "test.root";
  TFile *file = TFile::Open(outname,"RECREATE");
  TTree *tree = new TTree("t","Tracks");
  tree->Branch("testTracks","genfit::Track",&testTrack);





  tree->Fill();
  if (debug) std::cout<<"Write Tree ...";
  tree->Write();
  file->Close();*/


  unsigned int nFailed(0);
  const unsigned int nTests(100);

  for (unsigned int i=0; i<nTests; ++i) {

    if (!checkSetGetPosMom()) {
      std::cout << "failed checkSetGetPosMom nr" << i << "\n";
      ++nFailed;
    }

    if (!compareForthBackExtrapolation()) {
      std::cout << "failed compareForthBackExtrapolation nr" << i << "\n";
      ++nFailed;
    }

    if (!checkStopAtBoundary()) {
      std::cout << "failed checkStopAtBoundary nr" << i << "\n";
      ++nFailed;
    }

    if (!checkErrorPropagation()) {
      std::cout << "failed checkErrorPropagation nr" << i << "\n";
      ++nFailed;
    }

    if (!checkExtrapolateToLine()) {
      std::cout << "failed checkExtrapolateToLine nr" << i << "\n";
      ++nFailed;
    }

    if (!checkExtrapolateToPoint()) {
      std::cout << "failed checkExtrapolateToPoint nr" << i << "\n";
      ++nFailed;
    }

    if (!checkExtrapolateToCylinder()) {
      std::cout << "failed checkExtrapolateToCylinder nr" << i << "\n";
      ++nFailed;
    }

    if (!compareForthBackJacNoise()) {
      std::cout << "failed compareForthBackJacNoise nr" << i << "\n";
      ++nFailed;
    }

  }

  std::cout << "failed " << nFailed << " of " << nTests << " Tests." << std::endl;
  if (nFailed == 0) {
    std::cout << "passed all tests!" << std::endl;
  }




  return 0;
}


