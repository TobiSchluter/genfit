/* Copyright 2008-2015, Technische Universitaet Muenchen, Ludwig-Maximilians-Universität München
   Authors: Tobias Schlüter

   This file is part of GENFIT.

   GENFIT is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published
   by the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   GENFIT is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with GENFIT.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "RKTrackRepTime.h"
namespace genfit {
typedef struct TRKStep<8> RKStep;
typedef struct TExtrapStep<8> ExtrapStep;
}

#include <Exception.h>
#include <FieldManager.h>
#include <MaterialEffects.h>
#include <MeasuredStateOnPlane.h>
#include <MeasurementOnPlane.h>

#include <TDecompLU.h>
#include <TMath.h>
#include <TDatabasePDG.h>
#include <TGeoManager.h>

#include <iomanip>
#include <algorithm>

#define MINSTEP 0.001   // minimum step [cm] for Runge Kutta and iteration to POCA

namespace genfit {


RKTrackRepTime::RKTrackRepTime() :
  AbsTrackRep(),
  lastStartState_(this),
  lastEndState_(this),
  RKStepsFXStart_(0),
  RKStepsFXStop_(0),
  fJacobian_(6,6),
  fNoise_(6),
  useCache_(false),
  cachePos_(0)
{
  initArrays();
}


RKTrackRepTime::RKTrackRepTime(int pdgCode, char propDir) :
  AbsTrackRep(pdgCode, propDir),
  lastStartState_(this),
  lastEndState_(this),
  RKStepsFXStart_(0),
  RKStepsFXStop_(0),
  fJacobian_(6,6),
  fNoise_(6),
  useCache_(false),
  cachePos_(0)
{
  initArrays();
}


RKTrackRepTime::~RKTrackRepTime() {
  ;
}


double RKTrackRepTime::extrapolateToPlane(StateOnPlane& state,
    const SharedPlanePtr& plane,
    bool stopAtBoundary,
    bool calcJacobianNoise) const {

  if (debugLvl_ > 0) {
    std::cout << "RKTrackRepTime::extrapolateToPlane()\n";
  }


  if (state.getPlane() == plane) {
    if (debugLvl_ > 0) {
      std::cout << "state is already defined at plane. Do nothing! \n";
    }
    return 0;
  }

  checkCache(state, &plane);

  // to 8D
  M1x8 stateGlobal = {{0, 0, 0, 0, 0, 0, 0, 0}};
  getStateGlobal(state, stateGlobal);

  TMatrixDSym* covPtr(NULL);
  bool fillExtrapSteps(false);
  if (dynamic_cast<MeasuredStateOnPlane*>(&state) != NULL) {
    covPtr = &(static_cast<MeasuredStateOnPlane*>(&state)->getCov());
    fillExtrapSteps = true;
  }
  else if (calcJacobianNoise)
    fillExtrapSteps = true;

  // actual extrapolation
  bool isAtBoundary(false);
  double flightTime( 0. );
  double coveredDistance( Extrap(*(state.getPlane()), *plane, getCharge(state), getMass(state), isAtBoundary, stateGlobal, flightTime, fillExtrapSteps, covPtr, false, stopAtBoundary) );

  if (stopAtBoundary && isAtBoundary) {
    state.setPlane(SharedPlanePtr(new DetPlane(TVector3(stateGlobal[0], stateGlobal[1], stateGlobal[2]),
                                               TVector3(stateGlobal[3], stateGlobal[4], stateGlobal[5]))));
  }
  else {
    state.setPlane(plane);
  }

  // back to 6D
  getStateLocal(state, stateGlobal);
  setTime(state, getTime(state) + flightTime);
  //std::cout << getTime(state) << " " << stateGlobal[7] << std::endl;
  lastEndState_ = state;

  if (debugLvl_ > 0) {
    std::cout << "RKTrackRepTime::extrapolateToPlane(): coveredDistance = " << coveredDistance << std::endl;
  }

  return coveredDistance;
}

double RKTrackRepTime::extrapolateToLine(StateOnPlane& state,
    const TVector3& linePoint,
    const TVector3& lineDirection,
    bool stopAtBoundary,
    bool calcJacobianNoise) const {

  if (debugLvl_ > 0) {
    std::cout << "RKTrackRepTime::extrapolateToLine()\n";
  }

  checkCache(state, NULL);

  static const unsigned int maxIt(1000);

  // to 7D
  M1x8 stateGlobal;
  getStateGlobal(state, stateGlobal);

  bool fillExtrapSteps(false);
  if (dynamic_cast<MeasuredStateOnPlane*>(&state) != NULL) {
    fillExtrapSteps = true;
  }
  else if (calcJacobianNoise)
    fillExtrapSteps = true;

  double step(0.), lastStep(0.), maxStep(1.E99), angle(0), distToPoca(0), tracklength(0);
  double charge = getCharge(state);
  double mass = getMass(state);
  double flightTime = 0;
  TVector3 dir(stateGlobal[3], stateGlobal[4], stateGlobal[5]);
  TVector3 lastDir(0,0,0);
  TVector3 poca, poca_onwire;
  bool isAtBoundary(false);

  DetPlane startPlane(*(state.getPlane()));
  SharedPlanePtr plane(new DetPlane(linePoint, dir.Cross(lineDirection), lineDirection));
  unsigned int iterations(0);

  while(true){
    if(++iterations == maxIt) {
      Exception exc("RKTrackRepTime::extrapolateToLine ==> extrapolation to line failed, maximum number of iterations reached",__LINE__,__FILE__);
      exc.setFatal();
      throw exc;
    }

    lastStep = step;
    lastDir = dir;

    step = this->Extrap(startPlane, *plane, charge, mass, isAtBoundary, stateGlobal, flightTime, false, NULL, true, stopAtBoundary, maxStep);
    tracklength += step;

    dir.SetXYZ(stateGlobal[3], stateGlobal[4], stateGlobal[5]);
    poca.SetXYZ(stateGlobal[0], stateGlobal[1], stateGlobal[2]);
    poca_onwire = pocaOnLine(linePoint, lineDirection, poca);

    // check break conditions
    if (stopAtBoundary && isAtBoundary) {
      plane->setON(dir, poca);
      break;
    }

    angle = fabs(dir.Angle((poca_onwire-poca))-TMath::PiOver2()); // angle between direction and connection to point - 90 deg
    distToPoca = (poca_onwire-poca).Mag();
    if (angle*distToPoca < 0.1*MINSTEP) break;

    // if lastStep and step have opposite sign, the real normal vector lies somewhere between the last two normal vectors (i.e. the directions)
    // -> try mean value of the two (normalization not needed)
    if (lastStep*step < 0){
      dir += lastDir;
      maxStep = 0.5*fabs(lastStep); // make it converge!
    }

    startPlane = *plane;
    plane->setU(dir.Cross(lineDirection));
  }

  if (fillExtrapSteps) { // now do the full extrapolation with covariance matrix
    // make use of the cache
    lastEndState_.setPlane(plane);
    getStateLocal(lastEndState_, stateGlobal);

    tracklength = extrapolateToPlane(state, plane, false, true);
    lastEndState_.getAuxInfo()(1) = state.getAuxInfo()(1); // Flight time
  }
  else {
    state.setPlane(plane);
    getStateLocal(state, stateGlobal);
    state.getAuxInfo()(1) += flightTime;
  }

  if (debugLvl_ > 0) {
    std::cout << "RKTrackRepTime::extrapolateToLine(): Reached POCA after " << iterations+1 << " iterations. Distance: " << (poca_onwire-poca).Mag() << " cm. Angle deviation: " << dir.Angle((poca_onwire-poca))-TMath::PiOver2() << " rad \n";
  }

  lastEndState_ = state;

  return tracklength;
}


double RKTrackRepTime::extrapToPoint(StateOnPlane& state,
    const TVector3& point,
    const TMatrixDSym* G,
    bool stopAtBoundary,
    bool calcJacobianNoise) const {

  if (debugLvl_ > 0) {
    std::cout << "RKTrackRepTime::extrapolateToPoint()\n";
  }

  checkCache(state, NULL);

  static const unsigned int maxIt(1000);

  // to 7D
  M1x8 stateGlobal;
  getStateGlobal(state, stateGlobal);

  bool fillExtrapSteps(false);
  if (dynamic_cast<MeasuredStateOnPlane*>(&state) != NULL) {
    fillExtrapSteps = true;
  }
  else if (calcJacobianNoise)
    fillExtrapSteps = true;

  double step(0.), lastStep(0.), maxStep(1.E99), angle(0), distToPoca(0), tracklength(0);
  TVector3 dir(stateGlobal[3], stateGlobal[4], stateGlobal[5]);
  if (G != NULL) {
    if(G->GetNrows() != 3) {
      Exception exc("RKTrackRepTime::extrapolateToLine ==> G is not 3x3",__LINE__,__FILE__);
      exc.setFatal();
      throw exc;
    }
    dir = TMatrix(*G) * dir;
  }
  TVector3 lastDir(0,0,0);

  TVector3 poca;
  bool isAtBoundary(false);

  DetPlane startPlane(*(state.getPlane()));
  SharedPlanePtr plane(new DetPlane(point, dir));
  unsigned int iterations(0);
  double charge = getCharge(state);
  double mass = getMass(state);
  double flightTime = 0;

  while(true){
    if(++iterations == maxIt) {
      Exception exc("RKTrackRepTime::extrapolateToPoint ==> extrapolation to point failed, maximum number of iterations reached",__LINE__,__FILE__);
      exc.setFatal();
      throw exc;
    }

    lastStep = step;
    lastDir = dir;

    step = this->Extrap(startPlane, *plane, charge, mass, isAtBoundary, stateGlobal, flightTime, false, NULL, true, stopAtBoundary, maxStep);
    tracklength += step;

    dir.SetXYZ(stateGlobal[3], stateGlobal[4], stateGlobal[5]);
    if (G != NULL) {
      dir = TMatrix(*G) * dir;
    }
    poca.SetXYZ(stateGlobal[0], stateGlobal[1], stateGlobal[2]);

    // check break conditions
    if (stopAtBoundary && isAtBoundary) {
      plane->setON(dir, poca);
      break;
    }

    angle = fabs(dir.Angle((point-poca))-TMath::PiOver2()); // angle between direction and connection to point - 90 deg
    distToPoca = (point-poca).Mag();
    if (angle*distToPoca < 0.1*MINSTEP) break;

    // if lastStep and step have opposite sign, the real normal vector lies somewhere between the last two normal vectors (i.e. the directions)
    // -> try mean value of the two
    if (lastStep*step < 0){
      if (G != NULL) { // after multiplication with G, dir has not length 1 anymore in general
        dir.SetMag(1.);
        lastDir.SetMag(1.);
      }
      dir += lastDir;
      maxStep = 0.5*fabs(lastStep); // make it converge!
    }

    startPlane = *plane;
    plane->setNormal(dir);
  }

  if (fillExtrapSteps) { // now do the full extrapolation with covariance matrix
    // make use of the cache
    lastEndState_.setPlane(plane);
    getStateLocal(lastEndState_, stateGlobal);

    tracklength = extrapolateToPlane(state, plane, false, true);
    lastEndState_.getAuxInfo()(1) = state.getAuxInfo()(1); // Flight time
  }
  else {
    state.setPlane(plane);
    getStateLocal(state, stateGlobal);
    state.getAuxInfo()(1) += flightTime;
  }


  if (debugLvl_ > 0) {
    std::cout << "RKTrackRepTime::extrapolateToPoint(): Reached POCA after " << iterations+1 << " iterations. Distance: " << (point-poca).Mag() << " cm. Angle deviation: " << dir.Angle((point-poca))-TMath::PiOver2() << " rad \n";
  }

  lastEndState_ = state;

  return tracklength;
}


double RKTrackRepTime::extrapolateToCylinder(StateOnPlane& state,
    double radius,
    const TVector3& linePoint,
    const TVector3& lineDirection,
    bool stopAtBoundary,
    bool calcJacobianNoise) const {

  if (debugLvl_ > 0) {
    std::cout << "RKTrackRepTime::extrapolateToCylinder()\n";
  }

  checkCache(state, NULL);

  static const unsigned int maxIt(1000);

  // to 7D
  M1x8 stateGlobal;
  getStateGlobal(state, stateGlobal);

  bool fillExtrapSteps(false);
  if (dynamic_cast<MeasuredStateOnPlane*>(&state) != NULL) {
    fillExtrapSteps = true;
  }
  else if (calcJacobianNoise)
    fillExtrapSteps = true;

  double tracklength(0.), maxStep(1.E99);

  TVector3 dest, pos, dir;

  bool isAtBoundary(false);

  DetPlane startPlane(*(state.getPlane()));
  SharedPlanePtr plane(new DetPlane());
  unsigned int iterations(0);
  double charge = getCharge(state);
  double mass = getMass(state);
  double flightTime = 0;

  while(true){
    if(++iterations == maxIt) {
      Exception exc("RKTrackRepTime::extrapolateToCylinder ==> maximum number of iterations reached",__LINE__,__FILE__);
      exc.setFatal();
      throw exc;
    }

    pos.SetXYZ(stateGlobal[0], stateGlobal[1], stateGlobal[2]);
    dir.SetXYZ(stateGlobal[3], stateGlobal[4], stateGlobal[5]);

    // solve quadratic equation
    TVector3 AO = (pos - linePoint);
    TVector3 AOxAB = (AO.Cross(lineDirection));
    TVector3 VxAB  = (dir.Cross(lineDirection));
    float ab2    = (lineDirection * lineDirection);
    float a      = (VxAB * VxAB);
    float b      = 2 * (VxAB * AOxAB);
    float c      = (AOxAB * AOxAB) - (radius*radius * ab2);
    double arg = b*b - 4.*a*c;
    if(arg < 0) {
      Exception exc("RKTrackRepTime::extrapolateToCylinder ==> cannot solve",__LINE__,__FILE__);
      exc.setFatal();
      throw exc;
    }
    double term = sqrt(arg);
    double k1, k2;
    if (b<0) {
      k1 = (-b + term)/(2.*a);
      k2 = 2.*c/(-b + term);
    }
    else {
      k1 = 2.*c/(-b - term);
      k2 = (-b - term)/(2.*a);
    }

    // select smallest absolute solution -> closest cylinder surface
    double k = k1;
    if (fabs(k2)<fabs(k))
    k = k2;

    if (debugLvl_ > 0) {
      std::cout << "RKTrackRepTime::extrapolateToCylinder(); k = " << k << "\n";
    }

    dest = pos + k * dir;

    plane->setO(dest);
    plane->setUV((dest-linePoint).Cross(lineDirection), lineDirection);

    tracklength += this->Extrap(startPlane, *plane, charge, mass, isAtBoundary, stateGlobal, flightTime, false, NULL, true, stopAtBoundary, maxStep);

    // check break conditions
    if (stopAtBoundary && isAtBoundary) {
      pos.SetXYZ(stateGlobal[0], stateGlobal[1], stateGlobal[2]);
      dir.SetXYZ(stateGlobal[3], stateGlobal[4], stateGlobal[5]);
      plane->setO(pos);
      plane->setUV((pos-linePoint).Cross(lineDirection), lineDirection);
      break;
    }

    if(fabs(k)<MINSTEP) break;

    startPlane = *plane;

  }

  if (fillExtrapSteps) { // now do the full extrapolation with covariance matrix
    // make use of the cache
    lastEndState_.setPlane(plane);
    getStateLocal(lastEndState_, stateGlobal);

    tracklength = extrapolateToPlane(state, plane, false, true);
    lastEndState_.getAuxInfo()(1) = state.getAuxInfo()(1); // Flight time
  }
  else {
    state.setPlane(plane);
    getStateLocal(state, stateGlobal);
    state.getAuxInfo()(1) += flightTime;
  }

  lastEndState_ = state;

  return tracklength;
}



double RKTrackRepTime::extrapolateToCone(StateOnPlane& state,
    double openingAngle,
    const TVector3& conePoint,
    const TVector3& coneDirection,
    bool stopAtBoundary,
    bool calcJacobianNoise) const {

  if (debugLvl_ > 0) {
    std::cout << "RKTrackRepTime::extrapolateToCone()\n";
  }

  checkCache(state, NULL);

  static const unsigned int maxIt(1000);

  // to 7D
  M1x8 stateGlobal;
  getStateGlobal(state, stateGlobal);

  bool fillExtrapSteps(false);
  if (dynamic_cast<MeasuredStateOnPlane*>(&state) != NULL) {
    fillExtrapSteps = true;
  }
  else if (calcJacobianNoise)
    fillExtrapSteps = true;

  double tracklength(0.), maxStep(1.E99);

  TVector3 dest, pos, dir;

  bool isAtBoundary(false);

  DetPlane startPlane(*(state.getPlane()));
  SharedPlanePtr plane(new DetPlane());
  unsigned int iterations(0);
  double charge = getCharge(state);
  double mass = getMass(state);
  double flightTime = 0;

  while(true){
    if(++iterations == maxIt) {
      Exception exc("RKTrackRepTime::extrapolateToCone ==> maximum number of iterations reached",__LINE__,__FILE__);
      exc.setFatal();
      throw exc;
    }

    pos.SetXYZ(stateGlobal[0], stateGlobal[1], stateGlobal[2]);
    dir.SetXYZ(stateGlobal[3], stateGlobal[4], stateGlobal[5]);

    // solve quadratic equation a k^2 + 2 b k + c = 0
    // a = (U . D)^2 - cos^2 alpha * U^2
    // b = (Delta . D) * (U . D) - cos^2 alpha * (U . Delta)
    // c = (Delta . D)^2 - cos^2 alpha * Delta^2
    // Delta = P - V, P track point, U track direction, V cone point, D cone direction, alpha opening angle of cone
    TVector3 cDirection = coneDirection.Unit();
    TVector3 Delta = (pos - conePoint);
    double DirDelta = cDirection * Delta;
    double Delta2 = Delta*Delta;
    double UDir = dir * cDirection;
    double UDelta = dir * Delta;
    double U2 = dir * dir;
    double cosAngle2 = cos(openingAngle)*cos(openingAngle);
    double a = UDir*UDir - cosAngle2*U2;
    double b = UDir*DirDelta - cosAngle2*UDelta;
    double c = DirDelta*DirDelta - cosAngle2*Delta2;
    
    double arg = b*b - a*c;
    if(arg < -1e-9) {
      Exception exc("RKTrackRepTime::extrapolateToCone ==> cannot solve",__LINE__,__FILE__);
      exc.setFatal();
      throw exc;
    } else if(arg < 0) {
      arg = 0;
    }

    double term = sqrt(arg);
    double k1, k2;
    k1 = (-b + term) / a;
    k2 = (-b - term) / a;

    // select smallest absolute solution -> closest cone surface
    double k = k1;
    if(fabs(k2) < fabs(k)) {
      k = k2;
    }

    if (debugLvl_ > 0) {
      std::cout << "RKTrackRepTime::extrapolateToCone(); k = " << k << "\n";
    }

    dest = pos + k * dir;
    // std::cout << "In cone extrapolation ";
    // dest.Print();

    plane->setO(dest);
    plane->setUV((dest-conePoint).Cross(coneDirection), dest-conePoint);

    tracklength += this->Extrap(startPlane, *plane, charge, mass, isAtBoundary, stateGlobal, flightTime, false, NULL, true, stopAtBoundary, maxStep);

    // check break conditions
    if (stopAtBoundary && isAtBoundary) {
      pos.SetXYZ(stateGlobal[0], stateGlobal[1], stateGlobal[2]);
      dir.SetXYZ(stateGlobal[3], stateGlobal[4], stateGlobal[5]);
      plane->setO(pos);
      plane->setUV((pos-conePoint).Cross(coneDirection), pos-conePoint);
      break;
    }

    if(fabs(k)<MINSTEP) break;

    startPlane = *plane;

  }

  if (fillExtrapSteps) { // now do the full extrapolation with covariance matrix
    // make use of the cache
    lastEndState_.setPlane(plane);
    getStateLocal(lastEndState_, stateGlobal);

    tracklength = extrapolateToPlane(state, plane, false, true);
    lastEndState_.getAuxInfo()(1) = state.getAuxInfo()(1); // Flight time
  }
  else {
    state.setPlane(plane);
    getStateLocal(state, stateGlobal);
    state.getAuxInfo()(1) += flightTime;
  }

  lastEndState_ = state;

  return tracklength;
}


double RKTrackRepTime::extrapolateToSphere(StateOnPlane& state,
    double radius,
    const TVector3& point, // center
    bool stopAtBoundary,
    bool calcJacobianNoise) const {

  if (debugLvl_ > 0) {
    std::cout << "RKTrackRepTime::extrapolateToSphere()\n";
  }

  checkCache(state, NULL);

  static const unsigned int maxIt(1000);

  // to 7D
  M1x8 stateGlobal;
  getStateGlobal(state, stateGlobal);

  bool fillExtrapSteps(false);
  if (dynamic_cast<MeasuredStateOnPlane*>(&state) != NULL) {
    fillExtrapSteps = true;
  }
  else if (calcJacobianNoise)
    fillExtrapSteps = true;

  double tracklength(0.), maxStep(1.E99);

  TVector3 dest, pos, dir;

  bool isAtBoundary(false);

  DetPlane startPlane(*(state.getPlane()));
  SharedPlanePtr plane(new DetPlane());
  unsigned int iterations(0);
  double charge = getCharge(state);
  double mass = getMass(state);
  double flightTime = 0;

  while(true){
    if(++iterations == maxIt) {
      Exception exc("RKTrackRepTime::extrapolateToSphere ==> maximum number of iterations reached",__LINE__,__FILE__);
      exc.setFatal();
      throw exc;
    }

    pos.SetXYZ(stateGlobal[0], stateGlobal[1], stateGlobal[2]);
    dir.SetXYZ(stateGlobal[3], stateGlobal[4], stateGlobal[5]);

    // solve quadratic equation
    TVector3 AO = (pos - point);
    double dirAO = dir * AO;
    double arg = dirAO*dirAO - AO*AO + radius*radius;
    if(arg < 0) {
      Exception exc("RKTrackRepTime::extrapolateToSphere ==> cannot solve",__LINE__,__FILE__);
      exc.setFatal();
      throw exc;
    }
    double term = sqrt(arg);
    double k1, k2;
    k1 = -dirAO + term;
    k2 = -dirAO - term;

    // select smallest absolute solution -> closest cylinder surface
    double k = k1;
    if (fabs(k2)<fabs(k))
    k = k2;

    if (debugLvl_ > 0) {
      std::cout << "RKTrackRepTime::extrapolateToSphere(); k = " << k << "\n";
    }

    dest = pos + k * dir;

    plane->setON(dest, dest-point);

    tracklength += this->Extrap(startPlane, *plane, charge, mass, isAtBoundary, stateGlobal, flightTime, false, NULL, true, stopAtBoundary, maxStep);

    // check break conditions
    if (stopAtBoundary && isAtBoundary) {
      pos.SetXYZ(stateGlobal[0], stateGlobal[1], stateGlobal[2]);
      dir.SetXYZ(stateGlobal[3], stateGlobal[4], stateGlobal[5]);
      plane->setON(pos, pos-point);
      break;
    }

    if(fabs(k)<MINSTEP) break;

    startPlane = *plane;

  }

  if (fillExtrapSteps) { // now do the full extrapolation with covariance matrix
    // make use of the cache
    lastEndState_.setPlane(plane);
    getStateLocal(lastEndState_, stateGlobal);

    tracklength = extrapolateToPlane(state, plane, false, true);
    lastEndState_.getAuxInfo()(1) = state.getAuxInfo()(1); // Flight time
  }
  else {
    state.setPlane(plane);
    getStateLocal(state, stateGlobal);
    state.getAuxInfo()(1) += flightTime;
  }

  lastEndState_ = state;

  return tracklength;
}


double RKTrackRepTime::extrapolateBy(StateOnPlane& state,
    double step,
    bool stopAtBoundary,
    bool calcJacobianNoise) const {

  if (debugLvl_ > 0) {
    std::cout << "RKTrackRepTime::extrapolateBy()\n";
  }

  checkCache(state, NULL);

  static const unsigned int maxIt(1000);

  // to 7D
  M1x8 stateGlobal;
  getStateGlobal(state, stateGlobal);

  bool fillExtrapSteps(false);
  if (dynamic_cast<MeasuredStateOnPlane*>(&state) != NULL) {
    fillExtrapSteps = true;
  }
  else if (calcJacobianNoise)
    fillExtrapSteps = true;

  double tracklength(0.);

  TVector3 dest, pos, dir;

  bool isAtBoundary(false);

  DetPlane startPlane(*(state.getPlane()));
  SharedPlanePtr plane(new DetPlane());
  unsigned int iterations(0);
  double charge = getCharge(state);
  double mass = getMass(state);
  double flightTime = 0;

  while(true){
    if(++iterations == maxIt) {
      Exception exc("RKTrackRepTime::extrapolateBy ==> maximum number of iterations reached",__LINE__,__FILE__);
      exc.setFatal();
      throw exc;
    }

    pos.SetXYZ(stateGlobal[0], stateGlobal[1], stateGlobal[2]);
    dir.SetXYZ(stateGlobal[3], stateGlobal[4], stateGlobal[5]);

    dest = pos + 1.5*(step-tracklength) * dir;

    plane->setON(dest, dir);

    tracklength += this->Extrap(startPlane, *plane, charge, mass, isAtBoundary, stateGlobal, flightTime, false, NULL, true, stopAtBoundary, (step-tracklength));

    // check break conditions
    if (stopAtBoundary && isAtBoundary) {
      pos.SetXYZ(stateGlobal[0], stateGlobal[1], stateGlobal[2]);
      dir.SetXYZ(stateGlobal[3], stateGlobal[4], stateGlobal[5]);
      plane->setON(pos, dir);
      break;
    }

    if (fabs(tracklength-step) < MINSTEP) {
      if (debugLvl_ > 0) {
        std::cout << "RKTrackRepTime::extrapolateBy(): reached after " << iterations << " iterations. \n";
      }
      pos.SetXYZ(stateGlobal[0], stateGlobal[1], stateGlobal[2]);
      dir.SetXYZ(stateGlobal[3], stateGlobal[4], stateGlobal[5]);
      plane->setON(pos, dir);
      break;
    }

    startPlane = *plane;

  }

  if (fillExtrapSteps) { // now do the full extrapolation with covariance matrix
    // make use of the cache
    lastEndState_.setPlane(plane);
    getStateLocal(lastEndState_, stateGlobal);

    tracklength = extrapolateToPlane(state, plane, false, true);
    lastEndState_.getAuxInfo()(1) = state.getAuxInfo()(1); // Flight time
  }
  else {
    state.setPlane(plane);
    getStateLocal(state, stateGlobal);
    state.getAuxInfo()(1) += flightTime;
  }

  lastEndState_ = state;

  return tracklength;
}


TVector3 RKTrackRepTime::getPos(const StateOnPlane& state) const {
  M1x8 stateGlobal;
  getStateGlobal(state, stateGlobal);

  return TVector3(stateGlobal[0], stateGlobal[1], stateGlobal[2]);
}


TVector3 RKTrackRepTime::getMom(const StateOnPlane& state) const {
  M1x8 stateGlobal;
  getStateGlobal(state, stateGlobal);

  TVector3 mom(stateGlobal[3], stateGlobal[4], stateGlobal[5]);
  mom.SetMag(getCharge(state)/stateGlobal[6]);
  return mom;
}


void RKTrackRepTime::getPosMom(const StateOnPlane& state, TVector3& pos, TVector3& mom) const {
  M1x8 stateGlobal;
  getStateGlobal(state, stateGlobal);

  pos.SetXYZ(stateGlobal[0], stateGlobal[1], stateGlobal[2]);
  mom.SetXYZ(stateGlobal[3], stateGlobal[4], stateGlobal[5]);
  mom.SetMag(getCharge(state)/stateGlobal[6]);
}

void RKTrackRepTime::getPosMomCov(const MeasuredStateOnPlane& state, TVector3& pos, TVector3& mom, TMatrixDSym& cov) const {
  getPosMom(state, pos, mom);
  cov.ResizeTo(6,6);
  transformPM6(state, *((M6x6*) cov.GetMatrixArray()));
}


TMatrixDSym RKTrackRepTime::get6DCov(const MeasuredStateOnPlane& state) const {
  TMatrixDSym cov(6);
  transformPM6(state, *((M6x6*) cov.GetMatrixArray()));

  return cov;
}


double RKTrackRepTime::getCharge(const StateOnPlane& state) const {

  if (dynamic_cast<const MeasurementOnPlane*>(&state) != NULL) {
    Exception exc("RKTrackRepTime::getCharge - cannot get charge from MeasurementOnPlane",__LINE__,__FILE__);
    exc.setFatal();
    throw exc;
  }

  double pdgCharge( this->getPDGCharge() );

  // return pdgCharge with sign of q/p
  if (state.getState()(0) == 0)
    return pdgCharge;
  return copysign(pdgCharge, state.getState()(0));
}


double RKTrackRepTime::getMomMag(const StateOnPlane& state) const {
  // p = q / qop
  double p = getCharge(state)/state.getState()(0);
  assert (p>=0);
  return p;
}


double RKTrackRepTime::getMomVar(const MeasuredStateOnPlane& state) const {

  if (dynamic_cast<const MeasurementOnPlane*>(&state) != NULL) {
    Exception exc("RKTrackRepTime::getMomVar - cannot get momVar from MeasurementOnPlane",__LINE__,__FILE__);
    exc.setFatal();
    throw exc;
  }

  // p(qop) = q/qop
  // dp/d(qop) = - q / (qop^2)
  // (delta p) = (delta qop) * |dp/d(qop)| = delta qop * |q / (qop^2)|
  // (var p) = (var qop) * q^2 / (qop^4)

  // delta means sigma
  // cov(0,0) is sigma^2

  return state.getCov()(0,0) * pow(getCharge(state), 2)  / pow(state.getState()(0), 4);
}


double RKTrackRepTime::getSpu(const StateOnPlane& state) const {

  if (dynamic_cast<const MeasurementOnPlane*>(&state) != NULL) {
    Exception exc("RKTrackRepTime::getSpu - cannot get spu from MeasurementOnPlane",__LINE__,__FILE__);
    exc.setFatal();
    throw exc;
  }

  const TVectorD& auxInfo = state.getAuxInfo();
  if (auxInfo.GetNrows() == 2
      || auxInfo.GetNrows() == 1) // backwards compatibility with old RKTrackRepTime
    return state.getAuxInfo()(0);
  else
    return 1.;
}

double RKTrackRepTime::getTime(const StateOnPlane& state) const {

  const TVectorD& auxInfo = state.getAuxInfo();
  if (auxInfo.GetNrows() == 2)
    return state.getAuxInfo()(1);
  else
    return 0.;
}


void RKTrackRepTime::calcForwardJacobianAndNoise(const M1x8& startStateGlobal, const DetPlane& startPlane,
					     const M1x8& destStateGlobal, const DetPlane& destPlane) const {

  if (debugLvl_ > 0) {
    std::cout << "RKTrackRepTime::calcForwardJacobianAndNoise " << std::endl;
  }

  if (ExtrapSteps_.size() == 0) {
    Exception exc("RKTrackRepTime::calcForwardJacobianAndNoise ==> cache is empty. Extrapolation must run with a MeasuredStateOnPlane.",__LINE__,__FILE__);
    throw exc;
  }

  // The Jacobians returned from RKutta are transposed.
  TMatrixD jac(TMatrixD::kTransposed, TMatrixD(8, 8, ExtrapSteps_.back().jac_.begin()));
  // FIXME noise matrix resizing is awkward
  TMatrixDSym noise(7, ExtrapSteps_.back().noise_.begin());
  noise.ResizeTo(8, 8);
  for (int i = ExtrapSteps_.size() - 2; i >= 0; --i) {
    TMatrixDSym stepNoise(7, ExtrapSteps_[i].noise_.begin());
    stepNoise.ResizeTo(8, 8);
    noise += stepNoise.Similarity(jac);
    jac *= TMatrixD(TMatrixD::kTransposed, TMatrixD(8, 8, ExtrapSteps_[i].jac_.begin()));
  }

  // Project into 6x6 space.
  M1x3 pTilde = {{startStateGlobal[3], startStateGlobal[4], startStateGlobal[5]}};
  const TVector3& normal = startPlane.getNormal();
  double pTildeW = pTilde[0] * normal.X() + pTilde[1] * normal.Y() + pTilde[2] * normal.Z();
  double spu = pTildeW > 0 ? 1 : -1;
  for (unsigned int i=0; i<3; ++i) {
    pTilde[i] *= spu/pTildeW; // | pTilde * W | has to be 1 (definition of pTilde)
  }
  M6x8 J_pM;
  calcJ_pM_6x8(J_pM, startPlane.getU(), startPlane.getV(), pTilde, spu);
  M8x6 J_Mp;
  calcJ_Mp_8x6(J_Mp, destPlane.getU(), destPlane.getV(), *((M1x3*) &destStateGlobal[3]));
  jac.Transpose(jac); // Because the helper function wants transposed input.
  RKTools::J_pMTTxJ_MMTTxJ_MpTT(J_Mp, *(M8x8 *)jac.GetMatrixArray(),
				J_pM, *(M6x6 *)fJacobian_.GetMatrixArray());
  RKTools::J_MpTxnoise7xJ_Mp(J_Mp, *(M7x7 *)noise.GetMatrixArray(),
			   *(M6x6 *)fNoise_.GetMatrixArray());

  if (debugLvl_ > 0) {
    std::cout << "total jacobian : "; fJacobian_.Print();
    std::cout << "total noise : "; fNoise_.Print();
  }

}


void RKTrackRepTime::getForwardJacobianAndNoise(TMatrixD& jacobian, TMatrixDSym& noise, TVectorD& deltaState) const {

  jacobian.ResizeTo(6,6);
  jacobian = fJacobian_;

  noise.ResizeTo(6,6);
  noise = fNoise_;

  // lastEndState_ = jacobian * lastStartState_  + deltaState
  deltaState.ResizeTo(6);
  // Calculate this without temporaries:
  //deltaState = lastEndState_.getState() - jacobian * lastStartState_.getState()
  deltaState = lastStartState_.getState();
  deltaState *= jacobian;
  deltaState -= lastEndState_.getState();
  deltaState *= -1;

  if (debugLvl_ > 0) {
    std::cout << "delta state : "; deltaState.Print();
  }
}


void RKTrackRepTime::getBackwardJacobianAndNoise(TMatrixD& jacobian, TMatrixDSym& noise, TVectorD& deltaState) const {

  if (debugLvl_ > 0) {
    std::cout << "RKTrackRepTime::getBackwardJacobianAndNoise " << std::endl;
  }

  if (ExtrapSteps_.size() == 0) {
    Exception exc("RKTrackRepTime::getBackwardJacobianAndNoise ==> cache is empty. Extrapolation must run with a MeasuredStateOnPlane.",__LINE__,__FILE__);
    throw exc;
  }

  jacobian.ResizeTo(6,6);
  jacobian = fJacobian_;
  bool status = TDecompLU::InvertLU(jacobian, 0.0);
  if(status == 0){
    Exception e("cannot invert matrix, status = 0", __LINE__,__FILE__);
    e.setFatal();
    throw e;
  }

  noise.ResizeTo(6,6);
  noise = fNoise_;
  noise.Similarity(jacobian);

  // lastStartState_ = jacobian * lastEndState_  + deltaState
  deltaState.ResizeTo(6);
  deltaState = lastStartState_.getState() - jacobian * lastEndState_.getState();
}


std::vector<genfit::MatStep> RKTrackRepTime::getSteps() const {

  // Todo: test

  if (RKSteps_.size() == 0) {
    Exception exc("RKTrackRepTime::getSteps ==> cache is empty.",__LINE__,__FILE__);
    throw exc;
  }

  std::vector<MatStep> retVal;
  retVal.reserve(RKSteps_.size());

  for (unsigned int i = 0; i<RKSteps_.size(); ++i) {
    retVal.push_back(RKSteps_[i].matStep_);
  }

  return retVal;
}


double RKTrackRepTime::getRadiationLenght() const {

  // Todo: test

  if (RKSteps_.size() == 0) {
    Exception exc("RKTrackRepTime::getRadiationLenght ==> cache is empty.",__LINE__,__FILE__);
    throw exc;
  }

  double radLen(0);

  for (unsigned int i = 0; i<RKSteps_.size(); ++i) {
    radLen += RKSteps_.at(i).matStep_.stepSize_ / RKSteps_.at(i).matStep_.materialProperties_.getRadLen();
  }

  return radLen;
}



void RKTrackRepTime::setPosMom(StateOnPlane& state, const TVector3& pos, const TVector3& mom) const {

  if (state.getRep() != this){
    Exception exc("RKTrackRepTime::setPosMom ==> state is defined wrt. another TrackRep",__LINE__,__FILE__);
    throw exc;
  }

  if (dynamic_cast<MeasurementOnPlane*>(&state) != NULL) {
    Exception exc("RKTrackRepTime::setPosMom - cannot set pos/mom of a MeasurementOnPlane",__LINE__,__FILE__);
    exc.setFatal();
    throw exc;
  }

  if (mom.Mag2() == 0) {
    Exception exc("RKTrackRepTime::setPosMom - momentum is 0",__LINE__,__FILE__);
    exc.setFatal();
    throw exc;
  }

  // init auxInfo if that has not yet happened
  TVectorD& auxInfo = state.getAuxInfo();
  if (auxInfo.GetNrows() != 2) {
    bool alreadySet = auxInfo.GetNrows() == 1;  // backwards compatibility: don't overwrite old setting
    auxInfo.ResizeTo(2);
    if (!alreadySet)
      setSpu(state, 1.);
  }

  if (state.getPlane() != NULL && state.getPlane()->distance(pos) < MINSTEP) { // pos is on plane -> do not change plane!

    M1x8 stateGlobal;

    stateGlobal[0] = pos.X();
    stateGlobal[1] = pos.Y();
    stateGlobal[2] = pos.Z();

    stateGlobal[3] = mom.X();
    stateGlobal[4] = mom.Y();
    stateGlobal[5] = mom.Z();

    // normalize dir
    double norm = 1. / sqrt(stateGlobal[3]*stateGlobal[3] + stateGlobal[4]*stateGlobal[4] + stateGlobal[5]*stateGlobal[5]);
    for (unsigned int i=3; i<6; ++i)
      stateGlobal[i] *= norm;

    stateGlobal[6] = getCharge(state) * norm;
    stateGlobal[7] = 0;

    getStateLocal(state, stateGlobal);

  }
  else { // pos is not on plane -> create new plane!

    // TODO: Raise Warning that a new plane has been created!
    SharedPlanePtr plane(new DetPlane(pos, mom));
    state.setPlane(plane);

    TVectorD& state6(state.getState());

    state6(0) = getCharge(state)/mom.Mag(); // q/p
    state6(1) = 0.; // u'
    state6(2) = 0.; // v'
    state6(3) = 0.; // u
    state6(4) = 0.; // v
    state6(5) = 0.; // time

    setSpu(state, 1.);
  }

}


void RKTrackRepTime::setPosMom(StateOnPlane& state, const TVectorD& state6) const {
  if (state6.GetNrows()!=6){
    Exception exc("RKTrackRepTime::setPosMom ==> state has to be 6d (x, y, z, px, py, pz)",__LINE__,__FILE__);
    throw exc;
  }
  setPosMom(state, TVector3(state6(0), state6(1), state6(2)), TVector3(state6(3), state6(4), state6(5)));
}


void RKTrackRepTime::setPosMomErr(MeasuredStateOnPlane& state, const TVector3& pos, const TVector3& mom, const TVector3& posErr, const TVector3& momErr) const {

  // TODO: test!

  setPosMom(state, pos, mom);

  const TVector3& U(state.getPlane()->getU());
  const TVector3& V(state.getPlane()->getV());
  TVector3 W(state.getPlane()->getNormal());

  double pw = mom * W;
  double pu = mom * U;
  double pv = mom * V;

  TMatrixDSym& cov(state.getCov());

  cov(0,0) = pow(getCharge(state), 2) / pow(mom.Mag(), 6) *
             (mom.X()*mom.X() * momErr.X()*momErr.X()+
              mom.Y()*mom.Y() * momErr.Y()*momErr.Y()+
              mom.Z()*mom.Z() * momErr.Z()*momErr.Z());

  cov(1,1) = pow((U.X()/pw - W.X()*pu/(pw*pw)),2.) * momErr.X()*momErr.X() +
             pow((U.Y()/pw - W.Y()*pu/(pw*pw)),2.) * momErr.Y()*momErr.Y() +
             pow((U.Z()/pw - W.Z()*pu/(pw*pw)),2.) * momErr.Z()*momErr.Z();

  cov(2,2) = pow((V.X()/pw - W.X()*pv/(pw*pw)),2.) * momErr.X()*momErr.X() +
             pow((V.Y()/pw - W.Y()*pv/(pw*pw)),2.) * momErr.Y()*momErr.Y() +
             pow((V.Z()/pw - W.Z()*pv/(pw*pw)),2.) * momErr.Z()*momErr.Z();

  cov(3,3) = posErr.X()*posErr.X() * U.X()*U.X() +
             posErr.Y()*posErr.Y() * U.Y()*U.Y() +
             posErr.Z()*posErr.Z() * U.Z()*U.Z();

  cov(4,4) = posErr.X()*posErr.X() * V.X()*V.X() +
             posErr.Y()*posErr.Y() * V.Y()*V.Y() +
             posErr.Z()*posErr.Z() * V.Z()*V.Z();

}




void RKTrackRepTime::setPosMomCov(MeasuredStateOnPlane& state, const TVector3& pos, const TVector3& mom, const TMatrixDSym& cov6x6) const {

  if (cov6x6.GetNcols()!=6 || cov6x6.GetNrows()!=6){
    Exception exc("RKTrackRepTime::setPosMomCov ==> cov has to be 6x6 (x, y, z, px, py, pz)",__LINE__,__FILE__);
    throw exc;
  }

  setPosMom(state, pos, mom); // charge does not change!

  M1x8 stateGlobal;
  getStateGlobal(state, stateGlobal);

  const M6x6& cov6x6_( *((M6x6*) cov6x6.GetMatrixArray()) );

  // FIXME
  transformM6P(cov6x6_, stateGlobal, state);

}

void RKTrackRepTime::setPosMomCov(MeasuredStateOnPlane& state, const TVectorD& state6, const TMatrixDSym& cov6x6) const {

  if (state6.GetNrows()!=6){
    Exception exc("RKTrackRepTime::setPosMomCov ==> state has to be 6d (x, y, z, px, py, pz)",__LINE__,__FILE__);
    throw exc;
  }

  if (cov6x6.GetNcols()!=6 || cov6x6.GetNrows()!=6){
    Exception exc("RKTrackRepTime::setPosMomCov ==> cov has to be 6x6 (x, y, z, px, py, pz)",__LINE__,__FILE__);
    throw exc;
  }

  TVector3 pos(state6(0), state6(1), state6(2));
  TVector3 mom(state6(3), state6(4), state6(5));
  setPosMom(state, pos, mom); // charge does not change!

  M1x8 stateGlobal;
  getStateGlobal(state, stateGlobal);

  const M6x6& cov6x6_( *((M6x6*) cov6x6.GetMatrixArray()) );

  // FIXME
  transformM6P(cov6x6_, stateGlobal, state);

}


void RKTrackRepTime::setChargeSign(StateOnPlane& state, double charge) const {

  if (dynamic_cast<MeasurementOnPlane*>(&state) != NULL) {
    Exception exc("RKTrackRepTime::setChargeSign - cannot set charge of a MeasurementOnPlane",__LINE__,__FILE__);
    exc.setFatal();
    throw exc;
  }

  if (state.getState()(0) * charge < 0) {
    state.getState()(0) *= -1.;
  }
}


void RKTrackRepTime::setSpu(StateOnPlane& state, double spu) const {
  state.getAuxInfo().ResizeTo(2);
  (state.getAuxInfo())(0) = spu;
}

void RKTrackRepTime::setTime(StateOnPlane& state, double time) const {
  state.getAuxInfo().ResizeTo(2);
  (state.getAuxInfo())(1) = time;

  state.getState()(6) = time;
}


void RKTrackRepTime::derive(const double lambda, const M1x3& T,
                            const double E, const double dEdx, const double d2EdxdE,
                            const double B[3],
                            double& dlambda, M1x3& dT, double& dTime,
                            M5x5* pA = 0) const
{
  // Assumes |q| == 1
  const double c = 29.9792458; // speed of light in cm/ns
  const double kappa = c * 1e-5; // speed of light + Tesla
  const double H[3] = { kappa*B[0], kappa*B[1], kappa*B[2] };

  dTime = fabs(lambda) * E / c;  // == 1 / (c * beta)

  // dEdx is positive in our definition, dlambda should have the same
  // sign as lambda hence no minus in the following line, unlike
  // Bugge et al.
  dlambda = E*pow(lambda, 3) * dEdx; /* *q^-2 omitted */
  // Lorentz force
  dT[0] = lambda * (T[1]*H[2] - T[2]*H[1]);
  dT[1] = lambda * (T[2]*H[0] - T[0]*H[2]);
  dT[2] = lambda * (T[0]*H[1] - T[1]*H[0]);

  if (pA) {
    M5x5& A = *pA;
    std::fill(A.begin(), A.end(), 0);
    A(0,0) = 0; A(0,1) =  lambda*H[2]; A(0,2) = -lambda*H[1]; A(0,3) = T[1]*H[2] - T[2]*H[1];
    A(1,0) = -lambda*H[2]; A(1,1) = 0; A(1,2) =  lambda*H[0]; A(1,3) = T[2]*H[0] - T[0]*H[2];
    A(2,0) =  lambda*H[1]; A(2,1) = -lambda*H[0]; A(2,2) = 0; A(2,3) = T[0]*H[1] - T[1]*H[0];

    // (3.12) in Bugge et al., the derivative of (3.11).  The
    // different choice in units doesn't matter (our lambda doesn't
    // contain kappa).  That, or their units are confused, but I don't
    // want to redo the math with their choice.  Simplified, also
    // avoids dividing by zero if dEdx = 0.
    A(3,3) = dlambda/lambda*(3 - pow(lambda*E, -2)) + d2EdxdE;

    // The derivative of dTime by dlambda.
    A(4,3) = copysign(E / c, lambda) - 1/(lambda*fabs(lambda)*E) / c;
  }
}


double RKTrackRepTime::RKstep(const M1x8& stateGlobal, const double h,
                                const MaterialProperties& mat,
                                M1x8& newStateGlobal,
                                M8x8* pJ = 0) const
{
  const double m = TDatabasePDG::Instance()->GetParticle(getPDG())->Mass();
  const double pdgCharge( this->getPDGCharge() );
  // return pdgCharge with sign of q/p
  const double charge = copysign(pdgCharge, stateGlobal[6]);

  M5x5 *pA1 = 0, *pA2 = 0, *pA3 = 0, *pA4 = 0;
  M5x5 A1, A2, A3, A4;
  if (pJ) {
    pA1 = &A1; pA2 = &A2, pA3 = &A3, pA4 = &A4;
  }

  const M1x3 rStart = {{ stateGlobal[0], stateGlobal[1], stateGlobal[2] }};
  const M1x3 TStart = {{ stateGlobal[3], stateGlobal[4], stateGlobal[5] }};
  const double timeStart = stateGlobal[7];

  MaterialEffects::getInstance()->getParticleParameters(getPDG());
  const double lambdaStart = stateGlobal[6];
  const double EStart = hypot(m, charge / lambdaStart);
  const double dEdxStart = MaterialEffects::getInstance()->dEdx(mat, EStart);
  const double d2EdxdEStart = MaterialEffects::getInstance()->d2EdxdE(mat, EStart);
  //std::cout << std::setprecision(7) << EStart << " " << dEdxStart << " " << m << std::endl;

  double BStart[3];
  FieldManager::getInstance()->getFieldVal(rStart.begin(), BStart);

  double dLambda1;
  M1x3 dT1;
  double dTime1;
  derive(lambdaStart, TStart, EStart, dEdxStart, d2EdxdEStart, BStart,
         dLambda1, dT1, dTime1, pA1);

  const double lambda2 = lambdaStart + h/2*dLambda1;
  const double E2 = hypot(m, charge / lambda2);
  const double dEdx2 = MaterialEffects::getInstance()->dEdx(mat, E2);
  const double d2EdxdE2 = MaterialEffects::getInstance()->d2EdxdE(mat, E2);

  const M1x3 T2 = TStart + h/2*dT1;
  const M1x3 rMiddle = rStart + h/2*TStart + h*h/8*dT1;
  //const double time2 = timeStart + h/2*dTime1;
  double BMiddle[3];
  FieldManager::getInstance()->getFieldVal(rMiddle.begin(), BMiddle);

  double dLambda2;
  M1x3 dT2;
  double dTime2;
  derive(lambda2, T2, E2, dEdx2, d2EdxdE2, BMiddle,
         dLambda2, dT2, dTime2, pA2);

  const double lambda3 = lambdaStart + h/2*dLambda2;
  const double E3 = hypot(m, charge / lambda3);
  const double dEdx3 = MaterialEffects::getInstance()->dEdx(mat, E3);
  const double d2EdxdE3 = MaterialEffects::getInstance()->d2EdxdE(mat, E3);

  const M1x3 T3 = TStart + h/2*dT2;
  //const double time3 = timeStart + h/2*dTime2;

  double dLambda3;
  M1x3 dT3;
  double dTime3;
  derive(lambda3, T3, E3, dEdx3, d2EdxdE3, BMiddle,
         dLambda3, dT3, dTime3, pA3);

  const double lambda4 = lambdaStart + h*dLambda3;
  const double E4 = hypot(m, charge / lambda4);
  const double dEdx4 = MaterialEffects::getInstance()->dEdx(mat, E4);
  const double d2EdxdE4 = MaterialEffects::getInstance()->d2EdxdE(mat, E4);

  const M1x3 T4 = TStart + h*dT3;
  const M1x3 rEnd = rStart + h*TStart + h*h/2*dT3;
  //const double time4 = timeStart + h*dTime3;
  double BEnd[3];
  FieldManager::getInstance()->getFieldVal(rEnd.begin(), BEnd);

  double dLambda4;
  M1x3 dT4;
  double dTime4;
  derive(lambda4, T4, E4, dEdx4, d2EdxdE4, BEnd,
         dLambda4, dT4, dTime4, pA4);

  // Put it together ...
  M1x3 rFinal;
  rFinal = rStart + h*TStart + h*h/6*(dT1 + dT2 + dT3);
  M1x3 TFinal;
  TFinal = TStart + h/6 * (dT1 + 2*dT2 + 2*dT3 + dT4);
  const double norm = hypot(hypot(TFinal[0], TFinal[1]), TFinal[2]);
  for (size_t i = 0; i < 3; ++i)
    TFinal[i] /= norm;
  const double lambdaFinal = lambdaStart + h/6 * (dLambda1 + 2*dLambda2 + 2*dLambda3 + dLambda4);
  const double timeFinal = timeStart + h/6 * (dTime1 + 2*dTime2 + 2*dTime3 + dTime4);

  // ... and put it into the final result
  for (size_t i = 0; i < 3; ++i) {
    newStateGlobal[i] = rFinal[i];
    newStateGlobal[i + 3] = TFinal[i];
  }
  newStateGlobal[6] = lambdaFinal;
  newStateGlobal[7] = timeFinal;

  double epsLambda = fabs(dLambda1 - dLambda2 - dLambda3 + dLambda4);
  double epsTime = fabs(dTime1 - dTime2 - dTime3 + dTime4);
  M1x3 epsT;
  for (size_t i = 0; i < 3; ++i)
    epsT[i] = fabs(dT1[i] - dT2[i] - dT3[i] + dT4[i]);
  double eps = std::max(epsLambda, std::max(epsTime,
                                            *std::max_element(epsT.begin(), epsT.end())));

  if (pJ) {
    // Build the 8x8 Jacobian matrix.  We don't keep the row, column
    // corresponding to \Lambda in the notation of Lund loc.cit. as it
    // does not make it into the final covariance matrices of the 7x7
    // states (everything else wouldn't make sense).  We also assume
    // that Lund's C equals 0 (i.e. no field gradients, no material
    // density gradients).
    M8x8 J;
    std::fill(J.begin(), J.end(), 0);
    for (int i = 0; i < 3; ++i) {
      J(i, i) = 1;
      for (int j = 0; j < 4; ++j) {
        J(i, j + 3) = h * (i == j) + h*h/6 * (A1(i, j) + A2(i, j) + A3(i, j));
      }
    }
    for (int i = 0; i < 5; ++i) {
      for (int j = 0; j < 5; ++j) {
        J(i + 3, j + 3) = (i == j) + h/6 * (A1(i, j) + 2*A2(i, j) + 2*A3(i, j) + A4(i, j));
      }
    }

    // Life is a bit miserable: we have to take into account the
    // normalization of T while putting together the final Jacobian.
    M8x8& Jnew = *pJ;
    Jnew = J;
    for (int iRow = 3; iRow < 6; ++iRow) {
      for (int iCol = 3; iCol < 6; ++iCol) {
        Jnew(iRow, iCol) = J(iRow, iCol) / norm;
        // add the derivative of the norm ...
        double sum = 0;
        for (int k = 3; k < 6; ++k) {
          sum += stateGlobal[k] * J(k, iCol);
        }
        Jnew(iRow, iCol) -= stateGlobal[iRow] * sum / norm;
      }
    }
  }

  return h*h*eps;
}



double RKTrackRepTime::RKPropagate(M1x8& stateGlobal,
                                   M8x8* jacobianT,
                                   M1x3& SA,
                                   double S,
                                   const MaterialProperties& mat) const
{
  // The algorithm is
  //  E Lund et al 2009 JINST 4 P04001 doi:10.1088/1748-0221/4/04/P04001
  //  "Track parameter propagation through the application of a new adaptive Runge-Kutta-Nyström method in the ATLAS experiment"
  //  http://inspirehep.net/search?ln=en&ln=en&p=10.1088/1748-0221/4/04/P04001&of=hb&action_search=Search&sf=earliestdate&so=d&rm=&rg=25&sc=0
  // where the transport of the Jacobian is described in
  //   L. Bugge, J. Myrheim  Nucl.Instrum.Meth. 160 (1979) 43-48
  //   "A Fast Runge-kutta Method For Fitting Tracks In A Magnetic Field"
  //   http://inspirehep.net/record/145692
  // and
  //   L. Bugge, J. Myrheim  Nucl.Instrum.Meth. 179 (1981) 365-381
  //   "Tracking And Track Fitting"
  //   http://inspirehep.net/record/160548

  M1x8 oldStateGlobal(stateGlobal);
  M1x8 newStateGlobal;
  M8x8 propJac;
  double est = RKstep(stateGlobal, S, mat, newStateGlobal, jacobianT ? &propJac : 0);
  M8x8 newJacT;
  if (jacobianT) {
    if (0) {
      // Numerically evaluate the Jacobian, compare
      // no science behind these values, I verified that forward and
      // backward propagation yield inverse matrices to good
      // approximation.  In order to avoid bad roundoff errors, the actual
      // step taken is determined below, separately for each direction.
      const double defaultStepX = 1.E-8;
      double stepX;

      M8x8 numJac;

      // Calculate derivative for all three dimensions successively.
      // The algorithm follows the one in TF1::Derivative() :
      //   df(x) = (4 D(h/2) - D(h)) / 3
      // with D(h) = (f(x + h) - f(x - h)) / (2 h).
      //
      // Could perhaps do better by also using f(x) which would be stB.
      M1x8 rightShort, rightFull;
      M1x8 leftShort, leftFull;
      for (size_t i = 0; i < 8; ++i) {
        {
          M1x8 stateCopy(stateGlobal);
          double temp = stateCopy[i] + defaultStepX / 2;
          // Find the actual size of the step, which will differ from
          // defaultStepX due to roundoff.  This is the step-size we will
          // use for this direction.  Idea taken from Numerical Recipes,
          // 3rd ed., section 5.7.
          //
          // Note that if a number is exactly representable, it's double
          // will also be exact.  Outside denormals, this also holds for
          // halving.  Unless the exponent changes (which it only will in
          // the vicinity of zero) adding or subtracing doesn't make a
          // difference.
          //
          // We determine the roundoff error for the half-step.  If this
          // is exactly representable, the full step will also be.
          stepX = 2 * (temp - stateCopy[i]);
          stateCopy[i] = temp;
          RKstep(stateCopy, S, mat, rightShort, 0);
        }
        {
          M1x8 stateCopy(stateGlobal);
          stateCopy[i] -= stepX / 2;
          RKstep(stateCopy, S, mat, leftShort, 0);
        }
        {
          M1x8 stateCopy(stateGlobal);
          stateCopy[i] += stepX;
          RKstep(stateCopy, S, mat, rightFull, 0);
        }
        {
          M1x8 stateCopy(stateGlobal);
          stateCopy[i] -= stepX;
          RKstep(stateCopy, S, mat, leftFull, 0);
        }

        // Calculate the derivatives for the individual components of
        // the track parameters.
        for (size_t j = 0; j < 8; ++j) {
          double derivFull = (rightFull[j] - leftFull[j]) / 2 / stepX;
          double derivShort = (rightShort[j] - leftShort[j]) / stepX;

          numJac(j, i) = 1./3.*(4*derivShort - derivFull);
        }
      }
      std::cout << "S = " << S << " semianalytical ";
      propJac.print();
      std::cout << "numerical ";
      numJac.print();
      //propJac = numJac;
    }

    for (int i = 0; i < 8; ++i) {
      for (int j = 0; j < 8; ++j) {
        double sum = 0;
        for (int k = 0; k < 8; ++k) {
          sum += propJac(i, k) * (*jacobianT)(j, k);
        }
        newJacT(j, i) = sum;
      }
    }
  }

  static const double DLT ( .0002 );           // max. deviation for approximation-quality test

  double EST = est; // FIXME : why over S?

  SA[0] = newStateGlobal[3] - oldStateGlobal[3];
  SA[1] = newStateGlobal[4] - oldStateGlobal[4];
  SA[2] = newStateGlobal[5] - oldStateGlobal[5];

  std::copy(newStateGlobal.begin(), newStateGlobal.end(), stateGlobal.begin());

  if (jacobianT) {
    std::copy(newJacT.begin(), newJacT.end(), jacobianT->begin());
  }
  /*
  std::cout << S << std::endl;
  oldState7.print();
  newState7.print();
  */
  if (debugLvl_ > 0) {
    std::cout << "    RKTrackRepTime::RKPropagate. Step = "<< S << "; quality EST = " << EST  << " \n";
  }

  // Prevent the step length increase from getting too large, this is
  // just the point where it becomes 10.
  if (EST < DLT*1e-4)
    return 10;

  // Step length increase for a fifth order Runge-Kutta, see e.g. 17.2
  // in Numerical Recipes.  FIXME: move to caller.
  return pow(DLT/EST, 1./4.);
}



void RKTrackRepTime::initArrays() const {
  std::fill(noiseProjection_.begin(), noiseProjection_.end(), 0);
  for (unsigned int i=0; i<7; ++i) // initialize as diagonal matrix
    noiseProjection_[i*8] = 1;

  fJacobian_.UnitMatrix();
  fNoise_.Zero();
  limits_.reset();

  RKSteps_.reserve(100);
  ExtrapSteps_.reserve(100);

  lastStartState_.getAuxInfo().ResizeTo(2);
  lastEndState_.getAuxInfo().ResizeTo(2);
}


void RKTrackRepTime::getStateGlobal(const StateOnPlane& state, M1x8& stateGlobal) const {

  if (dynamic_cast<const MeasurementOnPlane*>(&state) != NULL) {
    Exception exc("RKTrackRepTime::getState7 - cannot get pos or mom from a MeasurementOnPlane",__LINE__,__FILE__);
    exc.setFatal();
    throw exc;
  }

  const TVector3& U(state.getPlane()->getU());
  const TVector3& V(state.getPlane()->getV());
  const TVector3& O(state.getPlane()->getO());
  const TVector3& W(state.getPlane()->getNormal());

  assert(state.getState().GetNrows() == 6);
  const double* state6 = state.getState().GetMatrixArray();

  double spu = getSpu(state);

  stateGlobal[0] = O.X() + state6[3]*U.X() + state6[4]*V.X(); // x
  stateGlobal[1] = O.Y() + state6[3]*U.Y() + state6[4]*V.Y(); // y
  stateGlobal[2] = O.Z() + state6[3]*U.Z() + state6[4]*V.Z(); // z

  stateGlobal[3] = spu * (W.X() + state6[1]*U.X() + state6[2]*V.X()); // a_x
  stateGlobal[4] = spu * (W.Y() + state6[1]*U.Y() + state6[2]*V.Y()); // a_y
  stateGlobal[5] = spu * (W.Z() + state6[1]*U.Z() + state6[2]*V.Z()); // a_z

  // normalize dir
  double norm = 1. / sqrt(stateGlobal[3]*stateGlobal[3] + stateGlobal[4]*stateGlobal[4] + stateGlobal[5]*stateGlobal[5]);
  for (unsigned int i=3; i<6; ++i) stateGlobal[i] *= norm;

  stateGlobal[6] = state6[0]; // q/p
  stateGlobal[7] = state6[5]; // time
}


void RKTrackRepTime::getStateLocal(StateOnPlane& state, const M1x8& stateGlobal) const {

  // state5: (q/p, u', v'. u, v, time)

  double spu(1.);

  const TVector3& O(state.getPlane()->getO());
  const TVector3& U(state.getPlane()->getU());
  const TVector3& V(state.getPlane()->getV());
  const TVector3& W(state.getPlane()->getNormal());

  // force A to be in normal direction and set spu accordingly
  double AtW( stateGlobal[3]*W.X() + stateGlobal[4]*W.Y() + stateGlobal[5]*W.Z() );
  if (AtW < 0.) {
    //fDir *= -1.;
    //AtW *= -1.;
    spu = -1.;
  }

  double* state6 = state.getState().GetMatrixArray();

  state6[0] = stateGlobal[6]; // q/p
  state6[1] = (stateGlobal[3]*U.X() + stateGlobal[4]*U.Y() + stateGlobal[5]*U.Z()) / AtW; // u' = (A * U) / (A * W)
  state6[2] = (stateGlobal[3]*V.X() + stateGlobal[4]*V.Y() + stateGlobal[5]*V.Z()) / AtW; // v' = (A * V) / (A * W)
  state6[3] = ((stateGlobal[0]-O.X())*U.X() +
               (stateGlobal[1]-O.Y())*U.Y() +
               (stateGlobal[2]-O.Z())*U.Z()); // u = (pos - O) * U
  state6[4] = ((stateGlobal[0]-O.X())*V.X() +
               (stateGlobal[1]-O.Y())*V.Y() +
               (stateGlobal[2]-O.Z())*V.Z()); // v = (pos - O) * V
  state6[5] = stateGlobal[7]; // time

  setSpu(state, spu);

}



void RKTrackRepTime::transformPM8(const MeasuredStateOnPlane& state,
                                  M8x8& out8x8) const {

  // get vectors and aux variables
  const TVector3& U(state.getPlane()->getU());
  const TVector3& V(state.getPlane()->getV());
  const TVector3& W(state.getPlane()->getNormal());

  const TVectorD& state5(state.getState());
  double spu(getSpu(state));

  M1x3 pTilde;
  pTilde[0] = spu * (W.X() + state5(1)*U.X() + state5(2)*V.X()); // a_x
  pTilde[1] = spu * (W.Y() + state5(1)*U.Y() + state5(2)*V.Y()); // a_y
  pTilde[2] = spu * (W.Z() + state5(1)*U.Z() + state5(2)*V.Z()); // a_z

  M6x8 J_pM;
  calcJ_pM_6x8(J_pM, U, V, pTilde, spu);

  // since the Jacobian contains a lot of zeros, and the resulting cov has to be symmetric,
  // the multiplication can be done much faster directly on array level
  // out = J_pM^T * in5x5 * J_pM
  const M6x6& in6x6_ = *((M6x6*) state.getCov().GetMatrixArray());
  RKTools::J_pMTxcov6xJ_pM(J_pM, in6x6_, out8x8);
}


void RKTrackRepTime::calcJ_pM_6x8(M6x8& J_pM, const TVector3& U, const TVector3& V, const M1x3& pTilde, double spu) const {
  /*if (debugLvl_ > 1) {
    std::cout << "RKTrackRepTime::calcJ_pM_6x8 \n";
    std::cout << "  U = "; U.Print();
    std::cout << "  V = "; V.Print();
    std::cout << "  pTilde = "; RKTools::printDim(pTilde, 3,1);
    std::cout << "  spu = " << spu << "\n";
  }*/

  std::fill(J_pM.begin(), J_pM.end(), 0);

  const double pTildeMag = sqrt(pTilde[0]*pTilde[0] + pTilde[1]*pTilde[1] + pTilde[2]*pTilde[2]);
  const double pTildeMag2 = pTildeMag*pTildeMag;

  const double utpTildeOverpTildeMag2 = (U.X()*pTilde[0] + U.Y()*pTilde[1] + U.Z()*pTilde[2]) / pTildeMag2;
  const double vtpTildeOverpTildeMag2 = (V.X()*pTilde[0] + V.Y()*pTilde[1] + V.Z()*pTilde[2]) / pTildeMag2;

  //J_pM matrix is d(x,y,z,ax,ay,az,q/p) / d(q/p,u',v',u,v)   (out is 7x7)

   // d(x,y,z)/d(u)
  J_pM(3,0) = U.X();
  J_pM(3,1) = U.Y();
  J_pM(3,2) = U.Z();
  // d(x,y,z)/d(v)
  J_pM(4,0) = V.X();
  J_pM(4,1) = V.Y();
  J_pM(4,2) = V.Z();
  // d(q/p)/d(q/p)
  J_pM(0, 6) = 1.;
  // d(ax,ay,az)/d(u')
  double fact = spu / pTildeMag;
  J_pM(1,3) = fact * ( U.X() - pTilde[0]*utpTildeOverpTildeMag2 ); // [1][3]
  J_pM(1,4) = fact * ( U.Y() - pTilde[1]*utpTildeOverpTildeMag2 ); // [1][4]
  J_pM(1,5) = fact * ( U.Z() - pTilde[2]*utpTildeOverpTildeMag2 ); // [1][5]
  // d(ax,ay,az)/d(v')
  J_pM(2,3) = fact * ( V.X() - pTilde[0]*vtpTildeOverpTildeMag2 ); // [2][3]
  J_pM(2,4) = fact * ( V.Y() - pTilde[1]*vtpTildeOverpTildeMag2 ); // [2][4]
  J_pM(2,5) = fact * ( V.Z() - pTilde[2]*vtpTildeOverpTildeMag2 ); // [2][5]
  // d(time)/d(time)
  J_pM(5,7) = 1.;

  /*if (debugLvl_ > 1) {
    std::cout << "  J_pM = "; J_pM.print()
  }*/
}


void RKTrackRepTime::transformPM6(const MeasuredStateOnPlane& state,
                                  M6x6& out6x6) const {

  // get vectors and aux variables
  const TVector3& U(state.getPlane()->getU());
  const TVector3& V(state.getPlane()->getV());
  const TVector3& W(state.getPlane()->getNormal());

  const TVectorD& state5(state.getState());
  double spu(getSpu(state));

  M1x3 pTilde;
  pTilde[0] = spu * (W.X() + state5(1)*U.X() + state5(2)*V.X()); // a_x
  pTilde[1] = spu * (W.Y() + state5(1)*U.Y() + state5(2)*V.Y()); // a_y
  pTilde[2] = spu * (W.Z() + state5(1)*U.Z() + state5(2)*V.Z()); // a_z

  const double pTildeMag = sqrt(pTilde[0]*pTilde[0] + pTilde[1]*pTilde[1] + pTilde[2]*pTilde[2]);
  const double pTildeMag2 = pTildeMag*pTildeMag;

  const double utpTildeOverpTildeMag2 = (U.X()*pTilde[0] + U.Y()*pTilde[1] + U.Z()*pTilde[2]) / pTildeMag2;
  const double vtpTildeOverpTildeMag2 = (V.X()*pTilde[0] + V.Y()*pTilde[1] + V.Z()*pTilde[2]) / pTildeMag2;

  //J_pM matrix is d(x,y,z,px,py,pz) / d(q/p,u',v',u,v)       (out is 6x6)

  const double qop = state5(0);
  const double p = getCharge(state)/qop; // momentum

  M6x6 J_pM_6x6;
  std::fill(J_pM_6x6.begin(), J_pM_6x6.end(), 0);

  // d(px,py,pz)/d(q/p)
  double fact = -1. * p / (pTildeMag * qop);
  J_pM_6x6(0,3) = fact * pTilde[0];
  J_pM_6x6(0,4) = fact * pTilde[1];
  J_pM_6x6(0,5) = fact * pTilde[2];
  // d(px,py,pz)/d(u')
  fact = p * spu / pTildeMag;
  J_pM_6x6(1,3) = fact * ( U.X() - pTilde[0]*utpTildeOverpTildeMag2 );
  J_pM_6x6(1,4) = fact * ( U.Y() - pTilde[1]*utpTildeOverpTildeMag2 );
  J_pM_6x6(1,5) = fact * ( U.Z() - pTilde[2]*utpTildeOverpTildeMag2 );
  // d(px,py,pz)/d(v')
  J_pM_6x6(2,3) = fact * ( V.X() - pTilde[0]*vtpTildeOverpTildeMag2 );
  J_pM_6x6(2,4) = fact * ( V.Y() - pTilde[1]*vtpTildeOverpTildeMag2 );
  J_pM_6x6(2,5) = fact * ( V.Z() - pTilde[2]*vtpTildeOverpTildeMag2 );
  // d(x,y,z)/d(u)
  J_pM_6x6(3,0) = U.X();
  J_pM_6x6(3,1) = U.Y();
  J_pM_6x6(3,2) = U.Z();
  // d(x,y,z)/d(v)
  J_pM_6x6(4,0) = V.X();
  J_pM_6x6(4,1) = V.Y();
  J_pM_6x6(4,2) = V.Z();

  // do the transformation
  // out = J_pM^T * in6x6 * J_pM
  const M6x6& in6x6_ = *((M6x6*) state.getCov().GetMatrixArray());
  RKTools::J_pMTxcov6xJ_pM(J_pM_6x6, in6x6_, out6x6);

}


void RKTrackRepTime::transformM8P(const M8x8& in8x8,
                              const M1x8& stateGlobal,
                              MeasuredStateOnPlane& state) const { // plane must already be set!

  // get vectors and aux variables
  const TVector3& U(state.getPlane()->getU());
  const TVector3& V(state.getPlane()->getV());

  M1x3& A = *((M1x3*) &stateGlobal[3]);

  M8x6 J_Mp;
  calcJ_Mp_8x6(J_Mp, U, V, A);

  // since the Jacobian contains a lot of zeros, and the resulting cov has to be symmetric,
  // the multiplication can be done much faster directly on array level
  // out6x6 = J_Mp^T * in * J_Mp
  M6x6& out6x6_ = *((M6x6*) state.getCov().GetMatrixArray());
  RKTools::J_MpTxcov8xJ_Mp(J_Mp, in8x8, out6x6_);

}


void RKTrackRepTime::calcJ_Mp_8x6(M8x6& J_Mp, const TVector3& U, const TVector3& V, const M1x3& A) const {

  /*if (debugLvl_ > 1) {
    std::cout << "RKTrackRepTime::calcJ_Mp_7x5 \n";
    std::cout << "  U = "; U.Print();
    std::cout << "  V = "; V.Print();
    std::cout << "  A = "; RKTools::printDim(A, 3,1);
  }*/

  std::fill(J_Mp.begin(), J_Mp.end(), 0);

  TVector3 W = U.Cross(V);
  const double AtU = A[0]*U.X() + A[1]*U.Y() + A[2]*U.Z();
  const double AtV = A[0]*V.X() + A[1]*V.Y() + A[2]*V.Z();
  const double AtW = A[0]*W.X() + A[1]*W.Y() + A[2]*W.Z();

  // J_Mp matrix is d(q/p,u',v',u,v,time) / d(x,y,z,ax,ay,az,q/p,time)   (in is 8x8)

  // d(u')/d(ax,ay,az)
  double fact = 1./(AtW*AtW);
  J_Mp(3,1) = fact * (U.X()*AtW - W.X()*AtU);
  J_Mp(4,1) = fact * (U.Y()*AtW - W.Y()*AtU);
  J_Mp(5,1) = fact * (U.Z()*AtW - W.Z()*AtU);
  // d(v')/d(ax,ay,az)
  J_Mp(3,2) = fact * (V.X()*AtW - W.X()*AtV);
  J_Mp(4,2) = fact * (V.Y()*AtW - W.Y()*AtV);
  J_Mp(5,2) = fact * (V.Z()*AtW - W.Z()*AtV);
  // d(q/p)/d(q/p)
  J_Mp(6,0) = 1.;
  //d(u)/d(x,y,z)
  J_Mp(0,3) = U.X();
  J_Mp(1,3) = U.Y();
  J_Mp(2,3) = U.Z();
  //d(v)/d(x,y,z)
  J_Mp(0,4) = V.X();
  J_Mp(1,4) = V.Y();
  J_Mp(2,4) = V.Z();
  //d(time)/d(time)
  J_Mp(7,5) = 1.;

  /*if (debugLvl_ > 1) {
    std::cout << "  J_Mp_7x5_ = "; RKTools::printDim(J_Mp, 7,5);
  }*/

}


void RKTrackRepTime::transformM6P(const M6x6& in6x6,
                                  const M1x8& stateGlobal,
                                  MeasuredStateOnPlane& state) const { // plane and charge must already be set!

  // get vectors and aux variables
  const TVector3& U(state.getPlane()->getU());
  const TVector3& V(state.getPlane()->getV());
  const TVector3& W(state.getPlane()->getNormal());

  const double AtU = stateGlobal[3]*U.X() + stateGlobal[4]*U.Y() + stateGlobal[5]*U.Z();
  const double AtV = stateGlobal[3]*V.X() + stateGlobal[4]*V.Y() + stateGlobal[5]*V.Z();
  const double AtW = stateGlobal[3]*W.X() + stateGlobal[4]*W.Y() + stateGlobal[5]*W.Z();

  // J_Mp matrix is d(q/p,u',v',u,v,time) / d(x,y,z,px,py,pz)       (in is 6x6)

  const double qop = stateGlobal[6];
  const double p = getCharge(state)/qop; // momentum

  M6x6 J_Mp_6x6;
  std::fill(J_Mp_6x6.begin(), J_Mp_6x6.end(), 0);

  //d(u)/d(x,y,z)
  J_Mp_6x6(0,3) = U.X();
  J_Mp_6x6(1,3) = U.Y();
  J_Mp_6x6(2,3) = U.Z();
  //d(v)/d(x,y,z)
  J_Mp_6x6(0,4) = V.X();
  J_Mp_6x6(1,4) = V.Y();
  J_Mp_6x6(2,4) = V.Z();
  // d(q/p)/d(px,py,pz)
  double fact = (-1.) * qop / p;
  J_Mp_6x6(3,0) = fact * stateGlobal[3];
  J_Mp_6x6(4,0) = fact * stateGlobal[4];
  J_Mp_6x6(5,0) = fact * stateGlobal[5];
  // d(u')/d(px,py,pz)
  fact = 1./(p*AtW*AtW);
  J_Mp_6x6(3,1) = fact * (U.X()*AtW - W.X()*AtU);
  J_Mp_6x6(4,1) = fact * (U.Y()*AtW - W.Y()*AtU);
  J_Mp_6x6(5,1) = fact * (U.Z()*AtW - W.Z()*AtU);
  // d(v')/d(px,py,pz)
  J_Mp_6x6(3,2) = fact * (V.X()*AtW - W.X()*AtV);
  J_Mp_6x6(4,2) = fact * (V.Y()*AtW - W.Y()*AtV);
  J_Mp_6x6(5,2) = fact * (V.Z()*AtW - W.Z()*AtV);

  // do the transformation
  // out6x6 = J_Mp^T * in * J_Mp
  M6x6& out6x6 = *((M6x6*) state.getCov().GetMatrixArray());
  RKTools::J_MpTxcov6xJ_Mp(J_Mp_6x6, in6x6, out6x6);

  out6x6[5*6 + 5] = 1;
}


//
// Runge-Kutta method for tracking a particles through a magnetic field.
// Uses Nystroem algorithm (See Handbook Nat. Bur. of Standards, procedure 25.5.20)
// in the way described in
//  E Lund et al 2009 JINST 4 P04001 doi:10.1088/1748-0221/4/04/P04001
//  "Track parameter propagation through the application of a new adaptive Runge-Kutta-Nyström method in the ATLAS experiment"
//  http://inspirehep.net/search?ln=en&ln=en&p=10.1088/1748-0221/4/04/P04001&of=hb&action_search=Search&sf=earliestdate&so=d&rm=&rg=25&sc=0
//
// Input parameters:
//    SU     - plane parameters
//    SU[0]  - direction cosines normal to surface Ex
//    SU[1]  -          -------                    Ey
//    SU[2]  -          -------                    Ez; Ex*Ex+Ey*Ey+Ez*Ez=1
//    SU[3]  - distance to surface from (0,0,0) > 0 cm
//
//    state7 - initial parameters (coordinates(cm), direction,
//             charge/momentum (Gev-1)
//    cov      and derivatives this parameters  (7x7)
//
//    X         Y         Z         Ax        Ay        Az        q/P
//    state7[0] state7[1] state7[2] state7[3] state7[4] state7[5] state7[6]
//
//    dX/dp     dY/dp     dZ/dp     dAx/dp    dAy/dp    dAz/dp    d(q/P)/dp
//    cov[ 0]   cov[ 1]   cov[ 2]   cov[ 3]   cov[ 4]   cov[ 5]   cov[ 6]               d()/dp1
//
//    cov[ 7]   cov[ 8]   cov[ 9]   cov[10]   cov[11]   cov[12]   cov[13]               d()/dp2
//    ............................................................................    d()/dpND
//
// Authors: R.Brun, M.Hansroul, V.Perevoztchikov (Geant3)
//
bool RKTrackRepTime::RKutta(const M1x4& SU,
                        const DetPlane& plane,
                        double charge,
                        double mass,
                        M1x8& stateGlobal,
                        M8x8* jacobianT,
                        double& coveredDistance,
                        double& flightTime,
                        bool& checkJacProj,
                        M7x7& noiseProjection,
                        StepLimits& limits,
                        bool onlyOneStep) const {

  // limits, check-values, etc. Can be tuned!
  static const double Wmax           ( 3000. );           // max. way allowed [cm]
  static const double AngleMax       ( 6.3 );           // max. total angle change of momentum. Prevents extrapolating a curler round and round if no active plane is found.
  static const double Pmin           ( 4.E-3 );           // minimum momentum for propagation [GeV]
  static const unsigned int maxNumIt ( 1000 );    // maximum number of iterations in main loop
  // Aux parameters
  M1x3&   R          ( *((M1x3*) &stateGlobal[0]) );  // Start coordinates  [cm]  (x,  y,  z)
  M1x3&   A          ( *((M1x3*) &stateGlobal[3]) );  // Start directions         (ax, ay, az);   ax^2+ay^2+az^2=1
  M1x3    SA         = {{0.,0.,0.}};             // Start directions derivatives dA/S
  double  Way        ( 0. );                     // Sum of absolute values of all extrapolation steps [cm]
  double  momentum   ( fabs(charge/stateGlobal[6]) ); // momentum [GeV]
  double  relMomLoss ( 0 );                      // relative momentum loss in RKutta
  double  deltaAngle ( 0. );                     // total angle by which the momentum has changed during extrapolation
  double  An(0), S(0), Sl(0), CBA(0);

  if (debugLvl_ > 0) {
    std::cout << "RKTrackRepTime::RKutta \n";
    std::cout << "position: "; TVector3(R[0], R[1], R[2]).Print();
    std::cout << "direction: "; TVector3(A[0], A[1], A[2]).Print();
    std::cout << "momentum: " << momentum << " GeV\n";
    std::cout << "destination: "; plane.Print();
  }

  checkJacProj = false;

  // check momentum
  if(momentum < Pmin){
    std::ostringstream sstream;
    sstream << "RKTrackRepTime::RKutta ==> momentum too low: " << momentum*1000. << " MeV";
    Exception exc(sstream.str(),__LINE__,__FILE__);
    exc.setFatal();
    throw exc;
  }

  unsigned int counter(0);

  // Step estimation (signed)
  MaterialProperties matForStep;
  S = estimateStep(stateGlobal, SU, plane, charge, relMomLoss, limits, matForStep);

  //
  // Main loop of Runge-Kutta method
  //
  while (fabs(S) >= MINSTEP || counter == 0) {

    if(++counter > maxNumIt){
      Exception exc("RKTrackRepTime::RKutta ==> maximum number of iterations exceeded",__LINE__,__FILE__);
      exc.setFatal();
      throw exc;
    }

    if (debugLvl_ > 0) {
      std::cout << "------ RKutta main loop nr. " << counter-1 << " ------\n";
    }

    M1x3 ABefore = {{ A[0], A[1], A[2] }};
    RKPropagate(stateGlobal, jacobianT, SA, S, matForStep); // the actual Runge Kutta propagation

    // update paths
    coveredDistance += S;       // add stepsize to way (signed)
    Way  += fabs(S);

    double beta = 1/hypot(1, mass*stateGlobal[6]/charge);
    flightTime += S / beta / 29.9792458; // in ns
    std::cout << counter << "\t " << flightTime << " " << S / beta / 29.9792458 << " " << stateGlobal[7] << std::endl;

    // check way limit
    if(Way > Wmax){
      std::ostringstream sstream;
      sstream << "RKTrackRepTime::RKutta ==> Total extrapolation length is longer than length limit : " << Way << " cm !";
      Exception exc(sstream.str(),__LINE__,__FILE__);
      exc.setFatal();
      throw exc;
    }

    if (onlyOneStep) return(true);

    // if stepsize has been limited by material, break the loop and return. No linear extrapolation!
    if (limits.getLowestLimit().first == stp_momLoss) {
      if (debugLvl_ > 0) {
        std::cout<<" momLossExceeded -> return(true); \n";
      }
      return(true);
    }

    // if stepsize has been limited by material boundary, break the loop and return. No linear extrapolation!
    if (limits.getLowestLimit().first == stp_boundary) {
      if (debugLvl_ > 0) {
        std::cout<<" at boundary -> return(true); \n";
      }
      return(true);
    }


    // estimate Step for next loop or linear extrapolation
    Sl = S; // last S used
    limits.removeLimit(stp_fieldCurv);
    limits.removeLimit(stp_momLoss);
    limits.removeLimit(stp_boundary);
    limits.removeLimit(stp_plane);
    S = estimateStep(stateGlobal, SU, plane, charge, relMomLoss, limits, matForStep);

    if (limits.getLowestLimit().first == stp_plane &&
        fabs(S) < MINSTEP) {
      if (debugLvl_ > 0) {
        std::cout<<" (at Plane && fabs(S) < MINSTEP) -> break and do linear extrapolation \n";
      }
      break;
    }
    if (limits.getLowestLimit().first == stp_momLoss &&
        fabs(S) < MINSTEP) {
      if (debugLvl_ > 0) {
        std::cout<<" (momLossExceeded && fabs(S) < MINSTEP) -> return(true), no linear extrapolation; \n";
      }
      RKSteps_.erase(RKSteps_.end()-1);
      --RKStepsFXStop_;
      return(true); // no linear extrapolation!
    }

    // check if total angle is bigger than AngleMax. Can happen if a curler should be fitted and it does not hit the active area of the next plane.
    double arg = ABefore[0]*A[0] + ABefore[1]*A[1] + ABefore[2]*A[2];
    arg = arg > 1 ? 1 : arg;
    arg = arg < -1 ? -1 : arg;
    deltaAngle += acos(arg);
    if (fabs(deltaAngle) > AngleMax){
      std::ostringstream sstream;
      sstream << "RKTrackRepTime::RKutta ==> Do not get to an active plane! Already extrapolated " << deltaAngle * 180 / TMath::Pi() << "°.";
      Exception exc(sstream.str(),__LINE__,__FILE__);
      exc.setFatal();
      throw exc;
    }

    // check if we went back and forth multiple times -> we don't come closer to the plane!
    if (counter > 3){
      if (S                            *RKSteps_.at(counter-1).matStep_.stepSize_ < 0 &&
          RKSteps_.at(counter-1).matStep_.stepSize_*RKSteps_.at(counter-2).matStep_.stepSize_ < 0 &&
          RKSteps_.at(counter-2).matStep_.stepSize_*RKSteps_.at(counter-3).matStep_.stepSize_ < 0){
        Exception exc("RKTrackRepTime::RKutta ==> Do not get closer to plane!",__LINE__,__FILE__);
        exc.setFatal();
        throw exc;
      }
    }

  } //end of main loop


  //
  // linear extrapolation to plane
  //
  if (limits.getLowestLimit().first == stp_plane) {

    if (fabs(Sl) > 0.001*MINSTEP){
      if (debugLvl_ > 0) {
        std::cout << " RKutta - linear extrapolation to surface " << S <<"\n";
      }
      Sl = 1./Sl;        // Sl = inverted last Stepsize Sl

      // normalize SA
      SA[0]*=Sl;  SA[1]*=Sl;  SA[2]*=Sl; // SA/Sl = delta A / delta way; local derivative of A with respect to the length of the way

      // calculate A
      A[0] += SA[0]*S;    // S  = distance to surface
      A[1] += SA[1]*S;    // A = A + S * SA*Sl
      A[2] += SA[2]*S;

      // normalize A
      CBA = 1./sqrt(A[0]*A[0]+A[1]*A[1]+A[2]*A[2]);  // 1/|A|
      A[0] *= CBA; A[1] *= CBA; A[2] *= CBA;

      R[0] += S*(A[0]-0.5*S*SA[0]);    // R = R + S*(A - 0.5*S*SA); approximation for final point on surface
      R[1] += S*(A[1]-0.5*S*SA[1]);
      R[2] += S*(A[2]-0.5*S*SA[2]);


      coveredDistance += S;
      Way  += fabs(S);

      double beta = 1/hypot(1, mass*stateGlobal[6]/charge);
      flightTime += S / beta / 29.9792458; // in ns;
    }
    else if (debugLvl_ > 0)  {
      std::cout << " RKutta - last stepsize too small -> can't do linear extrapolation! \n";
    }

    //
    // Project Jacobian of extrapolation onto destination plane
    //
    if (jacobianT != NULL) {

      // projected jacobianT
      // x x x x x x 0
      // x x x x x x 0
      // x x x x x x 0
      // x x x x x x 0
      // x x x x x x 0
      // x x x x x x 0
      // x x x x x x 1

      if (checkJacProj && RKSteps_.size()>0){
        Exception exc("RKTrackRepTime::Extrap ==> covariance is projected onto destination plane again",__LINE__,__FILE__);
        throw exc;
      }

      if (debugLvl_ > 0) {
        //std::cout << "  Jacobian^T of extrapolation before Projection:\n";
        //RKTools::printDim(*jacobianT, 7,7);
        std::cout << "  Project Jacobian of extrapolation onto destination plane\n";
      }
      An = A[0]*SU[0] + A[1]*SU[1] + A[2]*SU[2];
      An = (fabs(An) > 1.E-7 ? 1./An : 0); // 1/A_normal
      double norm;

      M8x8& jacPtr = *jacobianT;
      for(int i = 0; i<8; ++i) {
        norm = (jacPtr(i,0)*SU[0] + jacPtr(i,1)*SU[1] + jacPtr(i,2)*SU[2]) * An;  // dR_normal / A_normal
        jacPtr(i,0) -= norm*A [0];   jacPtr(i,1) -= norm*A [1];   jacPtr(i,2) -= norm*A [2];
        jacPtr(i,3) -= norm*SA[0];   jacPtr(i,4) -= norm*SA[1];   jacPtr(i,5) -= norm*SA[2];
      }
      checkJacProj = true;

      if (debugLvl_ > 0) {
        //std::cout << "  Jacobian^T of extrapolation after Projection:\n";
        //RKTools::printDim(*jacobianT, 7,7);
      }

      for (int iRow = 0; iRow < 3; ++iRow) {
        for (int iCol = 0; iCol < 3; ++iCol) {
          noiseProjection(iRow, iCol)     = (iRow == iCol) - An * SU[iCol] * A[iRow];
          noiseProjection(iRow + 3, iCol) =                - An * SU[iCol] * SA[iRow];
        }
      }

      // noiseProjection will look like this:
      // x x x 0 0 0 0
      // x x x 0 0 0 0
      // x x x 0 0 0 0
      // x x x 1 0 0 0
      // x x x 0 1 0 0
      // x x x 0 0 1 0
      // 0 0 0 0 0 0 1
    }

  } // end of linear extrapolation to surface

  return(true);

}


double RKTrackRepTime::estimateStep(const M1x8& stateGlobal,
                                    const M1x4& SU,
                                    const DetPlane& plane,
                                    const double& charge,
                                    double& relMomLoss,
                                    StepLimits& limits,
                                    MaterialProperties& mat) const
{
  if (useCache_) {
    if (cachePos_ >= RKSteps_.size()) {
      useCache_ = false;
    }
    else {
      if (RKSteps_.at(cachePos_).limits_.getLowestLimit().first == stp_plane) {
        // we need to step exactly to the plane, so don't use the cache!
        useCache_ = false;
        RKSteps_.erase(RKSteps_.begin() + cachePos_, RKSteps_.end());
      }
      else {
        if (debugLvl_ > 0) {
          std::cout << " RKTrackRepTime::estimateStep: use stepSize " << cachePos_ << " from cache: " << RKSteps_.at(cachePos_).matStep_.stepSize_ << "\n";
        }
        //for(int n = 0; n < 1*7; ++n) RKSteps_[cachePos_].state7_[n] = state7[n];
        ++RKStepsFXStop_;
        limits = RKSteps_.at(cachePos_).limits_;
        mat = RKSteps_.at(cachePos_).matStep_.materialProperties_;
        return RKSteps_.at(cachePos_++).matStep_.stepSize_;
      }
    }
  }

  limits.setLimit(stp_sMax, 25.); // max. step allowed [cm]

  if (debugLvl_ > 0) {
    std::cout << " RKTrackRepTime::estimateStep \n";
    std::cout << "  position:  "; TVector3(stateGlobal[0], stateGlobal[1], stateGlobal[2]).Print();
    std::cout << "  direction: "; TVector3(stateGlobal[3], stateGlobal[4], stateGlobal[5]).Print();
  }

  // calculate SL distance to surface
  double Dist ( SU[3] - (stateGlobal[0]*SU[0] +
                         stateGlobal[1]*SU[1] +
                         stateGlobal[2]*SU[2])  );  // Distance between start coordinates and surface
  double An ( stateGlobal[3]*SU[0] +
              stateGlobal[4]*SU[1] +
              stateGlobal[5]*SU[2]   );              // An = dir * N;  component of dir normal to surface

  double SLDist; // signed
  if (fabs(An) > 1.E-10)
    SLDist = Dist/An;
  else {
    SLDist = Dist*1.E10;
    if (An<0) SLDist *= -1.;
  }

  limits.setLimit(stp_plane, SLDist);
  limits.setStepSign(SLDist);

  if (debugLvl_ > 0) {
    std::cout << "  Distance to plane: " << Dist << "\n";
    std::cout << "  SL distance to plane: " << SLDist << "\n";
    if (limits.getStepSign()>0) 
      std::cout << "  Direction is  pointing towards surface.\n";
    else  
      std::cout << "  Direction is pointing away from surface.\n";
  }
  // DONE calculate SL distance to surface

  //
  // Limit according to curvature and magnetic field inhomogenities
  // and improve stepsize estimation to reach plane
  //
  double fieldCurvLimit( limits.getLowestLimitSignedVal() ); // signed
  double remainingDist = 9e99;
  double stepTaken = 9e99;

  MaterialEffects::getInstance()->initTrack(stateGlobal[0] + 0.1*MINSTEP*copysign(stateGlobal[3], SLDist),
                                            stateGlobal[1] + 0.1*MINSTEP*copysign(stateGlobal[4], SLDist),
                                            stateGlobal[2] + 0.1*MINSTEP*copysign(stateGlobal[5], SLDist),
                                            copysign(stateGlobal[3], SLDist),
                                            copysign(stateGlobal[4], SLDist),
                                            copysign(stateGlobal[5], SLDist));
  MaterialEffects::getInstance()->getMaterialProperties(mat);
  double slDist = MaterialEffects::getInstance()->findNextBoundaryStraightLine(fieldCurvLimit);

  // Limit step to not look ahead too far.  Mainly prevents us from
  // extrapolating long distances even though we are in thin sensors.
  fieldCurvLimit = std::min(fieldCurvLimit, 2*slDist);

  RKTrackRepTime::propagator extrap(this, stateGlobal, mat);

  static const unsigned int maxNumIt = 10;
  unsigned int counter(0);
  while (fabs(fieldCurvLimit) > MINSTEP) {

    if(++counter > maxNumIt){
      // if max iterations are reached, take a safe value
      // (in previous iteration, fieldCurvLimit has been not more than doubled)
      // and break.
      fieldCurvLimit *= 0.5;
      break;
    }

    double posAfter[3];
    double dirAfter[3];
    double q ( extrap.extrapolateBy(fieldCurvLimit, posAfter, dirAfter) );

    if (debugLvl_ > 0) {
      std::cout << "  maxStepArg = " << fieldCurvLimit << "; q = " << q  << " \n";
    }

    // remember steps and resulting SL distances to plane for stepsize improvement
    // calculate distance to surface
    Dist = SU[3] - (posAfter[0] * SU[0] +
                    posAfter[1] * SU[1] +
                    posAfter[2] * SU[2]); // Distance between position and surface

    An = dirAfter[0] * SU[0] +
         dirAfter[1] * SU[1] +
         dirAfter[2] * SU[2];    // An = dir * N;  component of dir normal to surface

    if (fabs(Dist/An) < fabs(remainingDist)) {
      remainingDist = Dist/An;
      stepTaken = fieldCurvLimit;
    }

    // resize limit according to q never grow step size more than
    // two-fold to avoid infinite grow-shrink loops with strongly
    // inhomogeneous fields.
    if (q>2) {
      fieldCurvLimit *= 2;
      break;
    }

    fieldCurvLimit *= q * 0.95;

    if (fabs(q-1) < 0.25 || // good enough!
        fabs(fieldCurvLimit) > limits.getLowestLimitVal()) // other limits are lower!
      break;
  }
  if (fabs(fieldCurvLimit) < MINSTEP)
    limits.setLimit(stp_fieldCurv, MINSTEP);
  else
    limits.setLimit(stp_fieldCurv, fieldCurvLimit);

  double stepToPlane(limits.getLimitSigned(stp_plane));
  if (fabs(remainingDist) < 8.E99) {
    stepToPlane = stepTaken + remainingDist;
  }
  limits.setLimit(stp_plane, stepToPlane);


  //
  // Select direction
  //
  // auto select
  if (propDir_ == 0 || !plane.isFinite()){
    if (debugLvl_ > 0) {
      std::cout << "  auto select direction";
      if (!plane.isFinite()) std::cout << ", plane is not finite";
      std::cout << ".\n";
    }
  }
  // see if straight line approximation is ok
  else if ( limits.getLimit(stp_plane) < 0.2*limits.getLimit(stp_fieldCurv) ){
    if (debugLvl_ > 0) {
      std::cout << "  straight line approximation is fine.\n";
    }

    // if direction is pointing to active part of surface
    if( plane.isInActive(stateGlobal[0], stateGlobal[1], stateGlobal[2],  stateGlobal[3], stateGlobal[4], stateGlobal[5]) ) {
      if (debugLvl_ > 0) {
        std::cout << "  direction is pointing to active part of surface. \n";
      }
    }
    // if we are near the plane, but not pointing to the active area, make a big step!
    else {
      limits.removeLimit(stp_plane);
      limits.setStepSign(propDir_);
      if (debugLvl_ > 0) {
        std::cout << "  we are near the plane, but not pointing to the active area. make a big step! \n";
      }
    }
  }
  // propDir_ is set and we are not pointing to an active part of a plane -> propDir_ decides!
  else {
    if (limits.getStepSign() * propDir_ < 0){
      limits.removeLimit(stp_plane);
      limits.setStepSign(propDir_);
      if (debugLvl_ > 0) {
        std::cout << "  invert Step according to propDir_ and make a big step. \n";
      }
    }
  }


  // call stepper and reduce stepsize if step not too small
  static const RKStep defaultRKStep;
  RKSteps_.push_back( defaultRKStep );
  std::vector<RKStep>::iterator lastStep = RKSteps_.end() - 1;
  std::copy(stateGlobal.begin(), stateGlobal.end(), lastStep->state_.begin());
  ++RKStepsFXStop_;

  if(limits.getLowestLimitVal() > MINSTEP){ // only call stepper if step estimation big enough
    MaterialEffects::getInstance()->stepper(extrap,
                                            charge/stateGlobal[6], // |p|
                                            relMomLoss,
                                            pdgCode_,
                                            lastStep->matStep_.materialProperties_,
                                            limits);
  } else { //assume material has not changed
    if  (RKSteps_.size()>1) {
      lastStep->matStep_.materialProperties_ = (lastStep - 1)->matStep_.materialProperties_;
    }
  }

  if (debugLvl_ > 0) {
    std::cout << "   final limits:\n";
    limits.Print();
  }

  // With this RK implementation, mom loss is part of what is considered curvature here.
  limits.removeLimit(stp_momLoss);

  double finalStep = limits.getLowestLimitSignedVal();

  lastStep->matStep_.stepSize_ = finalStep;
  lastStep->limits_ = limits;

  if (debugLvl_ > 0) {
    std::cout << "  --> Step to be used: " << finalStep << "\n";
    //limits.Print();
  }

  return finalStep;
}


TVector3 RKTrackRepTime::pocaOnLine(const TVector3& linePoint, const TVector3& lineDirection, const TVector3& point) const {

  TVector3 retVal(lineDirection);

  double t = 1./(retVal.Mag2()) * ((point*retVal) - (linePoint*retVal));
  retVal *= t;
  retVal += linePoint;
  return retVal; // = linePoint + t*lineDirection

}


double RKTrackRepTime::Extrap(const DetPlane& startPlane,
                          const DetPlane& destPlane,
                          double charge,
                          double mass,
                          bool& isAtBoundary,
                          M1x8& stateGlobal,
                          double& flightTime,
                          bool fillExtrapSteps,
                          TMatrixDSym* cov, // 5D
                          bool onlyOneStep,
                          bool stopAtBoundary,
                          double maxStep) const
{

  static const unsigned int maxNumIt(500);
  unsigned int numIt(0);

  double coveredDistance(0.);

  const TVector3& W(destPlane.getNormal());
  M1x4 SU = {{W.X(), W.Y(), W.Z(), destPlane.distance(0., 0., 0.)}};

  // make SU vector point away from origin
  if (W*destPlane.getO() < 0) {
    SU[0] *= -1;
    SU[1] *= -1;
    SU[2] *= -1;
  }


  M1x8 startStateGlobal = stateGlobal;
  while(true){

    if (debugLvl_ > 0) {
      std::cout << "\n============ RKTrackRepTime::Extrap loop nr. " << numIt << " ============\n";
      std::cout << "Start plane: "; startPlane.Print();
      std::cout << "fillExtrapSteps " << fillExtrapSteps << "\n";
    }

    if(++numIt > maxNumIt){
      Exception exc("RKTrackRepTime::Extrap ==> maximum number of iterations exceeded",__LINE__,__FILE__);
      exc.setFatal();
      throw exc;
    }

    // initialize jacobianT with unit matrix
    M8x8 J_MMT;
    for(int i = 0; i < 8*8; ++i) J_MMT[i] = 0;
    for(int i=0; i<8; ++i) J_MMT(i,i) = 1.;

    isAtBoundary = false;

    // propagation
    bool checkJacProj = false;
    limits_.reset();
    limits_.setLimit(stp_sMaxArg, maxStep-fabs(coveredDistance));

    double pStart = fabs(charge / stateGlobal[6]);
    if( ! RKutta(SU, destPlane, charge, mass, stateGlobal, fillExtrapSteps ? &J_MMT : 0,
		 coveredDistance, flightTime, checkJacProj, noiseProjection_,
		 limits_, onlyOneStep) ) {
      Exception exc("RKTrackRepTime::Extrap ==> Runge Kutta propagation failed",__LINE__,__FILE__);
      exc.setFatal();
      throw exc;
    }

    bool atPlane(limits_.getLowestLimit().first == stp_plane);
    if (limits_.getLowestLimit().first == stp_boundary)
      isAtBoundary = true;


    if (debugLvl_ > 0) {
      std::cout<<"RKSteps \n";
      for (std::vector<RKStep>::iterator it = RKSteps_.begin(); it != RKSteps_.end(); ++it){
        std::cout << "stepSize = " << it->matStep_.stepSize_ << "\t";
        it->matStep_.materialProperties_.Print();
      }
      std::cout<<"\n";
    }



    // call MatFX
    M7x7* pNoise = NULL;
    M7x7 noise;
    if(fillExtrapSteps) {
      pNoise = &noise;
      for(int i = 0; i < 7*7; ++i) noise[i] = 0; // set noiseArray_ to 0
    }

    double momLoss = 0;
    unsigned int nPoints(RKStepsFXStop_ - RKStepsFXStart_);
    if (nPoints>0){
      // momLoss has a sign - negative loss means momentum gain
      momLoss += MaterialEffects::getInstance()->effects(RKSteps_,
                                                         RKStepsFXStart_,
                                                         RKStepsFXStop_,
                                                         pStart, // momentum
                                                         pdgCode_,
                                                         pNoise);

      RKStepsFXStart_ = RKStepsFXStop_;

      if (debugLvl_ > 0) {
        std::cout << "momLoss: " << momLoss << " GeV; relative: " << momLoss/fabs(charge/stateGlobal[6])
            << "; coveredDistance = " << coveredDistance << "\n";
        if (debugLvl_ > 1 && pNoise) {
          std::cout << "7D noise: \n";
          noise.print();
        }
      }
    } // finished MatFX


    // fill ExtrapSteps_
    if (fillExtrapSteps) {
      static const ExtrapStep defaultExtrapStep;
      ExtrapSteps_.push_back(defaultExtrapStep);
      std::vector<ExtrapStep>::iterator lastStep = ExtrapSteps_.end() - 1;

      // Store Jacobian of this step for final calculation.
      lastStep->jac_ = J_MMT;

      if( checkJacProj == true ){
        //project the noise onto the destPlane
        RKTools::Np_N_NpT(noiseProjection_, noise);

        if (debugLvl_ > 1) {
          std::cout << "7D noise projected onto plane: \n";
          noise.print();
        }
      }

      // Store this step's noise for final calculation.
      lastStep->noise_ = noise;

      if (debugLvl_ > 2) {
        std::cout<<"ExtrapSteps \n";
        for (std::vector<ExtrapStep>::iterator it = ExtrapSteps_.begin(); it != ExtrapSteps_.end(); ++it){
          std::cout << "8D Jacobian: "; it->jac_.print();
          std::cout << "7D noise:    "; it->noise_.print();
        }
        std::cout<<"\n";
      }
    }



    // check if at boundary
    if (stopAtBoundary) {
      if (debugLvl_ > 0) {
        std::cout << "stopAtBoundary -> break; \n ";
      }
      break;
    }

    if (onlyOneStep) {
      if (debugLvl_ > 0) {
        std::cout << "onlyOneStep -> break; \n ";
      }
      break;
    }

    //break if we arrived at destPlane
    if(atPlane) {
      if (debugLvl_ > 0) {
        std::cout << "arrived at destPlane with a distance of  " << destPlane.distance(stateGlobal[0], stateGlobal[1], stateGlobal[2]) << " cm left. ";
        if (destPlane.isInActive(stateGlobal[0], stateGlobal[1], stateGlobal[2],  stateGlobal[3], stateGlobal[4], stateGlobal[5]))
          std::cout << "In active area of destPlane. \n";
        else
          std::cout << "NOT in active area of plane. \n";

        std::cout << "  position:  "; TVector3(stateGlobal[0], stateGlobal[1], stateGlobal[2]).Print();
        std::cout << "  direction: "; TVector3(stateGlobal[3], stateGlobal[4], stateGlobal[5]).Print();
      }
      break;
    }

  }

  if (fillExtrapSteps) {
    // propagate cov and add noise
    calcForwardJacobianAndNoise(startStateGlobal, startPlane, stateGlobal, destPlane);

    if (cov != NULL) {
      cov->Similarity(fJacobian_);
      *cov += fNoise_;
    }

    if (debugLvl_ > 0) {
      if (cov != NULL) {
        std::cout << "final covariance matrix after Extrap: "; cov->Print();
      }
    }
  }

  return coveredDistance;
}


void RKTrackRepTime::checkCache(const StateOnPlane& state, const SharedPlanePtr* plane) const {

  if (state.getRep() != this){
    Exception exc("RKTrackRepTime::checkCache ==> state is defined wrt. another TrackRep",__LINE__,__FILE__);
    exc.setFatal();
    throw exc;
  }

  if (dynamic_cast<const MeasurementOnPlane*>(&state) != NULL) {
    Exception exc("RKTrackRepTime::checkCache - cannot extrapolate MeasurementOnPlane",__LINE__,__FILE__);
    exc.setFatal();
    throw exc;
  }

  cachePos_ = 0;
  RKStepsFXStart_ = 0;
  RKStepsFXStop_ = 0;
  ExtrapSteps_.clear();
  initArrays();


  if (plane &&
      lastStartState_.getPlane() &&
      lastEndState_.getPlane() &&
      state.getPlane() == lastStartState_.getPlane() &&
      state.getState() == lastStartState_.getState() &&
      (*plane)->distance(getPos(lastEndState_)) <= MINSTEP) {
    useCache_ = true;

    // clean up cache. Only use steps with same sign.
    double firstStep(0);
    for (unsigned int i=0; i<RKSteps_.size(); ++i) {
      if (i == 0) {
        firstStep = RKSteps_.at(0).matStep_.stepSize_;
        continue;
      }
      if (RKSteps_.at(i).matStep_.stepSize_ * firstStep < 0) {
        if (RKSteps_.at(i-1).matStep_.materialProperties_ == RKSteps_.at(i).matStep_.materialProperties_) {
          RKSteps_.at(i-1).matStep_.stepSize_ += RKSteps_.at(i).matStep_.stepSize_;
        }
        RKSteps_.erase(RKSteps_.begin()+i, RKSteps_.end());
      }
    }

    if (debugLvl_ > 0) {
        std::cout << "RKTrackRepTime::checkCache: use cached material and step values.\n";
    }
  }
  else {

    if (debugLvl_ > 0) {
      std::cout << "RKTrackRepTime::checkCache: can NOT use cached material and step values.\n";

      if (plane != NULL) {
        if (state.getPlane() != lastStartState_.getPlane()) {
          std::cout << "state.getPlane() != lastStartState_.getPlane()\n";
        }
        else {
          if (! (state.getState() == lastStartState_.getState())) {
            std::cout << "state.getState() != lastStartState_.getState()\n";
          }
          else if (lastEndState_.getPlane().get() != NULL) {
            std::cout << "distance " << (*plane)->distance(getPos(lastEndState_)) << "\n";
          }
        }
      }
    }

    useCache_ = false;
    RKSteps_.clear();

    lastStartState_.setStatePlane(state.getState(), state.getPlane());
  }
}


double RKTrackRepTime::momMag(const M1x8& stateGlobal) const {
  // FIXME given this interface this function cannot work for charge =/= +-1
  return fabs(1/stateGlobal[6]);
}


bool RKTrackRepTime::isSameType(const AbsTrackRep* other) {
  if (dynamic_cast<const RKTrackRepTime*>(other) == NULL)
    return false;

  return true;
}


bool RKTrackRepTime::isSame(const AbsTrackRep* other) {
  if (getPDG() != other->getPDG())
    return false;

  return isSameType(other);
}


void RKTrackRepTime::Streamer(TBuffer &R__b)
{
   // Stream an object of class genfit::RKTrackRepTime.

   //This works around a msvc bug and should be harmless on other platforms
   typedef ::genfit::RKTrackRepTime thisClass;
   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      ::genfit::AbsTrackRep::Streamer(R__b);
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      R__b.CheckByteCount(R__s, R__c, thisClass::IsA());
      lastStartState_.setRep(this);
      lastEndState_.setRep(this);
   } else {
      ::genfit::AbsTrackRep::Streamer(R__b);
      R__c = R__b.WriteVersion(thisClass::IsA(), kTRUE);
      R__b.SetByteCount(R__c, kTRUE);
   }
}

} /* End of namespace genfit */
