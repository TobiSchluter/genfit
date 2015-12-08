#include <framework/logging/Logger.h>
#include <tracking/modules/genfitter/Geant4MaterialInterface.h>
#include <geometry/GeometryManager.h>

#include "genfit/Exception.h"
#include "genfit/RKTrackRep.h"

#include <assert.h>
#include <math.h>

#include "G4ThreeVector.hh"
#include "G4Navigator.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Material.hh"
#include "G4TouchableHistory.hh"

static const bool debug = false;
//static const bool debug = true;

using namespace Belle2;

namespace Belle2 {
  class G4SafeNavigator {
    // Guards against leaving the physical volume.
    //
    // Not inheriting from G4Navigator because CheckNextStep is not
    // virtual.
  public:
    G4SafeNavigator() = default;
    ~G4SafeNavigator() = default;

    G4VPhysicalVolume* GetWorldVolume() const { return nav_.GetWorldVolume(); }
    void SetWorldVolume(G4VPhysicalVolume* pWorld)
    {
      nav_.SetWorldVolume(pWorld);
      worldsolid_ = pWorld->GetLogicalVolume()->GetSolid();
    }

    G4VPhysicalVolume* LocateGlobalPointAndSetup(const G4ThreeVector& point,
                                                 const G4ThreeVector* direction = 0,
                                                 const G4bool pRelativeSearch = true,
                                                 const G4bool ignoreDirection = true);

    G4double CheckNextStep(const G4ThreeVector& pGlobalPoint,
                           const G4ThreeVector& pDirection,
                           const G4double pCurrentProposedStepLength,
                           G4double& pNewSafety);

    G4VPhysicalVolume* ResetHierarchyAndLocate(const G4ThreeVector& point,
                                               const G4ThreeVector& direction,
                                               const G4TouchableHistory& h);
    G4TouchableHistory* CreateTouchableHistory() const
    {
      return nav_.CreateTouchableHistory();
    }

    G4double straightLineDist(const G4double sMax)
    {
      double safety;
      return CheckNextStep(lastpoint_, lastdir_, fabs(sMax), safety);
    }

    G4Navigator& getNav() { return nav_; }
  private:
    G4ThreeVector lastpoint_;
    G4ThreeVector lastdir_;
    G4VPhysicalVolume* lastvolume_{0};
    G4Navigator nav_;
    const G4VSolid* worldsolid_{0};
  };
}

G4VPhysicalVolume* G4SafeNavigator::LocateGlobalPointAndSetup(const G4ThreeVector& point,
    const G4ThreeVector* direction,
    const G4bool pRelativeSearch,
    const G4bool ignoreDirection)
{
  if (point == lastpoint_ && lastvolume_) {
    return lastvolume_;
  }
  //B2INFO("###  init: " << point);
  G4VPhysicalVolume* volume = nav_.LocateGlobalPointAndSetup(point, direction, pRelativeSearch, ignoreDirection);
  if (!volume) {
    volume = nav_.GetWorldVolume();
  }
  // remember last point to speed up setup if possible
  lastpoint_ = point;
  if (direction) lastdir_ = *direction;
  lastvolume_ = volume;
  return volume;
}

G4VPhysicalVolume* G4SafeNavigator::ResetHierarchyAndLocate(const G4ThreeVector& point,
                                                            const G4ThreeVector& direction,
                                                            const G4TouchableHistory& h)
{
  if (point == lastpoint_ && lastvolume_) {
    return lastvolume_;
  }
  //B2INFO("### reset: " << point);
  G4VPhysicalVolume* volume = nav_.ResetHierarchyAndLocate(point, direction, h);
  if (!volume) {
    volume = nav_.GetWorldVolume();
  }
  // remember last point to speed up setup if possible
  lastpoint_ = point;
  lastdir_ = direction;
  lastvolume_ = volume;
  return volume;
}

G4double G4SafeNavigator::CheckNextStep(const G4ThreeVector& point,
                                        const G4ThreeVector& direction,
                                        const G4double pCurrentProposedStepLength,
                                        G4double& pNewSafety)
{
  //make sure we're inside the world volume
  if (worldsolid_->Inside(point) == kOutside) {
    pNewSafety = worldsolid_->DistanceToIn(point);
    return worldsolid_->DistanceToIn(point, direction);
  }
  return nav_.CheckNextStep(point, direction, pCurrentProposedStepLength, pNewSafety);
}


Geant4MaterialInterface::Geant4MaterialInterface()
  : nav_(new G4SafeNavigator()), currentVolume_(0)
{
  G4VPhysicalVolume* world = geometry::GeometryManager::getInstance().getTopVolume();
  nav_->SetWorldVolume(world);
}

Geant4MaterialInterface::~Geant4MaterialInterface()
{
}


bool
Geant4MaterialInterface::initTrack(double posX, double posY, double posZ,
                                   double dirX, double dirY, double dirZ)
{
  G4ThreeVector pos(posX * CLHEP::cm, posY * CLHEP::cm, posZ * CLHEP::cm);
  G4ThreeVector dir(dirX, dirY, dirZ);
  const G4VPhysicalVolume* newVolume = nav_->LocateGlobalPointAndSetup(pos, &dir);
  bool volChanged = newVolume != currentVolume_;
  currentVolume_ = newVolume;

  if (volChanged && (nav_->getNav().EnteredDaughterVolume() || nav_->getNav().ExitedMotherVolume())) {
    bool valid = false;
    lastNormal_ = G4ThreeVector(dirX, dirY, dirZ);//nav_->getNav().GetGlobalExitNormal(pos, &valid);
  } else {
    lastNormal_ = G4ThreeVector(dirX, dirY, dirZ);
  }
  return volChanged;
}


void
Geant4MaterialInterface::getMaterialParameters(double& density,
                                               double& Z,
                                               double& A,
                                               double& radiationLength,
                                               double& mEE)
{
  assert(currentVolume_);

  const G4Material* mat = currentVolume_->GetLogicalVolume()->GetMaterial();

  if (mat->GetNumberOfElements() == 1) {
    Z = mat->GetZ();
    A = mat->GetA();
  } else {
    // Calculate weight-averaged A, Z
    A = Z = 0;
    for (unsigned i = 0; i < mat->GetNumberOfElements(); ++i) {
      const G4Element* element = (*mat->GetElementVector())[i];
      Z += element->GetZ() * mat->GetFractionVector()[i];
      A += element->GetA() * mat->GetFractionVector()[i];
    }
  }

  density = mat->GetDensity() / CLHEP::g * CLHEP::cm3;
  // Z has correct units
  A *= CLHEP::mole / CLHEP::g;
  radiationLength = mat->GetRadlen() / CLHEP::cm;
  mEE = mat->GetIonisation()->GetMeanExcitationEnergy() / CLHEP::eV;
}


void
Geant4MaterialInterface::getMaterialParameters(genfit::MaterialProperties& parameters)
{
  double density, Z, A, radLen, mEE;
  this->Geant4MaterialInterface::getMaterialParameters(density, Z, A, radLen, mEE);
  parameters.setMaterialProperties(density, Z, A, radLen, mEE);
}


double
Geant4MaterialInterface::findNextBoundary(const genfit::AbsTrackRep::internalExtrapolator& extrap,
                                          double sMax) // signed
{
  //bool debug =  true;
  // This assumes that sMax is small enough to take only a single RK
  // step.  This restriction comes about because RKPropagate only
  // takes one step.
  const double delta(1.E-2); // cm, distance limit beneath which straight-line steps are taken.
  const double epsilon(1.E-1); // cm, allowed upper bound on arch deviation from straight line

  double posOrig[3];
  double dirOrig[3];
  extrap.getInitialState(posOrig, dirOrig);

  double posLast[3] = {posOrig[0], posOrig[1], posOrig[2]};

  int stepSign(sMax < 0 ? -1 : 1);

  G4ThreeVector pointOld(posOrig[0] * CLHEP::cm,
                         posOrig[1] * CLHEP::cm,
                         posOrig[2] * CLHEP::cm);
  G4ThreeVector dirOld(stepSign * dirOrig[0],
                       stepSign * dirOrig[1],
                       stepSign * dirOrig[2]);

  double s = 0;  // trajectory length to boundary

  const unsigned maxIt = 300;
  unsigned it = 0;

  // Initialize the geometry to the current location (set by caller).
  double safety;
  double slDist = nav_->CheckNextStep(pointOld, dirOld, fabs(sMax) * CLHEP::cm, safety);
  if (slDist == kInfinity)
    slDist = fabs(sMax);
  else
    slDist /= CLHEP::cm;
  safety /= CLHEP::cm;

  // No boundary in sight?
  if (safety > fabs(sMax)) {
    if (debug)
      std::cout << "   next boundary is farther away than sMax \n";
    return stepSign * safety; // sMax
  }

  // Are we at the boundary?
  if (slDist < delta) {
    if (debug) {
      std::cout << "   very close to the boundary pre-loop -> return @ it " << it
                << " stepSign*slDist = "
                << stepSign << "*" << slDist << "\n";
    }
    nav_->getNav().ComputeStep(pointOld, dirOld, fabs(sMax) * CLHEP::cm, safety);
    if (nav_->getNav().EnteredDaughterVolume() || nav_->getNav().ExitedMotherVolume()) {
      bool valid = false;
      G4ThreeVector normal = nav_->getNav().GetGlobalExitNormal(pointOld, &valid);
      if (valid) lastNormal_ = normal;
    }
    return stepSign * slDist;
  }
  double step = slDist;

  while (1) {
    if (++it > maxIt) {
      genfit::Exception exc("Geant4MaterialInterface::findNextBoundary ==> maximum number of iterations exceeded", __LINE__, __FILE__);
      exc.setFatal();
      throw exc;
    }

    if (step < delta) {
      if (debug) {
        std::cout << "step < delta ";
        std::cout << nav_->getNav().EnteredDaughterVolume() << " " << nav_->getNav().ExitedMotherVolume() << std::endl;
      }

      nav_->getNav().ComputeStep(pointOld, dirOld, fabs(sMax) * CLHEP::cm, safety);
      if (nav_->getNav().EnteredDaughterVolume() || nav_->getNav().ExitedMotherVolume()) {
        bool valid = false;
        G4ThreeVector normal = nav_->getNav().GetGlobalExitNormal(pointOld, &valid);
        if (valid) lastNormal_ = normal;
      }
      return stepSign * (s + step);
    }
    // We have to find whether there's any boundary on our path.

    // Follow curved arch, then see if we may have missed a boundary.
    // Always propagate complete way from original start to avoid
    // inconsistent extrapolations.  This is always a single RK step.
    double posNew[3];
    double dirNew[3];
    extrap.extrapolateBy(stepSign * (s + step), posNew, dirNew);

    G4ThreeVector pos(posNew[0] * CLHEP::cm, posNew[1] * CLHEP::cm, posNew[2] * CLHEP::cm);
    G4ThreeVector dir(stepSign * dirNew[0], stepSign * dirNew[1], stepSign * dirNew[2]);

    // Straight line distance between extrapolation finish and
    // the end of the previously determined safe segment.
    double dist2 = (pow(posNew[0] - posLast[0], 2)
                    + pow(posNew[1] - posLast[1], 2)
                    + pow(posNew[2] - posLast[2], 2));

    // If we moved less than safety, the volume cannot possibly have
    // changed, so we skip further checks.
    if (dist2 > safety * safety) {

      // Do we need to try again with a shorter step?  There are two
      // possible reasons:
      //
      // 1. Too much curvature?

      // Maximal lateral deviation² of the curved path from the
      // straight line connecting beginning and end.
      double maxDeviation2 = 0.25 * (step * step - dist2);
      if (maxDeviation2 > epsilon * epsilon) {
        // Need to take a shorter step to reliably estimate material,
        // but only if we didn't move by safety.

        // Take a shorter step, but never shorter than safety.
        step = safety + 0.5 * (step - safety);

        continue;
      }

      // 2. Volume changed?
      //
      // Where are we after the step?
      std::unique_ptr<G4TouchableHistory> hist(nav_->CreateTouchableHistory());
      G4VPhysicalVolume* newVolume = nav_->LocateGlobalPointAndSetup(pos, &dir);

      if (newVolume != currentVolume_) {

        // Volume changed during the extrapolation.

        // Extrapolation may not take the exact step length we asked
        // for, so it can happen that a requested step < safety takes
        // us across the boundary.  This is then the best estimate we
        // can get of the distance to the boundary with the stepper.
        if (step <= safety) {
          if (debug) {
            std::cout << "step <= safety";
            std::cout << nav_->getNav().EnteredDaughterVolume() << " " << nav_->getNav().ExitedMotherVolume() << std::endl;
          }

          lastNormal_ = dir;
          return stepSign * (s + step);
        }

        // Move back to last good point, but looking in the actual
        // direction of the step.
        G4ThreeVector dirCloser(pos - pointOld);
        dirCloser.setMag(1.);
        nav_->ResetHierarchyAndLocate(pointOld, dirCloser, *hist);

        // Look along the secant of the actual trajectory instead of
        // the original direction.  There should be a crossing within
        // distance step.
        double secantDist = nav_->CheckNextStep(pointOld, dirCloser,
                                                step * CLHEP::cm, safety) / CLHEP::cm;
        safety /= CLHEP::cm;
        if (secantDist >= step) {
          // Cannot be.  Just take a shorter step, and hope that this
          // works.
          slDist = secantDist;
          step = std::max(0.9 * step, safety);
        } else {
          slDist = step = std::max(secantDist, safety);
        }

        // Are we at the boundary?
        if (slDist < delta) {
          if (debug) {
            std::cout << "   very close to the boundary in loop -> return @ it " << it
                      << " stepSign*(s + slDist) = "
                      << stepSign << "*(" << s + slDist << ")\n";
            std::cout << nav_->getNav().EnteredDaughterVolume() << " " << nav_->getNav().ExitedMotherVolume() << std::endl;

          }
          nav_->getNav().ComputeStep(pointOld, dirCloser, step * CLHEP::cm, safety);
          if (nav_->getNav().EnteredDaughterVolume() || nav_->getNav().ExitedMotherVolume()) {
            bool valid = false;
            G4ThreeVector normal = nav_->getNav().GetGlobalExitNormal(pointOld, &valid);
            if (valid) lastNormal_ = normal;
          }

          return stepSign * (s + slDist);
        }

        continue;
      }
    }

    // We're in the new place, the step was safe, advance.
    s += step;
    std::copy(posNew, posNew + 3, posLast);
    pointOld = pos;
    nav_->LocateGlobalPointAndSetup(pos, &dir);
    slDist = nav_->CheckNextStep(pos, dir,
                                 (fabs(sMax) - s) * CLHEP::cm, safety);
    if (slDist == kInfinity)
      slDist = fabs(sMax) - s;
    else
      slDist /= CLHEP::cm;
    safety /= CLHEP::cm;
    step = slDist;

    // No boundary in sight?
    if (s + safety > fabs(sMax)) {
      if (debug)
        std::cout << "   next boundary is farther away than sMax \n";
      lastNormal_ = dir;
      return stepSign * (s + safety); // sMax
    }

    // Are we at the boundary?
    if (slDist < delta) {
      if (debug) {
        std::cout << "   very close to the boundary near end -> return @ it " << it
                  << " stepSign*(s + slDist) = "
                  << stepSign << "*(" << s + slDist << ")\n";
        std::cout << nav_->getNav().EnteredDaughterVolume() << " " << nav_->getNav().ExitedMotherVolume() << std::endl;
      }
      nav_->getNav().ComputeStep(pos, dir, (fabs(sMax) - s) * CLHEP::cm, safety);
      if (nav_->getNav().EnteredDaughterVolume() || nav_->getNav().ExitedMotherVolume()) {
        bool valid = false;
        G4ThreeVector normal = nav_->getNav().GetGlobalExitNormal(pointOld, &valid);
        if (valid) lastNormal_ = normal;
      }

      return stepSign * s + slDist;
    }
  }
}


double Geant4MaterialInterface::findNextBoundaryStraightLine(double sMax)
{
  return nav_->straightLineDist(fabs(sMax));
}

#if 0
double findNextBoundary(const AbsTrackRep::internalExtrapolator& extrap,
                        double sMax,
                        double normal[3],
                        MaterialProperties& materialBefore,
                        MaterialProperties& materialAfter)
{
  double posOrig[3];
  double dirOrig[3];
  extrap.getInitialState(posOrig, dirOrig);

  initTrack(posOrig[0], posOrig[1], posOrig[2],
            dirOrig[0], dirOrig[1], dirOrig[2]);
  getMaterialParameters(materialBefore);

  double posLast[3] = {posOrig[0], posOrig[1], posOrig[2]};

  int stepSign(sMax < 0 ? -1 : 1);

  G4ThreeVector pointOld(posOrig[0] * CLHEP::cm,
                         posOrig[1] * CLHEP::cm,
                         posOrig[2] * CLHEP::cm);
  G4ThreeVector dirOld(stepSign * dirOrig[0],
                       stepSign * dirOrig[1],
                       stepSign * dirOrig[2]);

  double s = 0;  // trajectory length to boundary

  const unsigned maxIt = 300;
  unsigned it = 0;

  // Initialize the geometry to the current location (set by caller).
  double safety;
  double slDist = nav_->CheckNextStep(pointOld, dirOld, fabs(sMax) * CLHEP::cm, safety);
  if (slDist == kInfinity)
    slDist = fabs(sMax);
  else
    slDist /= CLHEP::cm;
  safety /= CLHEP::cm;

  // No boundary in sight?
  if (safety > fabs(sMax)) {
    if (debug)
      std::cout << "   next boundary is farther away than sMax \n";
    materialAfter = materialBefore;
    std::copy(dirOrig, dirOrig + 3, normal);
    return stepSign * safety; // sMax
  }

  // Are we at the boundary?
  if (slDist < delta) {
    nav_->getN
    if (debug) {
      std::cout << "   very close to the boundary pre-loop -> return @ it " << it
                << " stepSign*slDist = "
                << stepSign << "*" << slDist << "\n";
      std::cout << nav_->getNav().EnteredDaughterVolume() << " " << nav_->getNav().ExitedMotherVolume() << std::endl;
    }
    return stepSign * slDist;
  }
  double step = slDist;

  while (1) {
    if (++it > maxIt) {
      genfit::Exception exc("Geant4MaterialInterface::findNextBoundary ==> maximum number of iterations exceeded", __LINE__, __FILE__);
      exc.setFatal();
      throw exc;
    }

    if (step < delta) {
      if (debug) {
        std::cout << "step < delta ";
        std::cout << nav_->getNav().EnteredDaughterVolume() << " " << nav_->getNav().ExitedMotherVolume() << std::endl;
      }

      return stepSign * (s + step);
    }
    // We have to find whether there's any boundary on our path.

    // Follow curved arch, then see if we may have missed a boundary.
    // Always propagate complete way from original start to avoid
    // inconsistent extrapolations.  This is always a single RK step.
    double posNew[3];
    double dirNew[3];
    extrap.extrapolateBy(stepSign * (s + step), posNew, dirNew);

    G4ThreeVector pos(posNew[0] * CLHEP::cm, posNew[1] * CLHEP::cm, posNew[2] * CLHEP::cm);
    G4ThreeVector dir(stepSign * dirNew[0], stepSign * dirNew[1], stepSign * dirNew[2]);

    // Straight line distance between extrapolation finish and
    // the end of the previously determined safe segment.
    double dist2 = (pow(posNew[0] - posLast[0], 2)
                    + pow(posNew[1] - posLast[1], 2)
                    + pow(posNew[2] - posLast[2], 2));

    // If we moved less than safety, the volume cannot possibly have
    // changed, so we skip further checks.
    if (dist2 > safety * safety) {

      // Do we need to try again with a shorter step?  There are two
      // possible reasons:
      //
      // 1. Too much curvature?

      // Maximal lateral deviation² of the curved path from the
      // straight line connecting beginning and end.
      double maxDeviation2 = 0.25 * (step * step - dist2);
      if (maxDeviation2 > epsilon * epsilon) {
        // Need to take a shorter step to reliably estimate material,
        // but only if we didn't move by safety.

        // Take a shorter step, but never shorter than safety.
        step = safety + 0.5 * (step - safety);

        continue;
      }

      // 2. Volume changed?
      //
      // Where are we after the step?
      std::unique_ptr<G4TouchableHistory> hist(nav_->CreateTouchableHistory());
      G4VPhysicalVolume* newVolume = nav_->LocateGlobalPointAndSetup(pos, &dir);

      if (newVolume != currentVolume_) {

        // Volume changed during the extrapolation.

        // Extrapolation may not take the exact step length we asked
        // for, so it can happen that a requested step < safety takes
        // us across the boundary.  This is then the best estimate we
        // can get of the distance to the boundary with the stepper.
        if (step <= safety) {
          if (debug) {
            std::cout << "step <= safety";
            std::cout << nav_->getNav().EnteredDaughterVolume() << " " << nav_->getNav().ExitedMotherVolume() << std::endl;
          }

          return stepSign * (s + step);
        }

        // Move back to last good point, but looking in the actual
        // direction of the step.
        G4ThreeVector dirCloser(pos - pointOld);
        dirCloser.setMag(1.);
        nav_->ResetHierarchyAndLocate(pointOld, dirCloser, *hist);

        // Look along the secant of the actual trajectory instead of
        // the original direction.  There should be a crossing within
        // distance step.
        double secantDist = nav_->CheckNextStep(pointOld, dirCloser,
                                                step * CLHEP::cm, safety) / CLHEP::cm;
        safety /= CLHEP::cm;
        if (secantDist >= step) {
          // Cannot be.  Just take a shorter step, and hope that this
          // works.
          slDist = secantDist;
          step = std::max(0.9 * step, safety);
        } else {
          slDist = step = std::max(secantDist, safety);
        }

        // Are we at the boundary?
        if (slDist < delta) {
          if (debug) {
            std::cout << "   very close to the boundary in loop -> return @ it " << it
                      << " stepSign*(s + slDist) = "
                      << stepSign << "*(" << s + slDist << ")\n";
            std::cout << nav_->getNav().EnteredDaughterVolume() << " " << nav_->getNav().ExitedMotherVolume() << std::endl;

          }
          return stepSign * (s + slDist);
        }

        continue;
      }
    }

    // We're in the new place, the step was safe, advance.
    s += step;
    std::copy(posNew, posNew + 3, posLast);
    pointOld = pos;
    nav_->LocateGlobalPointAndSetup(pos, &dir);
    slDist = nav_->CheckNextStep(pos, dir,
                                 (fabs(sMax) - s) * CLHEP::cm, safety);
    if (slDist == kInfinity)
      slDist = fabs(sMax) - s;
    else
      slDist /= CLHEP::cm;
    safety /= CLHEP::cm;
    step = slDist;

    // No boundary in sight?
    if (s + safety > fabs(sMax)) {
      if (debug)
        std::cout << "   next boundary is farther away than sMax \n";
      return stepSign * (s + safety); // sMax
    }

    // Are we at the boundary?
    if (slDist < delta) {
      if (debug) {
        std::cout << "   very close to the boundary near end -> return @ it " << it
                  << " stepSign*(s + slDist) = "
                  << stepSign << "*(" << s + slDist << ")\n";
        std::cout << nav_->getNav().EnteredDaughterVolume() << " " << nav_->getNav().ExitedMotherVolume() << std::endl;
      }
      return stepSign * s + slDist;
    }
  }
}
#endif
