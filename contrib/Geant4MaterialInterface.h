/* Copyright 2014, Ludwig-Maximilians-Universität München,
   Authors: Tobias Schlüter

   Provided as part of the Belle II software framework basf2.  Its
   licenses apply.
*/

#pragma once

#include <memory>

#include "genfit/AbsMaterialInterface.h"

#include <G4ThreeVector.hh>

class G4VPhysicalVolume;

namespace Belle2 {


  /**
   * @brief AbsMaterialInterface implementation for use with ROOT's TGeoManager.
   */
  class Geant4MaterialInterface : public genfit::AbsMaterialInterface {

  public:

    Geant4MaterialInterface();
    ~Geant4MaterialInterface();

    /** @brief Initialize the navigator at given position and with given
        direction.  Returns true if the volume changed.
     */
    bool initTrack(double posX, double posY, double posZ,
                   double dirX, double dirY, double dirZ);

    void getLastNormal(double normal[3]) const
    {
      normal[0] = lastNormal_.x(); normal[1] = lastNormal_.y(); normal[2] = lastNormal_.z();
    }

    /** @brief Get material parameters in current material
     */
    void getMaterialParameters(double& density,
                               double& Z,
                               double& A,
                               double& radiationLength,
                               double& mEE);

    void getMaterialParameters(genfit::MaterialProperties& parameters);

    /** @brief Make a step (following the curvature) until step length
     * sMax or the next boundary is reached.  After making a step to a
     * boundary, the position has to be beyond the boundary, i.e. the
     * current material has to be that beyond the boundary.  The actual
     * step made is returned.
     */
    double findNextBoundary(const genfit::AbsTrackRep::internalExtrapolator& extrap,
                            double sMax);

    double findNextBoundaryStraightLine(double sMax);

  private:

    std::unique_ptr<class G4SafeNavigator> nav_;
    const class G4VPhysicalVolume* currentVolume_;
    G4ThreeVector lastNormal_;
  };

}
