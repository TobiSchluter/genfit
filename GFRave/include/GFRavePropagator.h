/* Copyright 2008-2010, Technische Universitaet Muenchen,
   Authors: Christian Hoeppner & Sebastian Neubert & Johannes Rauch

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

/**
 *  @author Johannes Rauch (Technische Universit&auml;t M&uuml;nchen, original author)
 *
 */

/** @addtogroup GFRave
 * @{
 */

#ifndef GFRAVEPROPAGATOR_H
#define GFRAVEPROPAGATOR_H

#include "rave/Propagator.h"
#include "rave/Track.h"
#include "rave/Plane.h"
#include "rave/Cylinder.h"

#include "AbsTrackRep.h"

#include <map>


namespace genfit {

/**
 * @brief GFRavePropagator class
 * Inherits from rave::Propagator. A map has to be provided,
 * containing pointers to clones of the Cardinal Reps of the GFtracks.
 * The GFRavePropagator uses the information of the rave track to set
 * the state and covariance of the corresponding AbsTrackRep,
 * extrapolates and then returns a new rave track with the
 * extrapolated state and covariance.
 *
 */

class GFRavePropagator : public rave::Propagator
{
  public:
    GFRavePropagator();
    virtual GFRavePropagator* copy() const;
    virtual rave::Track closestTo ( const rave::Track &,
                                    const rave::Point3D &, bool transverse ) const;
    virtual std::pair < rave::Track, double > to ( const rave::Track & orig,
                                                   const ravesurf::Plane & ) const;
    virtual std::pair < rave::Track, double > to ( const rave::Track & orig,
                                                   const ravesurf::Cylinder & ) const;

    virtual ~GFRavePropagator();

    void setIdGFMeasuredStateOnPlaneMap(std::map < int, MeasuredStateOnPlane* > * map);

  private:


    // check if everything is ok, otherwise throw Exception;
    // get MeasuredStateOnPlane from IdGFMeasuredStateOnPlaneMap_ and set MeasuredStateOnPlane state and cov from track
    MeasuredStateOnPlane* getMeasuredStateOnPlane(const rave::Track & track) const;

    // data members
    std::map < int, MeasuredStateOnPlane* > * IdGFMeasuredStateOnPlaneMap_; // pointers to GFAbsTracksReps via rave track ID
};


} /* End of namespace genfit */

#endif
