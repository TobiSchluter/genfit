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


#ifndef GFRAVEVERTEXFACTORY_H
#define GFRAVEVERTEXFACTORY_H

#include <vector>

#include "GFRaveVertex.h"

#include "rave/Propagator.h"
#include "rave/MagneticField.h"
#include "rave/VertexFactory.h"

#include "Track.h"

namespace genfit {


/**
 * @brief GFRaveVertexFactory
 *  The GFRaveVertexFactory is basically a wrapper around the rave::VertexFactory.
 *  It takes care of initializing the rave::VertexFactory, building the necessary maps,
 *  convert GENFIT to rave objects and vice versa.
 **/
class GFRaveVertexFactory {
 public:
  // constructors, destructors
  GFRaveVertexFactory(int verbosity = 0, bool useVacuumPropagator = false);
  ~GFRaveVertexFactory();

  // functions

  void findVertices ( std::vector <  genfit::GFRaveVertex* > *, const std::vector < genfit::Track* > &, bool use_beamspot=false );
  void findVertices ( std::vector <  genfit::GFRaveVertex* > *, const std::vector < genfit::MeasuredStateOnPlane* > &, bool use_beamspot=false );

  void setBeamspot(const TVector3 & pos, const TMatrixDSym & cov);

  /**
   * Set the reconstruction method. See http://projects.hepforge.org/rave/trac/wiki/RaveMethods
   * Smoothing has to be turned on! e.g. kalman-smoothing:1
   */
  void setMethod(const std::string & method);

 private:

  void clearMaps();

  // data members
  rave::VertexFactory* factory_;
  rave::MagneticField* magneticField_;
  rave::Propagator* propagator_;

  /**
   * bookkeeping of original genfit::Tracks for later assignment to GFVertices
   */
  std::map<int, genfit::Track*> * IdGFTrackMap_;
  /**
   * map of copies of the cardinal reps for the GFRavePropagator; ownership of trackrep clones is HERE!!!
   */
  std::map<int, genfit::MeasuredStateOnPlane*> * IdGFMeasuredStateOnPlaneMap_;


};

} /* End of namespace genfit */

#endif

/** @} */

