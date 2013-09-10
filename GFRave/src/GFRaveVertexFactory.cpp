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

#include <string>

#include "GFRaveVertexFactory.h"
#include "GFRaveConverters.h"
#include "GFRaveVertex.h"

#include "GFRaveMagneticField.h"
#include "GFRavePropagator.h"

#include "Exception.h"

#include "rave/Propagator.h"
#include "rave/MagneticField.h"
#include "rave/VertexFactory.h"
#include "rave/Vertex.h"
#include "rave/Ellipsoid3D.h"


namespace genfit {

GFRaveVertexFactory::GFRaveVertexFactory(int verbosity, bool useVacuumPropagator) {
  IdGFTrackMap_ = new std::map<int, genfit::Track*>;
  IdGFMeasuredStateOnPlaneMap_ = new std::map<int, genfit::MeasuredStateOnPlane*>;

  if (useVacuumPropagator) {
    propagator_ = new rave::VacuumPropagator();
  }
  else {
    propagator_ = new GFRavePropagator();
    (static_cast<GFRavePropagator*>(propagator_))->setIdGFMeasuredStateOnPlaneMap(IdGFMeasuredStateOnPlaneMap_);
  }

  magneticField_ = new GFRaveMagneticField();

  if (verbosity > 0) ++verbosity; // verbosity has to be >1 for rave

  factory_ = new rave::VertexFactory(*magneticField_, *propagator_, "kalman-smoothing:1", verbosity); // here copies of magneticField_ and propagator_ are made!
}


GFRaveVertexFactory::~GFRaveVertexFactory(){
  clearMaps();
  delete IdGFTrackMap_;
  delete IdGFMeasuredStateOnPlaneMap_;

  delete magneticField_;
  delete propagator_;
  delete factory_;
}


void
GFRaveVertexFactory::findVertices ( std::vector <  genfit::GFRaveVertex* > * GFvertices, const std::vector < genfit::Track* > & GFTracks, bool use_beamspot ){
  clearMaps();

  try{
    RaveToGFVertices(GFvertices,
                             factory_->create(GFTracksToTracks(GFTracks, IdGFTrackMap_, IdGFMeasuredStateOnPlaneMap_, 0),
                                              use_beamspot),
                             IdGFTrackMap_, IdGFMeasuredStateOnPlaneMap_);
  }
  catch(Exception & e){
    std::cerr << e.what();
    return;
  }
}


void
GFRaveVertexFactory::findVertices ( std::vector <  genfit::GFRaveVertex* > * GFvertices, const std::vector < genfit::MeasuredStateOnPlane* > & GFMeasuredStateOnPlanes, bool use_beamspot ){
  clearMaps();

  try{
    RaveToGFVertices(GFvertices,
                     factory_->create(MeasuredStateOnPlanesToTracks(GFMeasuredStateOnPlanes, IdGFTrackMap_, IdGFMeasuredStateOnPlaneMap_, 0),
                                      use_beamspot),
                     IdGFTrackMap_, IdGFMeasuredStateOnPlaneMap_);
  }
  catch(Exception & e){
    std::cerr << e.what();
    return;
  }
}


void
GFRaveVertexFactory::setBeamspot(const TVector3 & pos, const TMatrixDSym & cov){
  factory_->setBeamSpot(rave::Ellipsoid3D(TVector3ToPoint3D(pos),
                        TMatrixDSymToCovariance3D(cov)));
}


void
GFRaveVertexFactory::setMethod(const std::string & method){
  size_t found = method.find("smoothing:1");
  if (found==std::string::npos){
    std::cerr << "GFRaveVertexFactory::setMethod(" << method << ") ==> smoothing not turned on! GFRaveTrackParameters will be unsmoothed!" << std::endl;
  }
  factory_->setDefaultMethod(method);
  std::cout << "GFRaveVertexFactory::setMethod ==> set method to " << factory_->method() << std::endl;
}


void
GFRaveVertexFactory::clearMaps(){
  IdGFTrackMap_->clear();

  for (unsigned int i=0; i<IdGFMeasuredStateOnPlaneMap_->size(); ++i){
    // in here are copies or trackreps -> we have ownership
    delete (*IdGFMeasuredStateOnPlaneMap_)[i];
  }
  IdGFMeasuredStateOnPlaneMap_->clear();
}


} /* End of namespace genfit */
