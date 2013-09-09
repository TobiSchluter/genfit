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


#include "GFRavePropagator.h"
#include "GFRaveConverters.h"
#include "Exception.h"

#include <iostream>


namespace genfit {

GFRavePropagator::GFRavePropagator() :
    IdGFMeasuredStateOnPlaneMap_(NULL)
{}
  
GFRavePropagator*
GFRavePropagator::copy() const
{
  return new GFRavePropagator(*this);
}


GFRavePropagator::~GFRavePropagator()
{

}
    
    
std::pair < rave::Track, double >
GFRavePropagator::to ( const rave::Track & orig,
                       const ravesurf::Cylinder & rcyl ) const
{
  // todo to be implemented!!
  Exception exc("GFRavePropagator::to (cylinder) ==> not yet implemented!",__LINE__,__FILE__);
  throw exc;
}


std::pair < rave::Track, double >
GFRavePropagator::to ( const rave::Track & orig,
                       const ravesurf::Plane & rplane ) const
{
  MeasuredStateOnPlane* state = getMeasuredStateOnPlane(orig);

  // will throw Exception if extrapolation does not work
  double path = state->extrapolateToPlane(PlaneToGFDetPlane(rplane));

  std::pair < rave::Track, double > ret(MeasuredStateOnPlaneToTrack(state, orig), path);
  return ret;
}


rave::Track
GFRavePropagator::closestTo ( const rave::Track & orig,
                              const rave::Point3D & pt, bool transverse ) const
{

  if (transverse){
    Exception exc("GFRavePropagator::closestTo ==> transverse is true, not implemented!",__LINE__,__FILE__);
    throw exc;
  }

  MeasuredStateOnPlane* state = getMeasuredStateOnPlane(orig);

  TVector3 point(Point3DToTVector3(pt));
  state->extrapolateToPoint(point);

  return MeasuredStateOnPlaneToTrack(state, orig);
}


MeasuredStateOnPlane*
GFRavePropagator::getMeasuredStateOnPlane(const rave::Track & track) const {

  if (IdGFMeasuredStateOnPlaneMap_==NULL) {
    Exception exc("GFRavePropagator::getTrackRep ==> IdGFMeasuredStateOnPlaneMap_ is NULL, cannot access genfit::Tracks!",__LINE__,__FILE__);
    throw exc;
  }

  if (!(track.isValid())) {
    Exception exc("GFRavePropagator::getTrackRep ==> rave::Track is not valid!",__LINE__,__FILE__);
    throw exc;
  }

  //std::cerr<<"GFRavePropagator::getTrackRep track id: "<<track.id()<<std::endl;
  //std::cerr<<"  pos: "; Point3DToTVector3(track.state().position()).Print();
  //std::cerr<<"  mom: "; Vector3DToTVector3(track.state().momentum()).Print();

  if (IdGFMeasuredStateOnPlaneMap_->count(track.id()) == 0) {
    Exception exc("GFRavePropagator::getTrackRep ==> no entry in IdGFMeasuredStateOnPlaneMap_ corresponding to track id, cannot access corresponding state!",__LINE__,__FILE__);
    throw exc;
  }

  MeasuredStateOnPlane* state = IdGFMeasuredStateOnPlaneMap_->at(track.id());

  setData(track, state); // set state and cov

  return state;
}


void
GFRavePropagator::setIdGFMeasuredStateOnPlaneMap(std::map < int, MeasuredStateOnPlane* > * map){
  if (map==NULL) {
    Exception exc("GFRavePropagator::setIdGFMeasuredStateOnPlaneMap ==> map is NULL!",__LINE__,__FILE__);
    throw exc;
  }
  IdGFMeasuredStateOnPlaneMap_ = map;
  //std::cout<<"IdGFMeasuredStateOnPlaneMap_: " << (int)IdGFMeasuredStateOnPlaneMap_ << std::endl;
}

} /* End of namespace genfit */
