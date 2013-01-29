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
#include "GFDetPlane.h"

#include <assert.h>
#include <iostream>
#include <cmath>
#include <TMath.h>
#include <TPolyMarker3D.h>
#include <TPolyLine3D.h>

ClassImp(GFDetPlane)

GFDetPlane::GFDetPlane(const TVector3& o,
                       const TVector3& u,
                       const TVector3& v,
                       GFAbsFinitePlane* finite)
  :fO(o), fU(u), fV(v), fFinitePlane(finite)
{
  sane();
}


GFDetPlane::GFDetPlane(GFAbsFinitePlane* finite) 
  :fFinitePlane(finite)
{
  // default constructor
  fO.SetXYZ(0.,0.,0.);
  fU.SetXYZ(1.,0.,0.);
  fV.SetXYZ(0.,1.,0.);
  // sane() not needed here
}

GFDetPlane::GFDetPlane(const TVector3& o,
                       const TVector3& n,
                       GFAbsFinitePlane* finite)
  :fO(o), fFinitePlane(finite)
{
  setNormal(n);
}


GFDetPlane::~GFDetPlane(){
  if(fFinitePlane!=NULL) delete fFinitePlane;
}


GFDetPlane::GFDetPlane(const GFDetPlane& rhs) {
  if(rhs.fFinitePlane != NULL) fFinitePlane = rhs.fFinitePlane->clone();
  else fFinitePlane = NULL;
  fO = rhs.fO;
  fU = rhs.fU;
  fV = rhs.fV;
}


GFDetPlane& GFDetPlane::operator=(const GFDetPlane& rhs) {
  if (this == &rhs)
    return *this;
  if(fFinitePlane != NULL) {
    delete fFinitePlane;
  }
  if(rhs.fFinitePlane != NULL){
    fFinitePlane = rhs.fFinitePlane->clone();
  }
  else{
    fFinitePlane = NULL;
  }
  fO = rhs.fO;
  fU = rhs.fU;
  fV = rhs.fV;
  return *this;
}


void 
GFDetPlane::set(const TVector3& o,
                const TVector3& u,
                const TVector3& v)
{
  fO = o;
  fU = u;
  fV = v;
  sane();
}


void 
GFDetPlane::setO(const TVector3& o)
{
  fO = o;
}

void 
GFDetPlane::setO(double X,double Y,double Z)
{
  fO.SetXYZ(X,Y,Z);
}

void 
GFDetPlane::setU(const TVector3& u)
{
  fU = u;
  sane(); // sets fV perpendicular to fU
}

void 
GFDetPlane::setU(double X,double Y,double Z)
{
  fU.SetXYZ(X,Y,Z);
  sane(); // sets fV perpendicular to fU
}

void 
GFDetPlane::setV(const TVector3& v)
{
  fV = v;
  fU = getNormal().Cross(fV);
  fU *= -1.;
  sane();
}

void 
GFDetPlane::setV(double X,double Y,double Z)
{
  fV.SetXYZ(X,Y,Z);
  fU = getNormal().Cross(fV);
  fU *= -1.;
  sane();
}

void 
GFDetPlane::setUV(const TVector3& u,const TVector3& v)
{
  fU = u;
  fV = v;
  sane();
}

void
GFDetPlane::setON(const TVector3& o,const TVector3& n){
  fO = o;
  setNormal(n);
}


TVector3
GFDetPlane::getNormal() const
{
  return fU.Cross(fV);
}

void
GFDetPlane::setNormal(double X,double Y,double Z){
  setNormal( TVector3(X,Y,Z) );
}

void
GFDetPlane::setNormal(const TVector3& n){
  fU = n.Orthogonal();
  fU.SetMag(1.);
  fV = n.Cross(fU);
  fV.SetMag(1.);
}

void GFDetPlane::setNormal(const double& theta, const double& phi){
  setNormal( TVector3(TMath::Sin(theta)*TMath::Cos(phi),TMath::Sin(theta)*TMath::Sin(phi),TMath::Cos(theta)) );
}


TVector2
GFDetPlane::project(const TVector3& x)const
{
  return TVector2(fU*x, fV*x);
}


TVector2 
GFDetPlane::LabToPlane(const TVector3& x)const
{
  return project(x-fO);
}


TVector3
GFDetPlane::toLab(const TVector2& x)const
{
  TVector3 d(fO);
  d += x.X()*fU;
  d += x.Y()*fV;
  return d;
}


TVector3 
GFDetPlane::dist(const TVector3& x)const
{
  return toLab(LabToPlane(x)) - x;
}


void
GFDetPlane::sane(){
  assert(fU!=fV);

  // ensure unit vectors
  fU.SetMag(1.);
  fV.SetMag(1.);

  // check if already orthogonal
  if (fU.Dot(fV) < 1.E-5) return;

  // ensure orthogonal system
  fV = getNormal().Cross(fU);
}


void
GFDetPlane::Print(const Option_t* option) const
{
  std::cout<<"GFDetPlane: "
	   <<"O("<<fO.X()<<","<<fO.Y()<<","<<fO.Z()<<") "
	   <<"u("<<fU.X()<<","<<fU.Y()<<","<<fU.Z()<<") "
	   <<"v("<<fV.X()<<","<<fV.Y()<<","<<fV.Z()<<") "
	   <<"n("<<getNormal().X()<<","<<getNormal().Y()<<","<<getNormal().Z()<<") "
		   <<std::endl;
  std::cout << fFinitePlane << std::endl;
  if(fFinitePlane!=NULL) fFinitePlane->Print(option);
    
}


/*
  I could write pages of comments about correct equality checking for 
  floating point numbers, but: When two planes are as close as 10E-5 cm
  in all nine numbers that define the plane, this will be enough for all
  practical purposes
 */
bool operator== (const GFDetPlane& lhs, const GFDetPlane& rhs){
  if (&lhs == &rhs)
    return true;
  static const double detplaneEpsilon = 1.E-5;
  if(
     fabs( (lhs.fO.X()-rhs.fO.X()) ) > detplaneEpsilon  ||
     fabs( (lhs.fO.Y()-rhs.fO.Y()) ) > detplaneEpsilon  ||
     fabs( (lhs.fO.Z()-rhs.fO.Z()) ) > detplaneEpsilon
     ) return false;
  else if(
	  fabs( (lhs.fU.X()-rhs.fU.X()) ) > detplaneEpsilon  ||
	  fabs( (lhs.fU.Y()-rhs.fU.Y()) ) > detplaneEpsilon  ||
	  fabs( (lhs.fU.Z()-rhs.fU.Z()) ) > detplaneEpsilon
	  ) return false;
  else if(
	  fabs( (lhs.fV.X()-rhs.fV.X()) ) > detplaneEpsilon  ||
	  fabs( (lhs.fV.Y()-rhs.fV.Y()) ) > detplaneEpsilon  ||
	  fabs( (lhs.fV.Z()-rhs.fV.Z()) ) > detplaneEpsilon
	  ) return false;
  return true;
}

bool operator!= (const GFDetPlane& lhs, const GFDetPlane& rhs){
  return !(lhs==rhs);
}


double GFDetPlane::distance(const TVector3& point) const {
  // |(point - fO)*(fU x fV)|
  return fabs( (point.X()-fO.X()) * (fU.Y()*fV.Z() - fU.Z()*fV.Y()) +
               (point.Y()-fO.Y()) * (fU.Z()*fV.X() - fU.X()*fV.Z()) +
               (point.Z()-fO.Z()) * (fU.X()*fV.Y() - fU.Y()*fV.X()));
}

double GFDetPlane::distance(double x, double y, double z) const {
  // |(point - fO)*(fU x fV)|
  return fabs( (x-fO.X()) * (fU.Y()*fV.Z() - fU.Z()*fV.Y()) +
               (y-fO.Y()) * (fU.Z()*fV.X() - fU.X()*fV.Z()) +
               (z-fO.Z()) * (fU.X()*fV.Y() - fU.Y()*fV.X()));
}


TVector2 GFDetPlane::straightLineToPlane (const TVector3& point, const TVector3& dir) const {
  TVector3 dirNorm(dir.Unit());
  TVector3 normal = getNormal();
  double dirTimesN = dirNorm*normal;
  if(fabs(dirTimesN)<1.E-6){//straight line is parallel to plane, so return infinity
    return TVector2(1.E100,1.E100);
  }
  double t = 1./dirTimesN * ((fO-point)*normal);
  return project(point - fO + t * dirNorm);
}


void GFDetPlane::reset() {
  fO.SetXYZ(0.,0.,0.);
  fU.SetXYZ(1.,0.,0.);
  fV.SetXYZ(0.,1.,0.);
  if(fFinitePlane!=NULL) {
    delete fFinitePlane;
    fFinitePlane = NULL;
  }
}

