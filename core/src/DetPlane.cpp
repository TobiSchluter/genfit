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

#include <assert.h>
#include <iostream>
#include <cmath>
#include <TMath.h>

#include "DetPlane.h"

namespace genfit {

DetPlane::DetPlane(const TVector3& o,
                       const TVector3& u,
                       const TVector3& v,
                       AbsFinitePlane* finite)
  :o_(o), u_(u), v_(v), finitePlane_(finite)
{
  sane();
}


DetPlane::DetPlane(AbsFinitePlane* finite)
  :finitePlane_(finite)
{
  // default constructor
  o_.SetXYZ(0.,0.,0.);
  u_.SetXYZ(1.,0.,0.);
  v_.SetXYZ(0.,1.,0.);
  // sane() not needed here
}

DetPlane::DetPlane(const TVector3& o,
                       const TVector3& n,
                       AbsFinitePlane* finite)
  :o_(o), finitePlane_(finite)
{
  setNormal(n);
}


DetPlane::~DetPlane(){
  if(finitePlane_ != nullptr)
    delete finitePlane_;
}


DetPlane::DetPlane(const DetPlane& rhs) {
  if(rhs.finitePlane_ != nullptr)
    finitePlane_ = rhs.finitePlane_->clone();
  else finitePlane_ = nullptr;
  o_ = rhs.o_;
  u_ = rhs.u_;
  v_ = rhs.v_;
}


DetPlane& DetPlane::operator=(const DetPlane& rhs) {
  if (this == &rhs)
    return *this;
  if(finitePlane_ != nullptr) {
    delete finitePlane_;
  }
  if(rhs.finitePlane_ != nullptr){
    finitePlane_ = rhs.finitePlane_->clone();
  }
  else{
    finitePlane_ = nullptr;
  }
  o_ = rhs.o_;
  u_ = rhs.u_;
  v_ = rhs.v_;
  return *this;
}


void DetPlane::set(const TVector3& o,
                const TVector3& u,
                const TVector3& v)
{
  o_ = o;
  u_ = u;
  v_ = v;
  sane();
}


void DetPlane::setO(const TVector3& o)
{
  o_ = o;
}

void DetPlane::setO(double X,double Y,double Z)
{
  o_.SetXYZ(X,Y,Z);
}

void DetPlane::setU(const TVector3& u)
{
  u_ = u;
  sane(); // sets v_ perpendicular to u_
}

void DetPlane::setU(double X,double Y,double Z)
{
  u_.SetXYZ(X,Y,Z);
  sane(); // sets v_ perpendicular to u_
}

void DetPlane::setV(const TVector3& v)
{
  v_ = v;
  u_ = getNormal().Cross(v_);
  u_ *= -1.;
  sane();
}

void DetPlane::setV(double X,double Y,double Z)
{
  v_.SetXYZ(X,Y,Z);
  u_ = getNormal().Cross(v_);
  u_ *= -1.;
  sane();
}

void DetPlane::setUV(const TVector3& u,const TVector3& v)
{
  u_ = u;
  v_ = v;
  sane();
}

void DetPlane::setON(const TVector3& o,const TVector3& n){
  o_ = o;
  setNormal(n);
}


TVector3 DetPlane::getNormal() const
{
  return u_.Cross(v_);
}

void DetPlane::setNormal(double X,double Y,double Z){
  setNormal( TVector3(X,Y,Z) );
}

void DetPlane::setNormal(const TVector3& n){
  u_ = n.Orthogonal();
  u_.SetMag(1.);
  v_ = n.Cross(u_);
  v_.SetMag(1.);
}

void DetPlane::setNormal(const double& theta, const double& phi){
  setNormal( TVector3(TMath::Sin(theta)*TMath::Cos(phi),TMath::Sin(theta)*TMath::Sin(phi),TMath::Cos(theta)) );
}


TVector2 DetPlane::project(const TVector3& x)const
{
  return TVector2(u_*x, v_*x);
}


TVector2 DetPlane::LabToPlane(const TVector3& x)const
{
  return project(x-o_);
}


TVector3 DetPlane::toLab(const TVector2& x)const
{
  TVector3 d(o_);
  d += x.X()*u_;
  d += x.Y()*v_;
  return d;
}


TVector3 DetPlane::dist(const TVector3& x)const
{
  return toLab(LabToPlane(x)) - x;
}


void DetPlane::sane(){
  assert(u_!=v_);

  // ensure unit vectors
  u_.SetMag(1.);
  v_.SetMag(1.);

  // check if already orthogonal
  if (u_.Dot(v_) < 1.E-5) return;

  // ensure orthogonal system
  v_ = getNormal().Cross(u_);
}


void DetPlane::Print(const Option_t* option) const
{
  std::cout<<"DetPlane: "
     <<"O("<<o_.X()<<","<<o_.Y()<<","<<o_.Z()<<") "
     <<"u("<<u_.X()<<","<<u_.Y()<<","<<u_.Z()<<") "
     <<"v("<<v_.X()<<","<<v_.Y()<<","<<v_.Z()<<") "
     <<"n("<<getNormal().X()<<","<<getNormal().Y()<<","<<getNormal().Z()<<") "
       <<std::endl;
  std::cout << finitePlane_ << std::endl;
  if(finitePlane_ != nullptr)
    finitePlane_->Print(option);

}


/*
  I could write pages of comments about correct equality checking for
  floating point numbers, but: When two planes are as close as 10E-5 cm
  in all nine numbers that define the plane, this will be enough for all
  practical purposes
 */
bool operator== (const DetPlane& lhs, const DetPlane& rhs){
  if (&lhs == &rhs)
    return true;
  static const double detplaneEpsilon = 1.E-5;
  if(
     fabs( (lhs.o_.X()-rhs.o_.X()) ) > detplaneEpsilon  ||
     fabs( (lhs.o_.Y()-rhs.o_.Y()) ) > detplaneEpsilon  ||
     fabs( (lhs.o_.Z()-rhs.o_.Z()) ) > detplaneEpsilon
     ) return false;
  else if(
    fabs( (lhs.u_.X()-rhs.u_.X()) ) > detplaneEpsilon  ||
    fabs( (lhs.u_.Y()-rhs.u_.Y()) ) > detplaneEpsilon  ||
    fabs( (lhs.u_.Z()-rhs.u_.Z()) ) > detplaneEpsilon
    ) return false;
  else if(
    fabs( (lhs.v_.X()-rhs.v_.X()) ) > detplaneEpsilon  ||
    fabs( (lhs.v_.Y()-rhs.v_.Y()) ) > detplaneEpsilon  ||
    fabs( (lhs.v_.Z()-rhs.v_.Z()) ) > detplaneEpsilon
    ) return false;
  return true;
}

bool operator!= (const DetPlane& lhs, const DetPlane& rhs){
  return !(lhs==rhs);
}


double DetPlane::distance(const TVector3& point) const {
  // |(point - o_)*(u_ x v_)|
  return fabs( (point.X()-o_.X()) * (u_.Y()*v_.Z() - u_.Z()*v_.Y()) +
               (point.Y()-o_.Y()) * (u_.Z()*v_.X() - u_.X()*v_.Z()) +
               (point.Z()-o_.Z()) * (u_.X()*v_.Y() - u_.Y()*v_.X()));
}

double DetPlane::distance(double x, double y, double z) const {
  // |(point - o_)*(u_ x v_)|
  return fabs( (x-o_.X()) * (u_.Y()*v_.Z() - u_.Z()*v_.Y()) +
               (y-o_.Y()) * (u_.Z()*v_.X() - u_.X()*v_.Z()) +
               (z-o_.Z()) * (u_.X()*v_.Y() - u_.Y()*v_.X()));
}


TVector2 DetPlane::straightLineToPlane (const TVector3& point, const TVector3& dir) const {
  TVector3 dirNorm(dir.Unit());
  TVector3 normal = getNormal();
  double dirTimesN = dirNorm*normal;
  if(fabs(dirTimesN)<1.E-6){//straight line is parallel to plane, so return infinity
    return TVector2(1.E100,1.E100);
  }
  double t = 1./dirTimesN * ((o_-point)*normal);
  return project(point - o_ + t * dirNorm);
}


//! gives u,v coordinates of the intersection point of a straight line with plane
void DetPlane::straightLineToPlane(const double& posX, const double& posY, const double& posZ,
                                   const double& dirX, const double& dirY, const double& dirZ,
                                   double& u, double& v) const {

  TVector3 W = getNormal();
  double dirTimesN = dirX*W.X() + dirY*W.Y() + dirZ*W.Z();
  if(fabs(dirTimesN)<1.E-6){//straight line is parallel to plane, so return infinity
    u = 1.E100;
    v = 1.E100;
    return;
  }
  double t = 1./dirTimesN * ((o_.X()-posX)*W.X() +
                             (o_.Y()-posY)*W.Y() +
                             (o_.Z()-posZ)*W.Z());

  double posOnPlaneX = posX-o_.X() + t*dirX;
  double posOnPlaneY = posY-o_.Y() + t*dirY;
  double posOnPlaneZ = posZ-o_.Z() + t*dirZ;

  u = u_.X()*posOnPlaneX + u_.Y()*posOnPlaneY + u_.Z()*posOnPlaneZ;
  v = v_.X()*posOnPlaneX + v_.Y()*posOnPlaneY + v_.Z()*posOnPlaneZ;
}


void DetPlane::reset() {
  o_.SetXYZ(0.,0.,0.);
  u_.SetXYZ(1.,0.,0.);
  v_.SetXYZ(0.,1.,0.);
  if(finitePlane_ != nullptr) {
    delete finitePlane_;
    finitePlane_ = nullptr;
  }
}

} /* End of namespace genfit */
