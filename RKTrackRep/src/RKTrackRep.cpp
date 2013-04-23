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

#include "RKTrackRep.h"

#include <Exception.h>
#include <FieldManager.h>
#include <MaterialEffects.h>

#include <TDatabasePDG.h>
#include <TMath.h>


#define MINSTEP 0.001   // minimum step [cm] for Runge Kutta and iteration to POCA
//#define DEBUG


namespace genfit {


void RKTrackRep::initArrays(){
  memset(fJ_pM_5x7,0x00,5*7*sizeof(double));
  memset(fJ_pM_5x6,0x00,5*6*sizeof(double));
  memset(fJ_Mp_7x5,0x00,7*5*sizeof(double));
  memset(fJ_Mp_6x5,0x00,6*5*sizeof(double));
}


void RKTrackRep::setPosMom(StateOnPlane* state, const TVector3& pos, const TVector3& mom) const {

  // TODO: test!

  if (state->getRep() != this){
    Exception exc("RKTrackRep::setPosMom ==> state is defined wrt. another TrackRep",__LINE__,__FILE__);
    throw exc;
  }

  if (state->getPlane()->distance(pos) < MINSTEP) { // pos is on plane -> do not change plane!

    M1x7 state7;

    state7[0] = pos.X();
    state7[1] = pos.Y();
    state7[2] = pos.Z();

    state7[3] = mom.X();
    state7[4] = mom.Y();
    state7[5] = mom.Z();

    // normalize dir
    double norm = 1. / sqrt(state7[3]*state7[3] + state7[4]*state7[4] + state7[5]*state7[5]);
    for (unsigned int i=3; i<6; ++i)
      state7[i] *= norm;

    state7[6] = getCharge(state) * norm;

    getState5(state, state7);

  }
  else { // pos is not on plane -> create new plane!

    // TODO: Raise Warning that a new plane has been created!
    SharedPlanePtr plane(new DetPlane(pos, mom));
    state->setPlane(plane);

    TVectorD& state5(state->getState());

    state5(0) = getCharge(state)/mom.Mag(); // q/p
    state5(1) = 0.; // u'
    state5(2) = 0.; // v'
    state5(3) = 0.; // u
    state5(4) = 0.; // v

    setSpu(state, 1.);
  }

}


void RKTrackRep::setPosMomErr(MeasuredStateOnPlane* state, const TVector3& pos, const TVector3& mom, const TVector3& posErr, const TVector3& momErr) const {

  // TODO: test!

  setPosMom(state, pos, mom);

  const TVector3& U(state->getPlane()->getU());
  const TVector3& V(state->getPlane()->getV());
  TVector3 W(state->getPlane()->getNormal());

  double pw = mom * W;
  double pu = mom * U;
  double pv = mom * V;

  TMatrixDSym& cov(state->getCov());

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




void RKTrackRep::setPosMomCov(MeasuredStateOnPlane* state, const TVector3& pos, const TVector3& mom, const TMatrixDSym& cov6x6) const {

  if (cov6x6.GetNcols()!=6 || cov6x6.GetNrows()!=6){
    Exception exc("RKTrackRep::setPosMomCov ==> cov has to be 6x6 (x, y, z, px, py, pz)",__LINE__,__FILE__);
    throw exc;
  }

  if (fabs(getCharge(state)) < 0.001){
    Exception exc("RKTrackRep::setPosMomCov ==> charge is 0. setPosMomCov cannot work with charge == 0 ",__LINE__,__FILE__);
    throw exc;
  }

  setPosMom(state, pos, mom); // charge does not change!

  M1x7 state7;
  getState7(state, state7);

  const M6x6& cov6x6_( *((M6x6*) cov6x6.GetMatrixArray()) );

  transformM6P(cov6x6_, state7, state);

}



double RKTrackRep::RKPropagate(M1x7& state7,
                        M7x7* cov,
                        M1x3& SA,
                        double S,
                        bool varField) const {

  // important fixed numbers
  static const double EC     = 0.000149896229;  // c/(2*10^12) resp. c/2Tera
  static const double P3     = 1./3.;           // 1/3
  static const double DLT    = .0002;           // max. deviation for approximation-quality test
  static const double par = 1./3.081615;
  // Aux parameters
  M1x3&   R           = *((M1x3*) &state7[0]);       // Start coordinates  [cm]  (x,  y,  z)
  M1x3&   A           = *((M1x3*) &state7[3]);       // Start directions         (ax, ay, az);   ax^2+ay^2+az^2=1
  double  S3(0), S4(0), PS2(0);
  M1x3     H0 = {0.,0.,0.}, H1 = {0.,0.,0.}, H2 = {0.,0.,0.}, r = {0.,0.,0.};
  // Variables for RKutta main loop
  double   A0(0), A1(0), A2(0), A3(0), A4(0), A5(0), A6(0);
  double   B0(0), B1(0), B2(0), B3(0), B4(0), B5(0), B6(0);
  double   C0(0), C1(0), C2(0), C3(0), C4(0), C5(0), C6(0);

  bool calcCov(cov != NULL);

  //
  // Runge Kutta Extrapolation
  //
  S3 = P3*S;
  S4 = 0.25*S;
  PS2 = state7[6]*EC * S;

  // First point
  r[0] = R[0];           r[1] = R[1];           r[2]=R[2];
  TVector3 pos(r[0], r[1], r[2]);// vector of start coordinates R0  (x, y, z)
  TVector3 field(FieldManager::getFieldVal(pos));       // magnetic field in 10^-4 T = kGauss
  H0[0] = PS2*field.X(); H0[1] = PS2*field.Y(); H0[2] = PS2*field.Z();     // H0 is PS2*(Hx, Hy, Hz) @ R0
  A0 = A[1]*H0[2]-A[2]*H0[1]; B0 = A[2]*H0[0]-A[0]*H0[2]; C0 = A[0]*H0[1]-A[1]*H0[0]; // (ax, ay, az) x H0
  A2 = A[0]+A0              ; B2 = A[1]+B0              ; C2 = A[2]+C0              ; // (A0, B0, C0) + (ax, ay, az)
  A1 = A2+A[0]              ; B1 = B2+A[1]              ; C1 = C2+A[2]              ; // (A0, B0, C0) + 2*(ax, ay, az)

  // Second point
  if (varField) {
    r[0] += A1*S4;         r[1] += B1*S4;         r[2] += C1*S4;
    pos.SetXYZ(r[0], r[1], r[2]);
    field = FieldManager::getFieldVal(pos);
    H1[0] = field.X()*PS2; H1[1] = field.Y()*PS2; H1[2] = field.Z()*PS2; // H1 is PS2*(Hx, Hy, Hz) @ (x, y, z) + 0.25*S * [(A0, B0, C0) + 2*(ax, ay, az)]
  }
  else if (calcCov) memcpy(H1, H0, 3*sizeof(double));
  A3 = B2*H1[2]-C2*H1[1]+A[0]; B3 = C2*H1[0]-A2*H1[2]+A[1]; C3 = A2*H1[1]-B2*H1[0]+A[2]; // (A2, B2, C2) x H1 + (ax, ay, az)
  A4 = B3*H1[2]-C3*H1[1]+A[0]; B4 = C3*H1[0]-A3*H1[2]+A[1]; C4 = A3*H1[1]-B3*H1[0]+A[2]; // (A3, B3, C3) x H1 + (ax, ay, az)
  A5 = A4-A[0]+A4            ; B5 = B4-A[1]+B4            ; C5 = C4-A[2]+C4            ; //    2*(A4, B4, C4) - (ax, ay, az)

  // Last point
  if (varField) {
    r[0]=R[0]+S*A4;         r[1]=R[1]+S*B4;         r[2]=R[2]+S*C4;  //setup.Field(r,H2);
    pos.SetXYZ(r[0], r[1], r[2]);
    field = FieldManager::getFieldVal(pos);
    H2[0] = field.X()*PS2;  H2[1] = field.Y()*PS2;  H2[2] = field.Z()*PS2; // H2 is PS2*(Hx, Hy, Hz) @ (x, y, z) + 0.25*S * (A4, B4, C4)
  }
  else if (calcCov) memcpy(H2, H0, 3*sizeof(double));
  A6 = B5*H2[2]-C5*H2[1]; B6 = C5*H2[0]-A5*H2[2]; C6 = A5*H2[1]-B5*H2[0]; // (A5, B5, C5) x H2


  //
  // Derivatives of track parameters
  //
  if(calcCov){
    double   dA0(0), dA2(0), dA3(0), dA4(0), dA5(0), dA6(0);
    double   dB0(0), dB2(0), dB3(0), dB4(0), dB5(0), dB6(0);
    double   dC0(0), dC2(0), dC3(0), dC4(0), dC5(0), dC6(0);

    // d(x, y, z)/d(x, y, z) submatrix is unit matrix
    (*cov)[0] = 1;  (*cov)[8] = 1;  (*cov)[16] = 1;
    // d(ax, ay, az)/d(ax, ay, az) submatrix is 0
    // start with d(x, y, z)/d(ax, ay, az)
    for(int i=3*7; i<49; i+=7) {

      if(i==42) {(*cov)[i+3]*=state7[6]; (*cov)[i+4]*=state7[6]; (*cov)[i+5]*=state7[6];}

      //first point
      dA0 = H0[2]*(*cov)[i+4]-H0[1]*(*cov)[i+5];    // dA0/dp }
      dB0 = H0[0]*(*cov)[i+5]-H0[2]*(*cov)[i+3];    // dB0/dp  } = dA x H0
      dC0 = H0[1]*(*cov)[i+3]-H0[0]*(*cov)[i+4];    // dC0/dp }

      if(i==42) {dA0+=A0; dB0+=B0; dC0+=C0;}     // if last row: (dA0, dB0, dC0) := (dA0, dB0, dC0) + (A0, B0, C0)

      dA2 = dA0+(*cov)[i+3];        // }
      dB2 = dB0+(*cov)[i+4];        //  } = (dA0, dB0, dC0) + dA
      dC2 = dC0+(*cov)[i+5];        // }

      //second point
      dA3 = (*cov)[i+3]+dB2*H1[2]-dC2*H1[1];    // dA3/dp }
      dB3 = (*cov)[i+4]+dC2*H1[0]-dA2*H1[2];    // dB3/dp  } = dA + (dA2, dB2, dC2) x H1
      dC3 = (*cov)[i+5]+dA2*H1[1]-dB2*H1[0];    // dC3/dp }

      if(i==42) {dA3+=A3-A[0]; dB3+=B3-A[1]; dC3+=C3-A[2];} // if last row: (dA3, dB3, dC3) := (dA3, dB3, dC3) + (A3, B3, C3) - (ax, ay, az)

      dA4 = (*cov)[i+3]+dB3*H1[2]-dC3*H1[1];    // dA4/dp }
      dB4 = (*cov)[i+4]+dC3*H1[0]-dA3*H1[2];    // dB4/dp  } = dA + (dA3, dB3, dC3) x H1
      dC4 = (*cov)[i+5]+dA3*H1[1]-dB3*H1[0];    // dC4/dp }

      if(i==42) {dA4+=A4-A[0]; dB4+=B4-A[1]; dC4+=C4-A[2];} // if last row: (dA4, dB4, dC4) := (dA4, dB4, dC4) + (A4, B4, C4) - (ax, ay, az)

      //last point
      dA5 = dA4+dA4-(*cov)[i+3];      // }
      dB5 = dB4+dB4-(*cov)[i+4];      //  } =  2*(dA4, dB4, dC4) - dA
      dC5 = dC4+dC4-(*cov)[i+5];      // }

      dA6 = dB5*H2[2]-dC5*H2[1];      // dA6/dp }
      dB6 = dC5*H2[0]-dA5*H2[2];      // dB6/dp  } = (dA5, dB5, dC5) x H2
      dC6 = dA5*H2[1]-dB5*H2[0];      // dC6/dp }

      if(i==42) {dA6+=A6; dB6+=B6; dC6+=C6;}     // if last row: (dA6, dB6, dC6) := (dA6, dB6, dC6) + (A6, B6, C6)

      if(i==42) {
        (*cov)[i]   += (dA2+dA3+dA4)*S3/state7[6];  (*cov)[i+3] = (dA0+dA3+dA3+dA5+dA6)*P3/state7[6]; // dR := dR + S3*[(dA2, dB2, dC2) +   (dA3, dB3, dC3) + (dA4, dB4, dC4)]
        (*cov)[i+1] += (dB2+dB3+dB4)*S3/state7[6];  (*cov)[i+4] = (dB0+dB3+dB3+dB5+dB6)*P3/state7[6]; // dA :=     1/3*[(dA0, dB0, dC0) + 2*(dA3, dB3, dC3) + (dA5, dB5, dC5) + (dA6, dB6, dC6)]
        (*cov)[i+2] += (dC2+dC3+dC4)*S3/state7[6];  (*cov)[i+5] = (dC0+dC3+dC3+dC5+dC6)*P3/state7[6];
      }
      else {
        (*cov)[i]   += (dA2+dA3+dA4)*S3;  (*cov)[i+3] = (dA0+dA3+dA3+dA5+dA6)*P3; // dR := dR + S3*[(dA2, dB2, dC2) +   (dA3, dB3, dC3) + (dA4, dB4, dC4)]
        (*cov)[i+1] += (dB2+dB3+dB4)*S3;  (*cov)[i+4] = (dB0+dB3+dB3+dB5+dB6)*P3; // dA :=     1/3*[(dA0, dB0, dC0) + 2*(dA3, dB3, dC3) + (dA5, dB5, dC5) + (dA6, dB6, dC6)]
        (*cov)[i+2] += (dC2+dC3+dC4)*S3;  (*cov)[i+5] = (dC0+dC3+dC3+dC5+dC6)*P3;
      }
    }
  }

  //
  // Track parameters in last point
  //
  R[0] += (A2+A3+A4)*S3;   A[0] += (SA[0]=(A0+A3+A3+A5+A6)*P3-A[0]);  // R  = R0 + S3*[(A2, B2, C2) +   (A3, B3, C3) + (A4, B4, C4)]
  R[1] += (B2+B3+B4)*S3;   A[1] += (SA[1]=(B0+B3+B3+B5+B6)*P3-A[1]);  // A  =     1/3*[(A0, B0, C0) + 2*(A3, B3, C3) + (A5, B5, C5) + (A6, B6, C6)]
  R[2] += (C2+C3+C4)*S3;   A[2] += (SA[2]=(C0+C3+C3+C5+C6)*P3-A[2]);  // SA = A_new - A_old

  // normalize A
  double CBA = 1./sqrt(A[0]*A[0]+A[1]*A[1]+A[2]*A[2]); // 1/|A|
  A[0] *= CBA; A[1] *= CBA; A[2] *= CBA;


  // Test approximation quality on given step
  double EST = fabs((A1+A6)-(A3+A4))+fabs((B1+B6)-(B3+B4))+fabs((C1+C6)-(C3+C4));  // EST = ||(ABC1+ABC6)-(ABC3+ABC4)||_1  =  ||(axzy x H0 + ABC5 x H2) - (ABC2 x H1 + ABC3 x H1)||_1
#ifdef DEBUG
   std::cerr << "   RKTrackRep::RKPropagate quality EST " << EST  << " \n";
#endif
  return pow(DLT/EST, par);
}






void RKTrackRep::getState7(const StateOnPlane* state, M1x7& state7) const {

  const TVector3& U(state->getPlane()->getU());
  const TVector3& V(state->getPlane()->getV());
  const TVector3& O(state->getPlane()->getO());
  TVector3 W(state->getPlane()->getNormal());

  const TVectorD& state5(state->getState());

  double spu = getSpu(state);

  state7[0] = O.X() + state5(3)*U.X() + state5(4)*V.X(); // x
  state7[1] = O.Y() + state5(3)*U.Y() + state5(4)*V.Y(); // y
  state7[2] = O.Z() + state5(3)*U.Z() + state5(4)*V.Z(); // z

  state7[3] = spu * (W.X() + state5(1)*U.X() + state5(2)*V.X()); // a_x
  state7[4] = spu * (W.Y() + state5(1)*U.Y() + state5(2)*V.Y()); // a_y
  state7[5] = spu * (W.Z() + state5(1)*U.Z() + state5(2)*V.Z()); // a_z

  // normalize dir
  double norm = 1. / sqrt(state7[3]*state7[3] + state7[4]*state7[4] + state7[5]*state7[5]);
  for (unsigned int i=3; i<6; ++i) state7[i] *= norm;

  state7[6] = state5(0); // q/p
}


void RKTrackRep::getState5(StateOnPlane* state, const M1x7& state7) const {

  // TODO: Test!

  double spu(1.);

  const TVector3& O(state->getPlane()->getO());
  const TVector3& U(state->getPlane()->getU());
  const TVector3& V(state->getPlane()->getV());
  TVector3 W(state->getPlane()->getNormal());

  TVector3 posShift(state7[0], state7[1], state7[2]);
  posShift -= state->getPlane()->getO();

  // force A to be in normal direction and set spu accordingly
  double AtW = state7[3]*W.X() + state7[4]*W.Y() + state7[5]*W.Z();
  if (AtW < 0) {
    //fDir *= -1.;
    //AtW *= -1.;
    spu = -1.;
  }

  TVectorD& state5 = state->getState();

  state5(0) = state7[6]; // q/p
  state5(1) = (state7[3]*U.X() + state7[4]*U.Y() + state7[5]*U.Z()) / AtW; // u' = (dir * U) / (A * W)
  state5(2) = (state7[3]*V.X() + state7[4]*V.Y() + state7[5]*V.Z()) / AtW; // v' = (dir * V) / (A * W)
  state5(3) = ((state7[0]-O.X())*U.X() +
               (state7[1]-O.Y())*U.Y() +
               (state7[2]-O.Z())*U.Z()); // u = (pos - O) * U
  state5(4) = ((state7[0]-O.X())*V.X() +
               (state7[1]-O.Y())*V.Y() +
               (state7[2]-O.Z())*V.Z()); // v = (pos - O) * V

  setSpu(state, spu);

}



void RKTrackRep::transformPM7(const MeasuredStateOnPlane* state,
                              M7x7& out7x7,
                              TMatrixD* Jac) const {

  // get vectors and aux variables
  const TVector3& U(state->getPlane()->getU());
  const TVector3& V(state->getPlane()->getV());
  TVector3 W(state->getPlane()->getNormal());

  const TVectorD& state5(state->getState());
  double spu(getSpu(state));

  M1x3 pTilde;
  pTilde[0] = spu * (W.X() + state5(1)*U.X() + state5(2)*V.X()); // a_x
  pTilde[1] = spu * (W.Y() + state5(1)*U.Y() + state5(2)*V.Y()); // a_y
  pTilde[2] = spu * (W.Z() + state5(1)*U.Z() + state5(2)*V.Z()); // a_z

  const double pTildeMag = sqrt(pTilde[0]*pTilde[0] + pTilde[1]*pTilde[1] + pTilde[2]*pTilde[2]);
  const double pTildeMag2 = pTildeMag*pTildeMag;

  const double utpTildeOverpTildeMag2 = (U.X()*pTilde[0] + U.Y()*pTilde[1] + U.Z()*pTilde[2]) / pTildeMag2;
  const double vtpTildeOverpTildeMag2 = (V.X()*pTilde[0] + V.Y()*pTilde[1] + V.Z()*pTilde[2]) / pTildeMag2;

  //J_pM matrix is d(x,y,z,ax,ay,az,q/p) / d(q/p,u',v',u,v)   (out is 7x7)

   // d(x,y,z)/d(u)
  fJ_pM_5x7[21] = U.X(); // [3][0]
  fJ_pM_5x7[22] = U.Y(); // [3][1]
  fJ_pM_5x7[23] = U.Z(); // [3][2]
  // d(x,y,z)/d(v)
  fJ_pM_5x7[28] = V.X(); // [4][2]
  fJ_pM_5x7[29] = V.Y(); // [4][2]
  fJ_pM_5x7[30] = V.Z(); // [4][2]
  // d(q/p)/d(q/p)
  fJ_pM_5x7[6] = 1.; // not needed for array matrix multiplication
  // d(ax,ay,az)/d(u')
  double fact = spu / pTildeMag;
  fJ_pM_5x7[10] = fact * ( U.X() - pTilde[0]*utpTildeOverpTildeMag2 ); // [1][3]
  fJ_pM_5x7[11] = fact * ( U.Y() - pTilde[1]*utpTildeOverpTildeMag2 ); // [1][4]
  fJ_pM_5x7[12] = fact * ( U.Z() - pTilde[2]*utpTildeOverpTildeMag2 ); // [1][5]
  // d(ax,ay,az)/d(v')
  fJ_pM_5x7[17] = fact * ( V.X() - pTilde[0]*vtpTildeOverpTildeMag2 ); // [2][3]
  fJ_pM_5x7[18] = fact * ( V.Y() - pTilde[1]*vtpTildeOverpTildeMag2 ); // [2][4]
  fJ_pM_5x7[19] = fact * ( V.Z() - pTilde[2]*vtpTildeOverpTildeMag2 ); // [2][5]


  // since the Jacobian contains a lot of zeros, and the resulting cov has to be symmetric,
  // the multiplication can be done much faster directly on array level
  // out = J_pM^T * in5x5 * J_pM
  const M5x5& in5x5_ = *((M5x5*) state->getCov().GetMatrixArray());
  RKTools::J_pMTxcov5xJ_pM(fJ_pM_5x7, in5x5_, out7x7);

  if (Jac!=NULL){
    Jac->ResizeTo(5,7);
    *Jac = TMatrixD(5,7, &(fJ_pM_5x7[0]));
  }
}


void RKTrackRep::transformPM6(const MeasuredStateOnPlane* state,
                              M6x6& out6x6,
                              TMatrixD* Jac) const {

  // get vectors and aux variables
  const TVector3& U(state->getPlane()->getU());
  const TVector3& V(state->getPlane()->getV());
  TVector3 W(state->getPlane()->getNormal());

  const TVectorD& state5(state->getState());
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

  // d(px,py,pz)/d(q/p)
  double fact = -1. * p / (pTildeMag * qop);
  fJ_pM_5x6[3] = fact * pTilde[0]; // [0][3]
  fJ_pM_5x6[4] = fact * pTilde[1]; // [0][4]
  fJ_pM_5x6[5] = fact * pTilde[2]; // [0][5]
  // d(px,py,pz)/d(u')
  fact = p * spu / pTildeMag;
  fJ_pM_5x6[9]  = fact * ( U.X() - pTilde[0]*utpTildeOverpTildeMag2 ); // [1][3]
  fJ_pM_5x6[10] = fact * ( U.Y() - pTilde[1]*utpTildeOverpTildeMag2 ); // [1][4]
  fJ_pM_5x6[11] = fact * ( U.Z() - pTilde[2]*utpTildeOverpTildeMag2 ); // [1][5]
  // d(px,py,pz)/d(v')
  fJ_pM_5x6[15] = fact * ( V.X() - pTilde[0]*vtpTildeOverpTildeMag2 ); // [2][3]
  fJ_pM_5x6[16] = fact * ( V.Y() - pTilde[1]*vtpTildeOverpTildeMag2 ); // [2][4]
  fJ_pM_5x6[17] = fact * ( V.Z() - pTilde[2]*vtpTildeOverpTildeMag2 ); // [2][5]
  // d(x,y,z)/d(u)
  fJ_pM_5x6[18] = U.X(); // [3][0]
  fJ_pM_5x6[19] = U.Y(); // [3][1]
  fJ_pM_5x6[20] = U.Z(); // [3][2]
  // d(x,y,z)/d(v)
  fJ_pM_5x6[24] = V.X(); // [4][0]
  fJ_pM_5x6[25] = V.Y(); // [4][1]
  fJ_pM_5x6[26] = V.Z(); // [4][2]


  // do the transformation
  // out = J_pM^T * in5x5 * J_pM
  const M5x5& in5x5_ = *((M5x5*) state->getCov().GetMatrixArray());
  RKTools::J_pMTxcov5xJ_pM(fJ_pM_5x6, in5x5_, out6x6);

  if (Jac!=NULL){
    Jac->ResizeTo(5,6);
    *Jac = TMatrixD(5,6, &(fJ_pM_5x6[0]));
  }
}


void RKTrackRep::transformM7P(const M7x7& in7x7,
                              const M1x7& state7,
                              MeasuredStateOnPlane* state, // plane must already be set!
                              TMatrixD* Jac) const {

  // get vectors and aux variables
  const TVector3& U(state->getPlane()->getU());
  const TVector3& V(state->getPlane()->getV());
  TVector3 W(state->getPlane()->getNormal());

  const double AtU = state7[3]*U.X() + state7[4]*U.Y() + state7[5]*U.Z();
  const double AtV = state7[3]*V.X() + state7[4]*V.Y() + state7[5]*V.Z();
  const double AtW = state7[3]*W.X() + state7[4]*W.Y() + state7[5]*W.Z();

  // J_Mp matrix is d(q/p,u',v',u,v) / d(x,y,z,ax,ay,az,q/p)   (in is 7x7)

  // d(u')/d(ax,ay,az)
  double fact = 1./(AtW*AtW);
  fJ_Mp_7x5[16] = fact * (U.X()*AtW - W.X()*AtU); // [3][1]
  fJ_Mp_7x5[21] = fact * (U.Y()*AtW - W.Y()*AtU); // [4][1]
  fJ_Mp_7x5[26] = fact * (U.Z()*AtW - W.Z()*AtU); // [5][1]
  // d(v')/d(ax,ay,az)
  fJ_Mp_7x5[17] = fact * (V.X()*AtW - W.X()*AtV); // [3][2]
  fJ_Mp_7x5[22] = fact * (V.Y()*AtW - W.Y()*AtV); // [4][2]
  fJ_Mp_7x5[27] = fact * (V.Z()*AtW - W.Z()*AtV); // [5][2]
  // d(q/p)/d(q/p)
  fJ_Mp_7x5[30] = 1.; // [6][0]  - not needed for array matrix multiplication
  //d(u)/d(x,y,z)
  fJ_Mp_7x5[3]  = U.X(); // [0][3]
  fJ_Mp_7x5[8]  = U.Y(); // [1][3]
  fJ_Mp_7x5[13] = U.Z(); // [2][3]
  //d(v)/d(x,y,z)
  fJ_Mp_7x5[4]  = V.X(); // [0][4]
  fJ_Mp_7x5[9]  = V.Y(); // [1][4]
  fJ_Mp_7x5[14] = V.Z(); // [2][4]


  // since the Jacobian contains a lot of zeros, and the resulting cov has to be symmetric,
  // the multiplication can be done much faster directly on array level
  // out5x5 = J_Mp^T * in * J_Mp
  M5x5& out5x5_ = *((M5x5*) state->getCov().GetMatrixArray());
  RKTools::J_MpTxcov7xJ_Mp(fJ_Mp_7x5, in7x7, out5x5_);

  if (Jac!=NULL){
    Jac->ResizeTo(7,5);
    *Jac = TMatrixD(7,5, &(fJ_Mp_7x5[0]));
  }
}


void RKTrackRep::transformM6P(const M6x6& in6x6,
                              const M1x7& state7,
                              MeasuredStateOnPlane* state, // plane and charge must already be set!
                              TMatrixD* Jac) const {

  // get vectors and aux variables
  const TVector3& U(state->getPlane()->getU());
  const TVector3& V(state->getPlane()->getV());
  TVector3 W(state->getPlane()->getNormal());

  const double AtU = state7[3]*U.X() + state7[4]*U.Y() + state7[5]*U.Z();
  const double AtV = state7[3]*V.X() + state7[4]*V.Y() + state7[5]*V.Z();
  const double AtW = state7[3]*W.X() + state7[4]*W.Y() + state7[5]*W.Z();

  // J_Mp matrix is d(q/p,u',v',u,v) / d(x,y,z,px,py,pz)       (in is 6x6)

  const double qop = state7[6];
  const double p = getCharge(state)/qop; // momentum

  //d(u)/d(x,y,z)
  fJ_Mp_6x5[3]  = U.X(); // [0][3]
  fJ_Mp_6x5[8]  = U.Y(); // [1][3]
  fJ_Mp_6x5[13] = U.Z(); // [2][3]
  //d(v)/d(x,y,z)
  fJ_Mp_6x5[4]  = V.X(); // [0][4]
  fJ_Mp_6x5[9]  = V.Y(); // [1][4]
  fJ_Mp_6x5[14] = V.Z(); // [2][4]
  // d(q/p)/d(px,py,pz)
  double fact = (-1.) * qop / p;
  fJ_Mp_6x5[15] = fact * state7[3]; // [3][0]
  fJ_Mp_6x5[20] = fact * state7[4]; // [4][0]
  fJ_Mp_6x5[25] = fact * state7[5]; // [5][0]
  // d(u')/d(px,py,pz)
  fact = 1./(p*AtW*AtW);
  fJ_Mp_6x5[16] = fact * (U.X()*AtW - W.X()*AtU); // [3][1]
  fJ_Mp_6x5[21] = fact * (U.Y()*AtW - W.Y()*AtU); // [4][1]
  fJ_Mp_6x5[26] = fact * (U.Z()*AtW - W.Z()*AtU); // [5][1]
  // d(v')/d(px,py,pz)
  fJ_Mp_6x5[17] = fact * (V.X()*AtW - W.X()*AtV); // [3][2]
  fJ_Mp_6x5[22] = fact * (V.Y()*AtW - W.Y()*AtV); // [4][2]
  fJ_Mp_6x5[27] = fact * (V.Z()*AtW - W.Z()*AtV); // [5][2]

  // do the transformation
  // out5x5 = J_Mp^T * in * J_Mp
  M5x5& out5x5_ = *((M5x5*) state->getCov().GetMatrixArray());
  RKTools::J_MpTxcov6xJ_Mp(fJ_Mp_6x5, in6x6, out5x5_);

  if (Jac!=NULL){
    Jac->ResizeTo(6,5);
    *Jac = TMatrixD(6,5, &(fJ_Mp_6x5[0]));;
  }
}



} /* End of namespace genfit */
