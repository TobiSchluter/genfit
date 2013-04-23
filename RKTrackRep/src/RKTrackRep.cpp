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
#include <TDecompSVD.h>
#include <TMath.h>


#define MINSTEP 0.001   // minimum step [cm] for Runge Kutta and iteration to POCA
//#define DEBUG


namespace genfit {


void RKTrackRep::initArrays(){
  memset(fNoise,0x00,7*7*sizeof(double));
  memset(fOldCov,0x00,7*7*sizeof(double));

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
                        M7x7* jacobian,
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

  bool calcJacobian(jacobian != nullptr);

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
  else if (calcJacobian) memcpy(H1, H0, 3*sizeof(double));
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
  else if (calcJacobian) memcpy(H2, H0, 3*sizeof(double));
  A6 = B5*H2[2]-C5*H2[1]; B6 = C5*H2[0]-A5*H2[2]; C6 = A5*H2[1]-B5*H2[0]; // (A5, B5, C5) x H2


  //
  // Derivatives of track parameters
  //
  if(calcJacobian){
    double   dA0(0), dA2(0), dA3(0), dA4(0), dA5(0), dA6(0);
    double   dB0(0), dB2(0), dB3(0), dB4(0), dB5(0), dB6(0);
    double   dC0(0), dC2(0), dC3(0), dC4(0), dC5(0), dC6(0);

    // d(x, y, z)/d(x, y, z) submatrix is unit matrix
    (*jacobian)[0] = 1;  (*jacobian)[8] = 1;  (*jacobian)[16] = 1;
    // d(ax, ay, az)/d(ax, ay, az) submatrix is 0
    // start with d(x, y, z)/d(ax, ay, az)
    for(int i=3*7; i<49; i+=7) {

      if(i==42) {(*jacobian)[i+3]*=state7[6]; (*jacobian)[i+4]*=state7[6]; (*jacobian)[i+5]*=state7[6];}

      //first point
      dA0 = H0[2]*(*jacobian)[i+4]-H0[1]*(*jacobian)[i+5];    // dA0/dp }
      dB0 = H0[0]*(*jacobian)[i+5]-H0[2]*(*jacobian)[i+3];    // dB0/dp  } = dA x H0
      dC0 = H0[1]*(*jacobian)[i+3]-H0[0]*(*jacobian)[i+4];    // dC0/dp }

      if(i==42) {dA0+=A0; dB0+=B0; dC0+=C0;}     // if last row: (dA0, dB0, dC0) := (dA0, dB0, dC0) + (A0, B0, C0)

      dA2 = dA0+(*jacobian)[i+3];        // }
      dB2 = dB0+(*jacobian)[i+4];        //  } = (dA0, dB0, dC0) + dA
      dC2 = dC0+(*jacobian)[i+5];        // }

      //second point
      dA3 = (*jacobian)[i+3]+dB2*H1[2]-dC2*H1[1];    // dA3/dp }
      dB3 = (*jacobian)[i+4]+dC2*H1[0]-dA2*H1[2];    // dB3/dp  } = dA + (dA2, dB2, dC2) x H1
      dC3 = (*jacobian)[i+5]+dA2*H1[1]-dB2*H1[0];    // dC3/dp }

      if(i==42) {dA3+=A3-A[0]; dB3+=B3-A[1]; dC3+=C3-A[2];} // if last row: (dA3, dB3, dC3) := (dA3, dB3, dC3) + (A3, B3, C3) - (ax, ay, az)

      dA4 = (*jacobian)[i+3]+dB3*H1[2]-dC3*H1[1];    // dA4/dp }
      dB4 = (*jacobian)[i+4]+dC3*H1[0]-dA3*H1[2];    // dB4/dp  } = dA + (dA3, dB3, dC3) x H1
      dC4 = (*jacobian)[i+5]+dA3*H1[1]-dB3*H1[0];    // dC4/dp }

      if(i==42) {dA4+=A4-A[0]; dB4+=B4-A[1]; dC4+=C4-A[2];} // if last row: (dA4, dB4, dC4) := (dA4, dB4, dC4) + (A4, B4, C4) - (ax, ay, az)

      //last point
      dA5 = dA4+dA4-(*jacobian)[i+3];      // }
      dB5 = dB4+dB4-(*jacobian)[i+4];      //  } =  2*(dA4, dB4, dC4) - dA
      dC5 = dC4+dC4-(*jacobian)[i+5];      // }

      dA6 = dB5*H2[2]-dC5*H2[1];      // dA6/dp }
      dB6 = dC5*H2[0]-dA5*H2[2];      // dB6/dp  } = (dA5, dB5, dC5) x H2
      dC6 = dA5*H2[1]-dB5*H2[0];      // dC6/dp }

      if(i==42) {dA6+=A6; dB6+=B6; dC6+=C6;}     // if last row: (dA6, dB6, dC6) := (dA6, dB6, dC6) + (A6, B6, C6)

      if(i==42) {
        (*jacobian)[i]   += (dA2+dA3+dA4)*S3/state7[6];  (*jacobian)[i+3] = (dA0+dA3+dA3+dA5+dA6)*P3/state7[6]; // dR := dR + S3*[(dA2, dB2, dC2) +   (dA3, dB3, dC3) + (dA4, dB4, dC4)]
        (*jacobian)[i+1] += (dB2+dB3+dB4)*S3/state7[6];  (*jacobian)[i+4] = (dB0+dB3+dB3+dB5+dB6)*P3/state7[6]; // dA :=     1/3*[(dA0, dB0, dC0) + 2*(dA3, dB3, dC3) + (dA5, dB5, dC5) + (dA6, dB6, dC6)]
        (*jacobian)[i+2] += (dC2+dC3+dC4)*S3/state7[6];  (*jacobian)[i+5] = (dC0+dC3+dC3+dC5+dC6)*P3/state7[6];
      }
      else {
        (*jacobian)[i]   += (dA2+dA3+dA4)*S3;  (*jacobian)[i+3] = (dA0+dA3+dA3+dA5+dA6)*P3; // dR := dR + S3*[(dA2, dB2, dC2) +   (dA3, dB3, dC3) + (dA4, dB4, dC4)]
        (*jacobian)[i+1] += (dB2+dB3+dB4)*S3;  (*jacobian)[i+4] = (dB0+dB3+dB3+dB5+dB6)*P3; // dA :=     1/3*[(dA0, dB0, dC0) + 2*(dA3, dB3, dC3) + (dA5, dB5, dC5) + (dA6, dB6, dC6)]
        (*jacobian)[i+2] += (dC2+dC3+dC4)*S3;  (*jacobian)[i+5] = (dC0+dC3+dC3+dC5+dC6)*P3;
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

  if (Jac!=nullptr){
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

  if (Jac!=nullptr){
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

  if (Jac!=nullptr){
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

  if (Jac!=nullptr){
    Jac->ResizeTo(6,5);
    *Jac = TMatrixD(6,5, &(fJ_Mp_6x5[0]));;
  }
}


//
// Runge-Kutta method for tracking a particles through a magnetic field.
// Uses Nystroem algorithm (See Handbook Nat. Bur. of Standards, procedure 25.5.20)
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
bool RKTrackRep::RKutta (const DetPlane& plane,
                         double charge,
                         M1x7& state7,
                         M7x7* jacobian,
                         double& coveredDistance,
                         bool& checkJacProj,
                         TMatrixD& noiseProjection,
                         bool onlyOneStep,
                         double maxStep) const {

  // limits, check-values, etc. Can be tuned!
  static const double Wmax   = 3000.;           // max. way allowed [cm]
  static const double AngleMax = 6.3;           // max. total angle change of momentum. Prevents extrapolating a curler round and round if no active plane is found.
  static const double Pmin   = 4.E-3;           // minimum momentum for propagation [GeV]
  static const unsigned int maxNumIt = 1000;    // maximum number of iterations in main loop
  // Aux parameters
  M1x3&   R           = *((M1x3*) &state7[0]);  // Start coordinates  [cm]  (x,  y,  z)
  M1x3&   A           = *((M1x3*) &state7[3]);  // Start directions         (ax, ay, az);   ax^2+ay^2+az^2=1
  M1x3    SA          = {0.,0.,0.};             // Start directions derivatives dA/S
  double  Way         = 0.;                     // Sum of absolute values of all extrapolation steps [cm]
  bool    atPlane = false;                      // stepper thinks that the plane will be reached in that step -> linear extrapolation and projection of jacobian
  bool    momLossExceeded = false;              // stepper had to limit stepsize due to momentum loss -> no next RKutta loop, no linear extrapolation
  double  momentum   = fabs(charge/state7[6]);// momentum [GeV]
  double  relMomLoss = 0;                      // relative momentum loss in RKutta
  double  deltaAngle = 0.;                     // total angle by which the momentum has changed during extrapolation
  double  An(0), S(0), Sl(0), CBA(0);
  M1x4    SU = {0.,0.,0.,0.};


  #ifdef DEBUG
    std::cout << "RKTrackRep::RKutta \n";
    std::cout << "position: "; fPos.Print();
    std::cout << "direction: "; fDir.Print();
    std::cout << "momentum: " << momentum << " GeV\n";
    std::cout << "destination: "; plane.Print();
  #endif

  checkJacProj = false;

  // check momentum
  if(momentum < Pmin){
    std::ostringstream sstream;
    sstream << "RKTrackRep::RKutta ==> momentum too low: " << momentum*1000. << " MeV";
    Exception exc(sstream.str(),__LINE__,__FILE__);
    exc.setFatal();
    throw exc;
  }


  // make SU vector point away from origin
  const TVector3 W(plane.getNormal());
  if (W*plane.getO() > 0) {
    SU[0] =     W.X();
    SU[1] =     W.Y();
    SU[2] =     W.Z();
  }
  else {
    SU[0] = -1.*W.X();
    SU[1] = -1.*W.Y();
    SU[2] = -1.*W.Z();
  }

  SU[3] = plane.distance(0., 0., 0.);

  unsigned int counter(0);

  // Step estimation (signed)
  S = estimateStep(state7, SU, plane, charge, relMomLoss, momLossExceeded, atPlane, maxStep);
  if (fabs(S) < 0.001*MINSTEP) {
    #ifdef DEBUG
      std::cout << " RKutta - step too small -> break \n";
    #endif
    ++counter; // skip the main loop, go to linear extrapolation step (will be skipped) and just project jacobian
  }

  //
  // Main loop of Runge-Kutta method
  //
  while (fabs(S) >= MINSTEP || counter == 0) {

    if(++counter > maxNumIt){
      Exception exc("RKTrackRep::RKutta ==> maximum number of iterations exceeded",__LINE__,__FILE__);
      exc.setFatal();
      throw exc;
    }

    #ifdef DEBUG
      std::cout << "------ RKutta main loop nr. " << counter-1 << " ------\n";
    #endif

    // update paths
    coveredDistance += S;       // add stepsize to way (signed)
    Way  += fabs(S);

    // check way limit
    if(Way > Wmax){
      std::ostringstream sstream;
      sstream << "RKTrackRep::RKutta ==> Total extrapolation length is longer than length limit : " << Way << " cm !";
      Exception exc(sstream.str(),__LINE__,__FILE__);
      exc.setFatal();
      throw exc;
    }

    M1x3 ABefore{A[0], A[1], A[2]};
    RKPropagate(state7, jacobian, SA, S); // the actual Runge Kutta propagation
    deltaAngle += acos(ABefore[0]*A[0] + ABefore[1]*A[1] + ABefore[2]*A[2]);

    if (onlyOneStep) return(true);

    // if stepsize has been limited by material, break the loop and return. No linear extrapolation!
    if (momLossExceeded) {
      #ifdef DEBUG
        std::cout<<" momLossExceeded -> return(true); \n";
      #endif
      return(true);
    }

    // estimate Step for next loop or linear extrapolation
    Sl = S; // last S used
    S = estimateStep(state7, SU, plane, charge, relMomLoss, momLossExceeded, atPlane, maxStep);

    if (atPlane && fabs(S) < MINSTEP) {
      #ifdef DEBUG
        std::cout<<" (atPlane && fabs(S) < MINSTEP) -> break; \n";
      #endif
      break; // else if at plane: do linear extrapolation
    }
    if (momLossExceeded && fabs(S) < MINSTEP) {
      #ifdef DEBUG
        std::cout<<" (momLossExceeded && fabs(S) < MINSTEP) -> return(true); \n";
      #endif
      return(true); // no linear extrapolation!
    }

    // check if total angle is bigger than AngleMax. Can happen if a curler should be fitted and it does not hit the active area of the next plane.
    if (fabs(deltaAngle) > AngleMax){
      std::ostringstream sstream;
      sstream << "RKTrackRep::RKutta ==> Do not get to an active plane! Already extrapolated " << deltaAngle * 180 / TMath::Pi() << "°.";
      Exception exc(sstream.str(),__LINE__,__FILE__);
      exc.setFatal();
      throw exc;
    }

    // check if we went back and forth multiple times -> we don't come closer to the plane!
    if (counter > 3){
      if (S                                       *materials_[counter-1].getSegmentLength() < 0 &&
          materials_[counter-1].getSegmentLength()*materials_[counter-2].getSegmentLength() < 0 &&
          materials_[counter-2].getSegmentLength()*materials_[counter-3].getSegmentLength() < 0){
        Exception exc("RKTrackRep::RKutta ==> Do not get closer to plane!",__LINE__,__FILE__);
        exc.setFatal();
        throw exc;
      }
    }

  } //end of main loop


  //
  // linear extrapolation to surface
  //
  if (atPlane) {

    if (fabs(Sl) > 0.001*MINSTEP){
      #ifdef DEBUG
        std::cout << " RKutta - linear extrapolation to surface\n";
      #endif
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
    }
#ifdef DEBUG
    else {
      std::cout << " RKutta - last stepsize too small -> can't do linear extrapolation! \n";
    }
#endif

    //
    // Project Jacobian of extrapolation onto destination plane
    //
    if(jacobian != nullptr){
      if (checkJacProj && materials_.size()>0){
        Exception exc("RKTrackRep::Extrap ==> covariance is projected onto destination plane again",__LINE__,__FILE__);
        throw exc;
      }
      //save old jacobian
      double * covAsPtr = (double*)jacobian;

      noiseProjection.SetMatrixArray(covAsPtr);

      //std::cerr << "The current Jac is:" << std::endl;
      //RKTools::printDim(covAsPtr,7,7);
      //std::cerr << "This was filled into noiseProjection so it must be the same:" << std::endl;
//       noiseProjection.Print();
      checkJacProj = true;
      #ifdef DEBUG
        std::cout << "  Project Jacobian of extrapolation onto destination plane\n";
      #endif
      An = A[0]*SU[0] + A[1]*SU[1] + A[2]*SU[2];
      fabs(An) > 1.E-7 ? An=1./An : An = 0; // 1/A_normal
      double norm;
      for(int i=0; i<49; i+=7) {
        norm = ((*jacobian)[i]*SU[0] + (*jacobian)[i+1]*SU[1] + (*jacobian)[i+2]*SU[2])*An;  // dR_normal / A_normal
        (*jacobian)[i]   -= norm*A [0];   (*jacobian)[i+1] -= norm*A [1];   (*jacobian)[i+2] -= norm*A [2];
        (*jacobian)[i+3] -= norm*SA[0];   (*jacobian)[i+4] -= norm*SA[1];   (*jacobian)[i+5] -= norm*SA[2];
      }
      TMatrixD projectedJac(7,7);

      projectedJac.SetMatrixArray(covAsPtr);

//       std::cerr << "The projected Jac is filled into projectedJac so we have:" << std::endl;
//       projectedJac.Print();
      //Tools::invert(noiseProjection);

      TDecompSVD invertAlgo(noiseProjection);

  bool status = invertAlgo.Invert(noiseProjection);
//     std::cerr << "The inverse of the unprojected jac is:" << std::endl;
//       noiseProjection.Print();
      noiseProjection = projectedJac * noiseProjection;
//       std::cerr << "And finaly in RKutta. The noise projection matrix is:" << std::endl;
      //noiseProjection = projectedJac * noiseProjection;
//       noiseProjection.Print();
    }

    coveredDistance += S;
    Way  += fabs(S);
  } // end of linear extrapolation to surface

  return(true);

}


double RKTrackRep::estimateStep(const M1x7& state7,
                                const M1x4& SU,
                                const DetPlane& plane,
                                const double& charge,
                                double& relMomLoss,
                                bool& momLossExceeded,
                                bool& atPlane,
                                double maxStepArg) const {

  static const double Smax = 25.; // max. step allowed [cm]
  double Step(0);
  bool improveEstimation (true);

  momLossExceeded = false;
  atPlane = false;

  #ifdef DEBUG
    std::cout << " RKTrackRep::estimateStep \n";
    std::cout << "  position: "; TVector3(state7[0], state7[1], state7[2]).Print();
    std::cout << "  direction: "; TVector3(state7[3], state7[4], state7[5]).Print();
  #endif

  // check and limit maxStepArg argument
  if (maxStepArg < 0) maxStepArg *= -1;
  if (maxStepArg > Smax) maxStepArg = Smax;

  // calculate distance to surface
  double Dist = SU[3] - (state7[0]*SU[0] + state7[1]*SU[1] + state7[2]*SU[2]);  // Distance between start coordinates and surface
  double An = state7[3]*SU[0] + state7[4]*SU[1] + state7[5]*SU[2];              // An = dir * N;  component of dir normal to surface

  if (fabs(An) > 1.E-10) Step = Dist/An;
  else {
    Step = Dist*1.E10;
    if (An<0) Step *= -1.;
  }

  // see if dir points towards surface (1) or not (-1)
  double StepSign(1);
  if (Step<0) StepSign = -1;

  #ifdef DEBUG
    std::cout << "  Distance to plane: " << Dist << "\n";
    std::cout << "  guess for Step: " << Step << "\n";
    if (StepSign>0) std::cout << "  Direction is  pointing towards surface.\n";
    else  std::cout << "  Direction is pointing away from surface.\n";
  #endif


  //
  // Limit maxStepArg according to curvature and magnetic field inhomogenities
  //
  // limit maxStepArg to the straight line distance to plane times some margin
  // (only if plane is not finite, otherwise maxStepArg only limited by the field & curvature is needed to set Step = propDir_*maxStepArg)
  static const double margin(1.5);
  if (!plane.isFinite() && maxStepArg > fabs(margin*Step)) maxStepArg = fabs(margin*Step);

  while (maxStepArg > MINSTEP) {
    M1x7 state7_temp;
    memcpy(state7_temp, state7, sizeof(state7));
    M1x3 SA;

    double q = RKPropagate(state7_temp, nullptr, SA, StepSign*maxStepArg, true);
#ifdef DEBUG
    std::cerr << "  maxStepArg = " << maxStepArg << "; q = " << q  << " \n";
#endif
    maxStepArg *= q * 0.95;

    if (q >= 1) break;
  }

#ifdef DEBUG
  std::cerr << "  limit from curvature and magnetic field inhomogenities: " << maxStepArg << "\n";
#endif

  // limit Step to maxStepArg
  if (fabs(Step) > maxStepArg) {
    Step = StepSign*maxStepArg;
    improveEstimation = false;
    #ifdef DEBUG
      std::cerr << "  improveEstimation = false;\n";
    #endif
  }


  //
  // Select direction
  //
  // auto select
  if (propDir_ == 0 || !plane.isFinite()){
    #ifdef DEBUG
      std::cerr << "  auto select direction";
      if (!plane.isFinite()) std::cerr << ", plane is not finite";
      std::cerr << ".\n";
    #endif
  }
  // see if straight line approximation is ok
  else if ( fabs(Step) < 0.2*maxStepArg ){
    #ifdef DEBUG
      std::cerr << "  straight line approximation is fine.\n";
    #endif

    // if direction is pointing to active part of surface
    if( plane.isInActive(state7[0], state7[1], state7[2],  state7[3], state7[4], state7[5]) ) {
      #ifdef DEBUG
        std::cerr << "  direction is pointing to active part of surface. \n";
      #endif
    }
    // if we are near the plane, but not pointing to the active area, make a big step!
    else {
      Step = propDir_*maxStepArg;
      improveEstimation = false;
      #ifdef DEBUG
        std::cerr << "  we are near the plane, but not pointing to the active area. make a big step! \n";
      #endif
    }
  }
  // propDir_ decides!
  else {
    if (Step * propDir_ < 0){
      Step = propDir_*maxStepArg;
      improveEstimation = false;
      #ifdef DEBUG
        std::cerr << "  invert Step according to propDir_ and set Step to propDir_*maxStepArg. \n";
      #endif
    }
  }

  #ifdef DEBUG
    std::cerr << "  guess for Step (signed): " << Step << "\n";
  #endif

  // re-check sign of Step
  if (Step>=0) StepSign = 1;
  else StepSign = -1;



  // call stepper and reduce stepsize if step not too small
  if (/*!fNoMaterial*/ true){

    if(fabs(Step) > MINSTEP){ // only call stepper if step estimation big enough

      M1x7 state7_temp;
      memcpy(state7_temp, state7, sizeof(state7));
      for (unsigned int i=3; i<6; ++i)
        state7_temp[i] *= StepSign;

      materials_.push_back(MaterialProperties());
      double StepMat = MaterialEffects::getInstance()->stepper(this,
                                                               state7_temp,
                                                               maxStepArg,
                                                               charge/state7[6], // |p|
                                                               relMomLoss,
                                                               pdgCode_,
                                                               materials_.back(),
                                                               true);
      if (fabs(Step) > StepMat) {
        Step = StepSign*StepMat;
        momLossExceeded = true;
      }

      #ifdef DEBUG
        std::cerr << "  limit from stepper: " << Step << "\n";
      #endif
    }
  }


  if (!momLossExceeded && improveEstimation){
    atPlane = true;
   #ifdef DEBUG
     std::cerr << "  At plane. " << Step << "\n";
   #endif

    // improve step estimation to surface according to curvature
    if (fabs(Step) > 0.1*maxStepArg && fabs(Step) > MINSTEP){

      M1x7 state7_temp;
      memcpy(state7_temp, state7, sizeof(state7));
      M1x3 SA;

      RKPropagate(state7_temp, nullptr, SA, Step, true);

      // calculate distance to surface
      Dist = SU[3] - (state7_temp[0] * SU[0] +
                      state7_temp[1] * SU[1] +
                      state7_temp[2] * SU[2]); // Distance between position and surface

      An = state7_temp[3] * SU[0] +
           state7_temp[4] * SU[1] +
           state7_temp[5] * SU[2];    // An = dir * N;  component of dir normal to surface

      Step += Dist/An;

      #ifdef DEBUG
        std::cout << "  Improved step estimation taking curvature into account: " << Step << "\n";
      #endif
    }

  }

  materials_.back().setSegmentLength(Step);

  #ifdef DEBUG
    std::cout << "  --> Step to be used: " << Step << "\n";
  #endif

  return Step;

}


TVector3 RKTrackRep::poca2Line(const TVector3& extr1,const TVector3& extr2,const TVector3& point) const {

  TVector3 pocaOnLine(extr2);
  pocaOnLine -= extr1; // wireDir

  if(pocaOnLine.Mag()<1.E-8){
    Exception exc("RKTrackRep::poca2Line ==> try to find POCA between line and point, but the line is really just a point",__LINE__,__FILE__);
    throw exc;
  }

  double t = 1./(pocaOnLine.Mag2()) * ((point*pocaOnLine) + extr1.Mag2() - (extr1*extr2));
  pocaOnLine *= t;
  pocaOnLine += extr1;
  return pocaOnLine; // = extr1 + t*wireDir

}


double RKTrackRep::Extrap(const DetPlane& plane, double charge, M1x7& state7, M7x7* cov, bool onlyOneStep, double maxStep) const {

  static const unsigned int maxNumIt(500);
  unsigned int numIt(0);

  const bool calcCov(cov!=nullptr);
  double coveredDistance(0.);
  double sumDistance(0.);
  TMatrixD noiseProjection(7,7);

  while(true){

    #ifdef DEBUG
      std::cout << "\n============ RKTrackRep::Extrap loop nr. " << numIt << " ============\n";
    #endif

    if(numIt++ > maxNumIt){
      Exception exc("RKTrackRep::Extrap ==> maximum number of iterations exceeded",__LINE__,__FILE__);
      exc.setFatal();
      throw exc;
    }

    // initialize cov with unit matrix
    if(calcCov){
      memcpy(fOldCov, cov, 7*7*sizeof(double));
      memset(cov,0x00,49*sizeof(double));
      for(int i=0; i<7; ++i) (*cov)[8*i] = 1.;
    }

    // propagation
    bool checkJacProj = true;

    if( ! this->RKutta(plane, charge, state7, cov, coveredDistance, checkJacProj, noiseProjection, onlyOneStep, maxStep) ) {
      Exception exc("RKTrackRep::Extrap ==> Runge Kutta propagation failed",__LINE__,__FILE__);
      exc.setFatal();
      throw exc;
    }

    #ifdef DEBUG
      std::cout<<"Original points \n";
      for (unsigned int i=0; i<points.size(); ++i){
        points[i].Print();
      }
      std::cout<<"\n";
    #endif

    sumDistance+=coveredDistance;

    // filter points // TODO: test!!!
    if (/*!fNoMaterial*/ true) { // points are only filled if mat fx are on
      if(materials_.size() > 2){ // check if there are at least three points
        for (unsigned int i=materials_.size()-1; i>0; --i){
          if (materials_[i] == materials_[i-1] &&
              (fabs(materials_[i].getSegmentLength()) < MINSTEP ||
               fabs(materials_[i-1].getSegmentLength()) < MINSTEP) ){
            materials_[i-1].addToSegmentLength(materials_[i].getSegmentLength());
            materials_.erase(materials_.begin()+i);
          }
        }
      }
      #ifdef DEBUG
        std::cout<<"Filtered materials_ \n";
        for (unsigned int i=0; i<materials_.size(); ++i){
          materials_[i].Print();
        }
        std::cout<<"\n";
      #endif
    }


    if(calcCov) memset(fNoise,0x00,7*7*sizeof(double)); // set fNoise to 0


    // call MatFX
    unsigned int nPoints(materials_.size());
    if (/*!fNoMaterial*/ true && nPoints>0){
      // momLoss has a sign - negative loss means momentum gain
      double momLoss = MaterialEffects::getInstance()->effects(materials_,
                                                               fabs(charge/state7[6]), // momentum
                                                               pdgCode_,
                                                               fNoise,
                                                               (double *)cov);

      #ifdef DEBUG
        std::cout << "momLoss: " << momLoss << " GeV; relative: " << momLoss/fabs(charge/state7[6]) << "\n";
      #endif

      // do momLoss only for defined 1/momentum .ne.0
      if(fabs(state7[6])>1.E-10) state7[6] = charge/(fabs(charge/state7[6])-momLoss);
    } // finished MatFX

    if(calcCov){ // propagate cov and add noise
      // numerical check:
      for(unsigned int i=0; i<7*7; ++i){
        if(fabs((*cov)[i]) > 1.E100){
          Exception exc("RKTrackRep::Extrap ==> covariance matrix exceeds numerical limits",__LINE__,__FILE__);
          exc.setFatal();
          throw exc;
        }
      }

      // cov = Jac^T * oldCov * Jac;
      // last column of jac is [0,0,0,0,0,0,1]
      // cov is symmetric
      RKTools::J_MMTxcov7xJ_MM(*cov, fOldCov);
      memcpy(cov, fOldCov, 7*7*sizeof(double));
      if( checkJacProj == true ){

  // XXX std::cerr << "noise in 7D before it gets added" << std::endl;
  // XXX RKTools::printDim((double*)fNoise,7,7);
  // XXX std::cerr << "noise in 5D before it gets added" << std::endl;
  // XXX TMatrixDSym noise5D;
  // XXX transformM7P( fNoise,noise5D,plane,state7);
  // XXX noise5D.Print();
        //project the noise onto the measurment plane
         //std::cerr << "the current noise is " << std::endl;
         //RKTools::printDim(fNoise,7,7);
  TMatrixDSym projectedNoise(7);
  projectedNoise.SetMatrixArray(fNoise);
  //std::cerr << "projectedNoise is filled with the current noise: " << std::endl;
  //projectedNoise.Print();
  // XXX std::cerr << "projection matrix for noise:" << std::endl;
  // XXX RKTools::printDim(noiseProjection.GetMatrixArray(),7,7);
  // XXX TMatrixD P2 = noiseProjection * noiseProjection;
  // XXX std::cerr << "P^2:" << std::endl;
  // XXX RKTools::printDim(P2.GetMatrixArray(),7,7);
  projectedNoise.SimilarityT(noiseProjection);
  //std::cerr << "and finally the projecte noise is as root matrix: " << std::endl;
  //projectedNoise.Print();

  double* projectedNoisePtr = projectedNoise.GetMatrixArray();
  // XXX std::cerr << "projected noise in 7D before it gets added" << std::endl;
  // XXX RKTools::printDim(projectedNoisePtr,7,7);
  // XXX std::cerr << "projected noise in 5D before it gets added" << std::endl;
  //TMatrixDSym noise5D;
  // XXX transformM7P( ( M7x7&) (*projectedNoisePtr),noise5D,plane,state7);
  // XXX noise5D.Print();
  // add noise to cov
  // XXX std::cerr << "the covariance in 7D before the noiese gets added" << std::endl;
  // XXX RKTools::printDim((double*)cov,7,7);
  // XXX std::cerr << "the covariance in 5D before the noiese gets added" << std::endl;
  // XXX TMatrixDSym cov5D;
  // XXX transformM7P( *cov,cov5D,plane,state7);
  // XXX cov5D.Print();
//  for (int i=0; i<7*7; ++i) (*cov)[i] += projectedNoisePtr[i];

  // XXX std::cerr << "cov + noise in 7D" << std::endl;
  // XXX RKTools::printDim((double*)cov,7,7);
  // XXX std::cerr << "cov + noise in 5D" << std::endl;
  // XXX transformM7P( *cov,cov5D,plane,state7);
  // XXX cov5D.Print();
  for (int i=0; i<7*7; ++i) (*cov)[i] += fNoise[i];
      } else {
  // XXX std::cerr << "noise in 7D before it gets added" << std::endl;
  // XXX RKTools::printDim(fNoise,7,7);
  // XXX std::cerr << "cov in 7D before noise gets added" << std::endl;
  // XXX RKTools::printDim((double*)cov,7,7);
  for (int i=0; i<7*7; ++i) (*cov)[i] += fNoise[i];
  // XXX std::cerr << "cov + noise 7D" << std::endl;
  // XXX RKTools::printDim((double*)cov,7,7);
      }


     // std::cerr << "Noise was added to cov in RKTrackRep" << std::endl;
    } // finished propagate cov and add noise

    if (onlyOneStep) break;

    //we arrived at the destination plane, if we point to the active area of the plane (if it is finite), and the distance is below threshold
    if( plane.distance(state7[0], state7[1], state7[2]) < MINSTEP) {
      #ifdef DEBUG
        std::cerr << "arrived at plane with a distance of  " << plane.distance((state7[0], state7[1], state7[2]) << " cm left. ";
      #endif
      if (plane.isInActive(state7[0], state7[1], state7[2],  state7[3], state7[4], state7[5])) {
        #ifdef DEBUG
          std::cerr << "In active area of plane. \n";
        #endif
        // check if Jacobian has been projected onto plane; Otherwise make another iteration
        if (calcCov && !checkJacProj && nPoints>0){
          #ifdef DEBUG
            std::cout << "Jacobian was not projected onto destination plane -> one more iteration. \n";
          #endif
          continue;
        }
        break;
      }
      #ifdef DEBUG
      else {
       std::cerr << "NOT in active area of plane. \n";
      }
      #endif
    }

  }

  return sumDistance;
}



} /* End of namespace genfit */
