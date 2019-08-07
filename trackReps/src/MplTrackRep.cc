/* Copyright 2019, Belle II Collaboration
   Authors: Dmitrii Neverov

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

#include "MplTrackRep.h"

#include <FieldManager.h>
#include <TDatabasePDG.h>
#include <MeasurementOnPlane.h>
#include <Exception.h>

#include <math.h>

#include <TBuffer.h>

using namespace genfit;

MplTrackRep::MplTrackRep(int pdgCode, float magCharge, char propDir) :
  RKTrackRep(pdgCode, propDir),
  m_magCharge(magCharge),
  m_mass(TDatabasePDG::Instance()->GetParticle(pdgCode)->Mass()) // We could ofc use AbsTrackRep::getMass(state) but we have no state here to call on
{
}

MplTrackRep::~MplTrackRep()
{
}

double MplTrackRep::getCharge(const StateOnPlane& state) const {

  if (dynamic_cast<const MeasurementOnPlane*>(&state) != nullptr) {
    Exception exc("RKTrackRep::getCharge - cannot get charge from MeasurementOnPlane",__LINE__,__FILE__);
    exc.setFatal();
    throw exc;
  }

  double pdgCharge( m_magCharge * (this->getPDGCharge() > 0 ? 1.0 : -1.0));

  // return pdgCharge with sign of q/p
  if (state.getState()(0) * pdgCharge < 0)
    return -pdgCharge;
  else
    return pdgCharge;
}

double MplTrackRep::RKPropagate(M1x7& state7,
                                M7x7* jacobianT,
                                M1x3& SA,
                                double S,
                                bool varField,
                                bool calcOnlyLastRowOfJ) const
{
  // The algorithm is
  //  E Lund et al 2009 JINST 4 P04001 doi:10.1088/1748-0221/4/04/P04001
  //  "Track parameter propagation through the application of a new adaptive Runge-Kutta-Nyström method in the ATLAS experiment"
  //  http://inspirehep.net/search?ln=en&ln=en&p=10.1088/1748-0221/4/04/P04001&of=hb&action_search=Search&sf=earliestdate&so=d&rm=&rg=25&sc=0
  // where the transport of the Jacobian is described in
  //   L. Bugge, J. Myrheim  Nucl.Instrum.Meth. 160 (1979) 43-48
  //   "A Fast Runge-kutta Method For Fitting Tracks In A Magnetic Field"
  //   http://inspirehep.net/record/145692
  // and
  //   L. Bugge, J. Myrheim  Nucl.Instrum.Meth. 179 (1981) 365-381
  //   "Tracking And Track Fitting"
  //   http://inspirehep.net/record/160548

  // important fixed numbers
  static const double EC  ( 0.000149896229 );  // c/(2*10^12) resp. c/2Tera FIXME this 1/2 here is super sneaky
  static const double P3  ( 1./3. );           // 1/3
  static const double DLT ( .0002 );           // max. deviation for approximation-quality test
  double sign = state7[6] > 0 ? 1.0 : -1.0;
  // Aux parameters
  M1x3&   R           = *((M1x3*) &state7[0]);       // Start coordinates  [cm]  (x,  y,  z)
  M1x3&   A           = *((M1x3*) &state7[3]);       // Start directions         (ax, ay, az);   ax^2+ay^2+az^2=1
  double  S3(0), S4(0), PS2(0);
  M1x3     H0 = {{0.,0.,0.}}, H1 = {{0.,0.,0.}}, H2 = {{0.,0.,0.}};
  M1x3     r = {{0.,0.,0.}};
  // Variables for Runge Kutta solver
  double   A0(0), A1(0), A2(0), A3(0), A4(0), A5(0), A6(0);
  double   B0(0), B1(0), B2(0), B3(0), B4(0), B5(0), B6(0);
  double   C0(0), C1(0), C2(0), C3(0), C4(0), C5(0), C6(0);
  // Additional variables for momentum evolution FIXME these are all cryptic in accordance with the rest of the code around
  double   D0(0), D1(0), D2(0), D3(0), D4(0), D5(0);
  double   F0(0), F1(0), F2(0), F3(0);
  double   AH0(0), AH1(0), AH2(0), AH3(0);

  //
  // Runge Kutta Extrapolation
  //
  S3 = P3*S;
  S4 = 0.25*S;
  PS2 = m_magCharge * EC * S * sign;

 // First point
  r[0] = R[0];           r[1] = R[1];           r[2]=R[2];
  FieldManager::getInstance()->getFieldVal(r[0], r[1], r[2], H0[0], H0[1], H0[2]);       // magnetic field in 10^-1 T = kGauss
  H0[0] *= PS2; H0[1] *= PS2; H0[2] *= PS2;     // H0 is PS2*(Hx, Hy, Hz) @ R0; effectively this is h/2 * Force
  D0 = fabs(m_magCharge/state7[6]); // p_n
  F0 = std::sqrt(m_mass * m_mass + D0 * D0) / (D0 * D0); // E / p^2
  AH0 = A[0]*H0[0] + A[1]*H0[1] + A[2]*H0[2]; // A dot Force

  A0 = F0 * (H0[0] - A[0] * AH0); B0 = F0 * (H0[1] - A[1] * AH0); C0 = F0 * (H0[2] - A[2] * AH0); // h/2 * k_1
  A2 = A[0]+A0              ; B2 = A[1]+B0              ; C2 = A[2]+C0              ; // r'_n + h/2 * k_1
  A1 = A2+A[0]              ; B1 = B2+A[1]              ; C1 = C2+A[2]              ; // 2*r'_n + h/2 * k_1

  // Second point
  if (varField) {
    r[0] += A1*S4;         r[1] += B1*S4;         r[2] += C1*S4;
    FieldManager::getInstance()->getFieldVal(r[0], r[1], r[2], H1[0], H1[1], H1[2]);
    H1[0] *= PS2; H1[1] *= PS2; H1[2] *= PS2; // H1 is PS2*(Hx, Hy, Hz) @ (x, y, z) + 0.25*S * [(A0, B0, C0) + 2*(ax, ay, az)]
  }
  else { H1 = H0; };
  D1 = D0 + F0 * D0 * AH0; // p_n + h/2 * l_1
  F1 = std::sqrt(m_mass * m_mass + D1 * D1) / (D1 * D1); // E / p^2
  AH1 = A2*H1[0] + B2*H1[1] + C2*H1[2]; // A dot Force

  A3 = A[0] + F1*(H1[0] - A2*AH1); B3 = A[1] + F1*(H1[1] - B2*AH1); C3 = A[2] + F1*(H1[2] - C2*AH1); // r'_n + h/2 * k_2
  D2 = D0 + F1 * D1 * AH1; // p_n + h/2 * l_2
  F2 = std::sqrt(m_mass * m_mass + D2 * D2) / (D2 * D2); // E / p^2
  AH2 = A3*H1[0] + B3*H1[1] + C3*H1[2]; // A dot Force

  A4 = A[0] + F2*(H1[0] - A3*AH2); B4 = A[1] + F2*(H1[1] - B3*AH2); C4 = A[2] + F2*(H1[2] - C3*AH2); // r'_n + h/2 * k_3
  A5 = A4-A[0]+A4            ; B5 = B4-A[1]+B4            ; C5 = C4-A[2]+C4            ; //    r'_n + h * k_3
  D3 = D0 + 2.0 * F2 * D2 * AH2; // p_n + h * l_3
  F3 = std::sqrt(m_mass * m_mass + D3 * D3) / (D3 * D3); // E / p^2
  AH3 = A4*H1[0] + B4*H1[1] + C4*H1[2]; // A dot Force

  // Last point
  if (varField) {
    r[0]=R[0]+S*A4;         r[1]=R[1]+S*B4;         r[2]=R[2]+S*C4;  //setup.Field(r,H2);
    FieldManager::getInstance()->getFieldVal(r[0], r[1], r[2], H2[0], H2[1], H2[2]);
    H2[0] *= PS2; H2[1] *= PS2; H2[2] *= PS2; // H2 is PS2*(Hx, Hy, Hz) @ (x, y, z) + 0.25*S * (A4, B4, C4)
  }
  else { H2 = H0; };
  A6 = F3 * (H2[0] - A5*AH3); B6 = F3 * (H2[1] - B5*AH3); C6 = F3 * (H2[2] - C5*AH3); // h/2 * k_4
  D4 = F3 * D3 * AH3 - D0; // h/2 * l_4 - p_n
  D5 = P3*(D1 + 2*D2 + D3 + D4); //p_n+1

  //
  // Derivatives of track parameters
  //
  if(jacobianT != nullptr){

    // jacobianT         //FIXME seems in magnetic case there are no shortcuts?
    // 1 0 0 0 0 0 0  x
    // 0 1 0 0 0 0 0  y
    // 0 0 1 0 0 0 0  z
    // x x x x x x 0  a_x
    // x x x x x x 0  a_y
    // x x x x x x 0  a_z
    // x x x x x x 1  q/p
    M7x7& J = *jacobianT;

    double   dA0(0), dA2(0), dA3(0), dA4(0), dA5(0), dA6(0);
    double   dB0(0), dB2(0), dB3(0), dB4(0), dB5(0), dB6(0);
    double   dC0(0), dC2(0), dC3(0), dC4(0), dC5(0), dC6(0);
    double   dD0(0), dD1(0), dD2(0), dD3(0), dD4(0); 

    int start(0);

//     if (!calcOnlyLastRowOfJ) {

//       if (!varField) { // FIXME let's be honest and calculate everything everytime
//         // d(x, y, z)/d(x, y, z) submatrix is unit matrix
//         J(0, 0) = 1;  J(1, 1) = 1;  J(2, 2) = 1;
//         // d(ax, ay, az)/d(ax, ay, az) submatrix is 0
//         // start with d(x, y, z)/d(ax, ay, az)
//         start = 3;
//       }

      for(int i=start; i<7; ++i) { 

        //first point
        dD0 = -D0*D0/m_magCharge/sign*J(i,6);
        dA0 = (1/(F0*F0*D0*D0*D0) - 2/D0)*A0*dD0 - (D1-D0)/D0*J(i,3) - F0*A[0]*(J(i,3)*H0[0] + J(i,4)*H0[1] + J(i,5)*H0[2]); // FIXME A true marvel of clarity
        dB0 = (1/(F0*F0*D0*D0*D0) - 2/D0)*B0*dD0 - (D1-D0)/D0*J(i,4) - F0*A[1]*(J(i,3)*H0[0] + J(i,4)*H0[1] + J(i,5)*H0[2]);
        dC0 = (1/(F0*F0*D0*D0*D0) - 2/D0)*C0*dD0 - (D1-D0)/D0*J(i,5) - F0*A[2]*(J(i,3)*H0[0] + J(i,4)*H0[1] + J(i,5)*H0[2]);

        dD1 = dD0 + (1/(F0*F0*D0*D0*D0) - 1/D0)*(D1-D0)*dD0 + F0*D0*(J(i,3)*H0[0] + J(i,4)*H0[1] + J(i,5)*H0[2]);
        dA2 = dA0+J(i, 3);
        dB2 = dB0+J(i, 4);
        dC2 = dC0+J(i, 5);

        //second point
        dD2 = dD0 + (1/(F1*F1*D1*D1*D1) - 1/D1)*(D2-D0)*dD1 + F1*D1*(dA2*H1[0] + dB2*H1[1] + dC2*H1[2]);
        dA3 = J(i,3)+(1/(F1*F1*D1*D1*D1) - 2/D1)*(A2-A[0])*dD1 - (D2-D0)/D1*dA2 - F1*A2*(dA2*H1[0] + dB2*H1[1] + dC2*H1[2]); // FIXME it's only getting better
        dB3 = J(i,4)+(1/(F1*F1*D1*D1*D1) - 2/D1)*(B2-A[1])*dD1 - (D2-D0)/D1*dB2 - F1*B2*(dA2*H1[0] + dB2*H1[1] + dC2*H1[2]);
        dC3 = J(i,5)+(1/(F1*F1*D1*D1*D1) - 2/D1)*(C2-A[2])*dD1 - (D2-D0)/D1*dC2 - F1*C2*(dA2*H1[0] + dB2*H1[1] + dC2*H1[2]);

        dD3 = dD0 + 2*(1/(F2*F2*D2*D2*D2) - 1/D2)*(D3-D0)*dD2 + 2*F2*D2*(dA3*H1[0] + dB3*H1[1] + dC3*H1[2]);
        dA4 = J(i, 3)+(1/(F2*F2*D2*D2*D2) - 2/D2)*(A3-A[0])*dD2 - (D3-D0)/D2*dA3 - F2*A3*(dA3*H1[0] + dB3*H1[1] + dC3*H1[2]);
        dB4 = J(i, 4)+(1/(F2*F2*D2*D2*D2) - 2/D2)*(B3-A[1])*dD2 - (D3-D0)/D2*dB3 - F2*B3*(dA3*H1[0] + dB3*H1[1] + dC3*H1[2]);
        dC4 = J(i, 5)+(1/(F2*F2*D2*D2*D2) - 2/D2)*(C3-A[2])*dD2 - (D3-D0)/D2*dC3 - F2*C3*(dA3*H1[0] + dB3*H1[1] + dC3*H1[2]);

        //last point
        dA5 = dA4+dA4-J(i, 3);      // }
        dB5 = dB4+dB4-J(i, 4);      //  } =  2*(dA4, dB4, dC4) - dA
        dC5 = dC4+dC4-J(i, 5);      // }

        dD4 = -dD0+(1/(F3*F3*D3*D3*D3) - 1/D3)*(D4+D0)*dD3 + F3*D3*(dA5*H2[0] + dB5*H2[1] + dC5*H2[2]);
        dA6 = (1/(F3*F3*D3*D3*D3) - 2/D3)*(A4-A[0])*dD3 - (D4+D0)/D3*dA5 - F3*A5*(dA5*H2[0] + dB5*H2[1] + dC5*H2[2]);
        dB6 = (1/(F3*F3*D3*D3*D3) - 2/D3)*(B4-A[1])*dD3 - (D4+D0)/D3*dB5 - F3*B5*(dA5*H2[0] + dB5*H2[1] + dC5*H2[2]);
        dC6 = (1/(F3*F3*D3*D3*D3) - 2/D3)*(C4-A[2])*dD3 - (D4+D0)/D3*dC5 - F3*C5*(dA5*H2[0] + dB5*H2[1] + dC5*H2[2]);

        // this gives the same results as multiplying the old with the new Jacobian
        J(i, 0) += (dA2+dA3+dA4)*S3;  J(i, 3) = ((dA0+2.*dA3)+(dA5+dA6))*P3; // dR := dR + S3*[(dA2, dB2, dC2) +   (dA3, dB3, dC3) + (dA4, dB4, dC4)]
        J(i, 1) += (dB2+dB3+dB4)*S3;  J(i, 4) = ((dB0+2.*dB3)+(dB5+dB6))*P3; // dA :=     1/3*[(dA0, dB0, dC0) + 2*(dA3, dB3, dC3) + (dA5, dB5, dC5) + (dA6, dB6, dC6)]
        J(i, 2) += (dC2+dC3+dC4)*S3;  J(i, 5) = ((dC0+2.*dC3)+(dC5+dC6))*P3;
        J(i,6) = -m_magCharge*sign/D5/D5*P3*(dD1 + 2*dD2 + dD3 + dD4);
      }

//     } // end if (!calcOnlyLastRowOfJ)

  }
  //
  // Track parameters in last point
  //
  R[0] += (A2+A3+A4)*S3;   A[0] += (SA[0]=((A0+2.*A3)+(A5+A6))*P3-A[0]);  // R  = R0 + S3*[(A2, B2, C2) +   (A3, B3, C3) + (A4, B4, C4)]
  R[1] += (B2+B3+B4)*S3;   A[1] += (SA[1]=((B0+2.*B3)+(B5+B6))*P3-A[1]);  // A  =     1/3*[(A0, B0, C0) + 2*(A3, B3, C3) + (A5, B5, C5) + (A6, B6, C6)]
  R[2] += (C2+C3+C4)*S3;   A[2] += (SA[2]=((C0+2.*C3)+(C5+C6))*P3-A[2]);  // SA = A_new - A_old
  state7[6] = m_magCharge * sign / D5; // g / p_n+1

  // normalize A
  double CBA ( 1./sqrt(A[0]*A[0]+A[1]*A[1]+A[2]*A[2]) ); // 1/|A|
  A[0] *= CBA; A[1] *= CBA; A[2] *= CBA;


  // Test approximation quality on given step
  double EST ( fabs((A1+A6)-(A3+A4)) +
               fabs((B1+B6)-(B3+B4)) +
               fabs((C1+C6)-(C3+C4))  );  // EST = ||(ABC1+ABC6)-(ABC3+ABC4)||_1  =  ||(axzy x H0 + ABC5 x H2) - (ABC2 x H1 + ABC3 x H1)||_1
  if (debugLvl_ > 0) {
    debugOut << "    RKTrackRep::RKPropagate. Step = "<< S << "; quality EST = " << EST  << " \n";
  }

  // Prevent the step length increase from getting too large, this is
  // just the point where it becomes 10.
  if (EST < DLT*1e-5)
    return 10;

  // Step length increase for a fifth order Runge-Kutta, see e.g. 17.2
  // in Numerical Recipes.  FIXME: move to caller.
  return pow(DLT/EST, 1./5.);
}

void MplTrackRep::Streamer(TBuffer &R__b)
{
   // I guess I have to reimplement this since it can not be taken from RKTrackRep?
   typedef ::genfit::MplTrackRep thisClass;
   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      ::genfit::AbsTrackRep::Streamer(R__b);
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      R__b.CheckByteCount(R__s, R__c, thisClass::IsA());
      lastStartState_.setRep(this);
      lastEndState_.setRep(this);
   } else {
      ::genfit::AbsTrackRep::Streamer(R__b);
      R__c = R__b.WriteVersion(thisClass::IsA(), kTRUE);
      R__b.SetByteCount(R__c, kTRUE);
   }
}
