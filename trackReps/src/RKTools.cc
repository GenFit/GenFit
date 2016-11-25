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

#include <RKTools.h>
#include "IO.h"

#include <TMatrixD.h>
#include <TMatrixDSym.h>

#include <iostream>

namespace genfit {


void RKTools::J_pMTxcov5xJ_pM(const M5x7& J_pM, const M5x5& cov5, M7x7& out7){
  // 5D -> 7D

  // J_pM
  // 0 0 0 0 0 0 1
  // 0 0 0 x x x 0
  // 0 0 0 x x x 0
  // x x x 0 0 0 0
  // x x x 0 0 0 0

  // it is assumed that J_pM is only non-zero here:
  // [6] = 1

  // [1*7+3]
  // [1*7+4]
  // [1*7+5]

  // [2*7+3]
  // [2*7+4]
  // [2*7+5]

  // [3*7+0]
  // [3*7+1]
  // [3*7+2]

  // [4*7+0]
  // [4*7+1]
  // [4*7+2]

  double JTC0  = J_pM[21] * cov5[18] + J_pM[28] * cov5[23];
  double JTC1  = J_pM[21] * cov5[23] + J_pM[28] * cov5[24];
  double JTC2  = J_pM[21] * cov5[16] + J_pM[28] * cov5[21];
  double JTC3  = J_pM[21] * cov5[17] + J_pM[28] * cov5[22];
  double JTC4  = J_pM[22] * cov5[18] + J_pM[29] * cov5[23];
  double JTC5  = J_pM[22] * cov5[23] + J_pM[29] * cov5[24];
  double JTC6  = J_pM[22] * cov5[16] + J_pM[29] * cov5[21];
  double JTC7  = J_pM[22] * cov5[17] + J_pM[29] * cov5[22];
  double JTC8  = J_pM[23] * cov5[18] + J_pM[30] * cov5[23];
  double JTC9  = J_pM[23] * cov5[23] + J_pM[30] * cov5[24];
  double JTC10 = J_pM[23] * cov5[16] + J_pM[30] * cov5[21];
  double JTC11 = J_pM[23] * cov5[17] + J_pM[30] * cov5[22];
  double JTC12 = J_pM[10] * cov5[6]  + J_pM[17] * cov5[11];
  double JTC13 = J_pM[10] * cov5[11] + J_pM[17] * cov5[12];
  double JTC14 = J_pM[11] * cov5[6]  + J_pM[18] * cov5[11];
  double JTC15 = J_pM[11] * cov5[11] + J_pM[18] * cov5[12];

  // loops are vectorizable by the compiler!
  for (int i=0; i<3; ++i) out7[i]  = JTC0 * J_pM[21+i] + JTC1 * J_pM[28+i];
  for (int i=0; i<3; ++i) out7[3+i]  = JTC2 * J_pM[10+i] + JTC3 * J_pM[17+i];
  out7[6]  =  J_pM[21] * cov5[15] + J_pM[28] * cov5[20];

  for (int i=0; i<2; ++i) out7[8+i]  = JTC4 * J_pM[22+i] + JTC5 * J_pM[29+i];
  for (int i=0; i<3; ++i) out7[10+i] = JTC6 * J_pM[10+i] + JTC7 * J_pM[17+i];
  out7[13] = J_pM[22] * cov5[15] + J_pM[29] * cov5[20];

  out7[16] = JTC8 * J_pM[23] + JTC9 * J_pM[30];
  for (int i=0; i<3; ++i) out7[17+i] = JTC10 * J_pM[10+i] + JTC11 * J_pM[17+i];
  out7[20] = J_pM[23] * cov5[15] + J_pM[30] * cov5[20];

  for (int i=0; i<3; ++i) out7[24+i] = JTC12 * J_pM[10+i] + JTC13 * J_pM[17+i];
  out7[27] = J_pM[10] * cov5[5] + J_pM[17] * cov5[10];

  for (int i=0; i<2; ++i) out7[32+i] = JTC14 * J_pM[11+i] + JTC15 * J_pM[18+i];
  out7[34] = J_pM[11] * cov5[5] + J_pM[18] * cov5[10];

  out7[40] = (J_pM[12] * cov5[6] + J_pM[19] * cov5[11]) * J_pM[12] + (J_pM[12] * cov5[11] + J_pM[19] * cov5[12]) * J_pM[19];
  out7[41] = J_pM[12] * cov5[5] + J_pM[19] * cov5[10];

  out7[48] = cov5[0];

  // symmetric part
  out7[7]  = out7[1];
  out7[14] = out7[2];  out7[15] = out7[9];
  out7[21] = out7[3];  out7[22] = out7[10];  out7[23] = out7[17];
  out7[28] = out7[4];  out7[29] = out7[11];  out7[30] = out7[18];  out7[31] = out7[25];
  out7[35] = out7[5];  out7[36] = out7[12];  out7[37] = out7[19];  out7[38] = out7[26];  out7[39] = out7[33];
  out7[42] = out7[6];  out7[43] = out7[13];  out7[44] = out7[20];  out7[45] = out7[27];  out7[46] = out7[34];  out7[47] = out7[41];

}


void RKTools::J_pMTxcov5xJ_pM(const M5x6& J_pM, const M5x5& cov5, M6x6& out6){
  // 5D -> 6D

  // J_pM
  // 0 0 0 x x x
  // 0 0 0 x x x
  // 0 0 0 x x x
  // x x x 0 0 0
  // x x x 0 0 0

  // it is assumed that J_pM is only non-zero here:
  // [3]
  // [4]
  // [5]

  // [1*6+3]
  // [1*6+4]
  // [1*6+5]

  // [2*6+3]
  // [2*6+4]
  // [2*6+5]

  // [3*6+0]
  // [3*6+1]
  // [3*6+2]

  // [4*6+0]
  // [4*6+1]
  // [4*6+2]

  double JTC0  = J_pM[18] * cov5[15+3] + J_pM[24] * cov5[20+3];
  double JTC1  = J_pM[18] * cov5[20+3] + J_pM[24] * cov5[20+4];
  double JTC2  = J_pM[18] * cov5[15] + J_pM[24] * cov5[20];
  double JTC3  = J_pM[18] * cov5[15+1] + J_pM[24] * cov5[20+1];
  double JTC4  = J_pM[18] * cov5[15+2] + J_pM[24] * cov5[20+2];
  double JTC5  = J_pM[18+1] * cov5[15+3] + J_pM[24+1] * cov5[20+3];
  double JTC6  = J_pM[18+1] * cov5[20+3] + J_pM[24+1] * cov5[20+4];
  double JTC7  = J_pM[18+1] * cov5[15] + J_pM[24+1] * cov5[20];
  double JTC8  = J_pM[18+1] * cov5[15+1] + J_pM[24+1] * cov5[20+1];
  double JTC9  = J_pM[18+1] * cov5[15+2] + J_pM[24+1] * cov5[20+2];
  double JTC10 = J_pM[18+2] * cov5[15] + J_pM[24+2] * cov5[20];
  double JTC11 = J_pM[18+2] * cov5[15+1] + J_pM[24+2] * cov5[20+1];
  double JTC12 = J_pM[18+2] * cov5[15+2] + J_pM[24+2] * cov5[20+2];
  double JTC13 = J_pM[3] * cov5[0*5] + J_pM[6+3] * cov5[5] + J_pM[12+3] * cov5[10];
  double JTC14 = J_pM[3] * cov5[5] + J_pM[6+3] * cov5[5+1] + J_pM[12+3] * cov5[10+1];
  double JTC15 = J_pM[3] * cov5[10] + J_pM[6+3] * cov5[10+1] + J_pM[12+3] * cov5[10+2];
  double JTC16 = J_pM[4] * cov5[0*5] + J_pM[6+4] * cov5[5] + J_pM[12+4] * cov5[10];
  double JTC17 = J_pM[4] * cov5[5] + J_pM[6+4] * cov5[5+1] + J_pM[12+4] * cov5[10+1];
  double JTC18 = J_pM[4] * cov5[10] + J_pM[6+4] * cov5[10+1] + J_pM[12+4] * cov5[10+2];

  // loops are vectorizable by the compiler!
  for (int i=0; i<3; ++i) out6[i] = JTC0 * J_pM[18+i] + JTC1 * J_pM[24+i];
  for (int i=0; i<3; ++i) out6[3+i] = JTC2 * J_pM[3+i] + JTC3 * J_pM[9+i] + JTC4 * J_pM[15+i];

  for (int i=0; i<2; ++i) out6[7+i] = JTC5 * J_pM[19+i] + JTC6 * J_pM[25+i];
  for (int i=0; i<3; ++i) out6[9+i] = JTC7 * J_pM[3+i] + JTC8 * J_pM[9+i] + JTC9 * J_pM[15+i];

  out6[12+2] = (J_pM[18+2] * cov5[15+3] + J_pM[24+2] * cov5[20+3]) * J_pM[18+2] + (J_pM[18+2] * cov5[20+3] + J_pM[24+2] * cov5[20+4]) * J_pM[24+2];
  for (int i=0; i<3; ++i) out6[15+i] = JTC10 * J_pM[3+i] + JTC11 * J_pM[9+i] + JTC12 * J_pM[15+i];

  for (int i=0; i<3; ++i) out6[21+i] = JTC13 * J_pM[3+i] + JTC14 * J_pM[9+i] + JTC15 * J_pM[15+i];

  for (int i=0; i<3; ++i) out6[28+i] = JTC16 * J_pM[4+i] + JTC17 * J_pM[10+i] + JTC18 * J_pM[16+i];

  out6[30+5] = (J_pM[5] * cov5[0*5] + J_pM[6+5] * cov5[5] + J_pM[12+5] * cov5[10]) * J_pM[5] + (J_pM[5] * cov5[5] + J_pM[6+5] * cov5[5+1] + J_pM[12+5] * cov5[10+1]) * J_pM[6+5] + (J_pM[5] * cov5[10] + J_pM[6+5] * cov5[10+1] + J_pM[12+5] * cov5[10+2]) * J_pM[12+5];

  // symmetric part
  out6[6] = out6[1];
  out6[12] = out6[2];  out6[12+1] = out6[6+2];
  out6[18] = out6[3];  out6[18+1] = out6[6+3];  out6[18+2] = out6[12+3];
  out6[24] = out6[4];  out6[24+1] = out6[6+4];  out6[24+2] = out6[12+4];  out6[24+3] = out6[18+4];
  out6[30] = out6[5];  out6[30+1] = out6[6+5];  out6[30+2] = out6[12+5];  out6[30+3] = out6[18+5];  out6[30+4] = out6[24+5];

}


void RKTools::J_MpTxcov7xJ_Mp(const M7x5& J_Mp, const M7x7& cov7, M5x5& out5)
{
  // 7D -> 5D

  // J_Mp
  // 0 0 0 x x
  // 0 0 0 x x
  // 0 0 0 x x
  // 0 x x 0 0
  // 0 x x 0 0
  // 0 x x 0 0
  // 1 0 0 0 0

  // it is assumed that J_Mp is only non-zero here:
  // [3]
  // [1*7+3]
  // [2*7+3]

  // [4]
  // [1*7+4]
  // [2*7+4]

  // [3*7+1]
  // [4*7+1]
  // [5*7+1]

  // [3*7+2]
  // [4*7+2]
  // [5*7+2]

  // [6*7+0] = 1


  double JTC0  = (J_Mp[16] * cov7[24] + J_Mp[21] * cov7[31] + J_Mp[26] * cov7[38]);
  double JTC1  = (J_Mp[16] * cov7[31] + J_Mp[21] * cov7[32] + J_Mp[26] * cov7[39]);
  double JTC2  = (J_Mp[16] * cov7[38] + J_Mp[21] * cov7[39] + J_Mp[26] * cov7[40]);
  double JTC3  = (J_Mp[16] * cov7[21] + J_Mp[21] * cov7[28] + J_Mp[26] * cov7[35]);
  double JTC4  = (J_Mp[16] * cov7[22] + J_Mp[21] * cov7[29] + J_Mp[26] * cov7[36]);
  double JTC5  = (J_Mp[16] * cov7[23] + J_Mp[21] * cov7[30] + J_Mp[26] * cov7[37]);
  double JTC6  = (J_Mp[17] * cov7[21] + J_Mp[22] * cov7[28] + J_Mp[27] * cov7[35]);
  double JTC7  = (J_Mp[17] * cov7[22] + J_Mp[22] * cov7[29] + J_Mp[27] * cov7[36]);
  double JTC8  = (J_Mp[17] * cov7[23] + J_Mp[22] * cov7[30] + J_Mp[27] * cov7[37]);
  double JTC9  = (J_Mp[3] * cov7[0] + J_Mp[8] * cov7[7] + J_Mp[13] * cov7[14]);
  double JTC10 = (J_Mp[3] * cov7[7] + J_Mp[8] * cov7[8] + J_Mp[13] * cov7[15]);
  double JTC11 = (J_Mp[3] * cov7[14] + J_Mp[8] * cov7[15] + J_Mp[13] * cov7[16]);

  out5[0] = cov7[48];
  out5[1] = J_Mp[16] * cov7[45] + J_Mp[21] * cov7[46] + J_Mp[26] * cov7[47];
  out5[2] = J_Mp[17] * cov7[45] + J_Mp[22] * cov7[46] + J_Mp[27] * cov7[47];
  out5[3] = J_Mp[3] * cov7[42] + J_Mp[8] * cov7[43] + J_Mp[13] * cov7[44];
  out5[4] = J_Mp[4] * cov7[42] + J_Mp[9] * cov7[43] + J_Mp[14] * cov7[44];

  // loops are vectorizable by the compiler!
  for (int i=0; i<2; ++i) out5[6+i] = JTC0 * J_Mp[16+i] + JTC1 * J_Mp[21+i] + JTC2 * J_Mp[26+i];
  for (int i=0; i<2; ++i) out5[8+i] = JTC3 * J_Mp[3+i] + JTC4 * J_Mp[8+i] + JTC5 * J_Mp[13+i];

  out5[12] = (J_Mp[17] * cov7[24] + J_Mp[22] * cov7[31] + J_Mp[27] * cov7[38]) * J_Mp[17] + (J_Mp[17] * cov7[31] + J_Mp[22] * cov7[32] + J_Mp[27] * cov7[39]) * J_Mp[22] + (J_Mp[17] * cov7[38] + J_Mp[22] * cov7[39] + J_Mp[27] * cov7[40]) * J_Mp[27];
  for (int i=0; i<2; ++i) out5[13+i] = JTC6 * J_Mp[3+i] + JTC7 * J_Mp[8+i] + JTC8 * J_Mp[13+i];

  for (int i=0; i<2; ++i) out5[18+i] = JTC9 * J_Mp[3+i] + JTC10 * J_Mp[8+i] + JTC11 * J_Mp[13+i];

  out5[24] = (J_Mp[4] * cov7[0] + J_Mp[9] * cov7[7] + J_Mp[14] * cov7[14]) * J_Mp[4] + (J_Mp[4] * cov7[7] + J_Mp[9] * cov7[8] + J_Mp[14] * cov7[15]) * J_Mp[9] + (J_Mp[4] * cov7[14] + J_Mp[9] * cov7[15] + J_Mp[14] * cov7[16]) * J_Mp[14];

  // symmetric part
  out5[5] = out5[1];
  out5[10] = out5[2];  out5[11] = out5[7];
  out5[15] = out5[3];  out5[16] = out5[8];  out5[17] = out5[13];
  out5[20] = out5[4];  out5[21] = out5[9];  out5[22] = out5[14];  out5[23] = out5[19];

}


void RKTools::J_MpTxcov6xJ_Mp(const M6x5& J_Mp, const M6x6& cov6, M5x5& out5)
{
  // 6D -> 5D

  // J_Mp
  // 0 0 0 x x
  // 0 0 0 x x
  // 0 0 0 x x
  // x x x 0 0
  // x x x 0 0
  // x x x 0 0

  // it is assumed that J_Mp is only non-zero here:
  // [3]
  // [1*6+3]
  // [2*6+3]

  // [4]
  // [1*6+4]
  // [2*6+4]

  // [3*6+0]
  // [4*6+0]
  // [5*6+0]

  // [3*6+1]
  // [4*6+1]
  // [5*6+1]

  // [3*6+2]
  // [4*6+2]
  // [5*6+2]

  double JTC0  = (J_Mp[15] * cov6[18+3] + J_Mp[20] * cov6[24+3] + J_Mp[25] * cov6[30+3]);
  double JTC1  = (J_Mp[15] * cov6[24+3] + J_Mp[20] * cov6[24+4] + J_Mp[25] * cov6[30+4]);
  double JTC2  = (J_Mp[15] * cov6[30+3] + J_Mp[20] * cov6[30+4] + J_Mp[25] * cov6[30+5]);
  double JTC3  = (J_Mp[15] * cov6[18] + J_Mp[20] * cov6[24] + J_Mp[25] * cov6[30]);
  double JTC4  = (J_Mp[15] * cov6[18+1] + J_Mp[20] * cov6[24+1] + J_Mp[25] * cov6[30+1]);
  double JTC5  = (J_Mp[15] * cov6[18+2] + J_Mp[20] * cov6[24+2] + J_Mp[25] * cov6[30+2]);
  double JTC6  = (J_Mp[15+1] * cov6[18+3] + J_Mp[20+1] * cov6[24+3] + J_Mp[25+1] * cov6[30+3]);
  double JTC7  = (J_Mp[15+1] * cov6[24+3] + J_Mp[20+1] * cov6[24+4] + J_Mp[25+1] * cov6[30+4]);
  double JTC8  = (J_Mp[15+1] * cov6[30+3] + J_Mp[20+1] * cov6[30+4] + J_Mp[25+1] * cov6[30+5]);
  double JTC9  = (J_Mp[15+1] * cov6[18] + J_Mp[20+1] * cov6[24] + J_Mp[25+1] * cov6[30]);
  double JTC10 = (J_Mp[15+1] * cov6[18+1] + J_Mp[20+1] * cov6[24+1] + J_Mp[25+1] * cov6[30+1]);
  double JTC11 = (J_Mp[15+1] * cov6[18+2] + J_Mp[20+1] * cov6[24+2] + J_Mp[25+1] * cov6[30+2]);
  double JTC12 = (J_Mp[15+2] * cov6[18] + J_Mp[20+2] * cov6[24] + J_Mp[25+2] * cov6[30]);
  double JTC13 = (J_Mp[15+2] * cov6[18+1] + J_Mp[20+2] * cov6[24+1] + J_Mp[25+2] * cov6[30+1]);
  double JTC14 = (J_Mp[15+2] * cov6[18+2] + J_Mp[20+2] * cov6[24+2] + J_Mp[25+2] * cov6[30+2]);
  double JTC15 = (J_Mp[3] * cov6[0] + J_Mp[5+3] * cov6[6] + J_Mp[10+3] * cov6[12]);
  double JTC16 = (J_Mp[3] * cov6[6] + J_Mp[5+3] * cov6[6+1] + J_Mp[10+3] * cov6[12+1]);
  double JTC17 = (J_Mp[3] * cov6[12] + J_Mp[5+3] * cov6[12+1] + J_Mp[10+3] * cov6[12+2]);

  // loops are vectorizable by the compiler!
  for (int i=0; i<3; ++i) out5[i] = JTC0 * J_Mp[15+i] + JTC1 * J_Mp[20+i] + JTC2 * J_Mp[25+i];
  for (int i=0; i<2; ++i) out5[3+i] = JTC3 * J_Mp[3+i] + JTC4 * J_Mp[8+i] + JTC5 * J_Mp[13+i];

  for (int i=0; i<2; ++i) out5[6+i] = JTC6 * J_Mp[16+i] + JTC7 * J_Mp[21+i] + JTC8 * J_Mp[26+i];
  for (int i=0; i<2; ++i) out5[8+i] = JTC9 * J_Mp[3+i] + JTC10 * J_Mp[8+i] + JTC11 * J_Mp[13+i];

  out5[10+2] = (J_Mp[15+2] * cov6[18+3] + J_Mp[20+2] * cov6[24+3] + J_Mp[25+2] * cov6[30+3]) * J_Mp[15+2] + (J_Mp[15+2] * cov6[24+3] + J_Mp[20+2] * cov6[24+4] + J_Mp[25+2] * cov6[30+4]) * J_Mp[20+2] + (J_Mp[15+2] * cov6[30+3] + J_Mp[20+2] * cov6[30+4] + J_Mp[25+2] * cov6[30+5]) * J_Mp[25+2];
  for (int i=0; i<2; ++i) out5[13+i] = JTC12 * J_Mp[3+i] + JTC13 * J_Mp[8+i] + JTC14 * J_Mp[13+i];

  for (int i=0; i<2; ++i) out5[18+i] = JTC15 * J_Mp[3+i] + JTC16 * J_Mp[8+i] + JTC17 * J_Mp[13+i];

  out5[20+4] = (J_Mp[4] * cov6[0] + J_Mp[5+4] * cov6[6] + J_Mp[10+4] * cov6[12]) * J_Mp[4] + (J_Mp[4] * cov6[6] + J_Mp[5+4] * cov6[6+1] + J_Mp[10+4] * cov6[12+1]) * J_Mp[5+4] + (J_Mp[4] * cov6[12] + J_Mp[5+4] * cov6[12+1] + J_Mp[10+4] * cov6[12+2]) * J_Mp[10+4];

  // symmetric part
  out5[5] = out5[1];
  out5[10] = out5[2];  out5[10+1] = out5[5+2];
  out5[15] = out5[3];  out5[15+1] = out5[5+3];  out5[15+2] = out5[10+3];
  out5[20] = out5[4];  out5[20+1] = out5[5+4];  out5[20+2] = out5[10+4];  out5[20+3] = out5[15+4];

}

void RKTools::J_pMTTxJ_MMTTxJ_MpTT(const M7x5& J_pMT, const M7x7& J_MMT, const M5x7& J_MpT, M5x5& J_pp)
{
    // calculates  J_pp = J_pM * J_MM * J_Mp
    // input J_MMT is transposed version of actual jacobian J_MM
    // input J_pMT is transposed version of actual jacobian J_pM (Master to plane)
    // input J_MpT is transposed version of actual jacobian J_Mp (plane to Master)

  // J_pMT
  // 0 0 0 x x
  // 0 0 0 x x
  // 0 0 0 x x
  // 0 x x 0 0
  // 0 x x 0 0
  // 0 x x 0 0
  // 1 0 0 0 0

  // J_MMT if MMproj == true
  // x x x x x x 0
  // x x x x x x 0
  // x x x x x x 0
  // x x x x x x 0
  // x x x x x x 0
  // x x x x x x 0
  // x x x x x x x

  // J_MpT
  // 0 0 0 0 0 0 1
  // 0 0 0 x x x 0
  // 0 0 0 x x x 0
  // x x x 0 0 0 0
  // x x x 0 0 0 0


  J_pp[0*5+0] = J_MMT[6*7+6];
  J_pp[0*5+1] = 0;
  J_pp[0*5+2] = 0;
  J_pp[0*5+3] = 0;
  J_pp[0*5+4] = 0;

  J_pp[1*5+0] = J_pMT[3*5+1] * J_MMT[6*7+3] + J_pMT[4*5+1] * J_MMT[6*7+4] + J_pMT[5*5+1] * J_MMT[6*7+5];
  J_pp[1*5+1] = (  (J_pMT[3*5+1] * J_MMT[3*7+3] + J_pMT[4*5+1] * J_MMT[3*7+4] + J_pMT[5*5+1] * J_MMT[3*7+5]) * J_MpT[1*7+3] 
                 + (J_pMT[3*5+1] * J_MMT[4*7+3] + J_pMT[4*5+1] * J_MMT[4*7+4] + J_pMT[5*5+1] * J_MMT[4*7+5]) * J_MpT[1*7+4]
                 + (J_pMT[3*5+1] * J_MMT[5*7+3] + J_pMT[4*5+1] * J_MMT[5*7+4] + J_pMT[5*5+1] * J_MMT[5*7+5]) * J_MpT[1*7+5]);
  J_pp[1*5+2] = (  (J_pMT[3*5+1] * J_MMT[3*7+3] + J_pMT[4*5+1] * J_MMT[3*7+4] + J_pMT[5*5+1] * J_MMT[3*7+5]) * J_MpT[2*7+3]
                 + (J_pMT[3*5+1] * J_MMT[4*7+3] + J_pMT[4*5+1] * J_MMT[4*7+4] + J_pMT[5*5+1] * J_MMT[4*7+5]) * J_MpT[2*7+4]
                 + (J_pMT[3*5+1] * J_MMT[5*7+3] + J_pMT[4*5+1] * J_MMT[5*7+4] + J_pMT[5*5+1] * J_MMT[5*7+5]) * J_MpT[2*7+5]);
  J_pp[1*5+3] = (  (J_pMT[3*5+1] * J_MMT[0*7+3] + J_pMT[4*5+1] * J_MMT[0*7+4] + J_pMT[5*5+1] * J_MMT[0*7+5]) * J_MpT[3*7+0]
                 + (J_pMT[3*5+1] * J_MMT[1*7+3] + J_pMT[4*5+1] * J_MMT[1*7+4] + J_pMT[5*5+1] * J_MMT[1*7+5]) * J_MpT[3*7+1]
                 + (J_pMT[3*5+1] * J_MMT[2*7+3] + J_pMT[4*5+1] * J_MMT[2*7+4] + J_pMT[5*5+1] * J_MMT[2*7+5]) * J_MpT[3*7+2]);
  J_pp[1*5+4] = (  (J_pMT[3*5+1] * J_MMT[0*7+3] + J_pMT[4*5+1] * J_MMT[0*7+4] + J_pMT[5*5+1] * J_MMT[0*7+5]) * J_MpT[4*7+0]
                 + (J_pMT[3*5+1] * J_MMT[1*7+3] + J_pMT[4*5+1] * J_MMT[1*7+4] + J_pMT[5*5+1] * J_MMT[1*7+5]) * J_MpT[4*7+1]
                 + (J_pMT[3*5+1] * J_MMT[2*7+3] + J_pMT[4*5+1] * J_MMT[2*7+4] + J_pMT[5*5+1] * J_MMT[2*7+5]) * J_MpT[4*7+2]);

  J_pp[2*5+0] =  J_pMT[3*5+2] * J_MMT[6*7+3] + J_pMT[4*5+2] * J_MMT[6*7+4] + J_pMT[5*5+2] * J_MMT[6*7+5];
  J_pp[2*5+1] = (  (J_pMT[3*5+2] * J_MMT[3*7+3] + J_pMT[4*5+2] * J_MMT[3*7+4] + J_pMT[5*5+2] * J_MMT[3*7+5]) * J_MpT[1*7+3]
                 + (J_pMT[3*5+2] * J_MMT[4*7+3] + J_pMT[4*5+2] * J_MMT[4*7+4] + J_pMT[5*5+2] * J_MMT[4*7+5]) * J_MpT[1*7+4]
                 + (J_pMT[3*5+2] * J_MMT[5*7+3] + J_pMT[4*5+2] * J_MMT[5*7+4] + J_pMT[5*5+2] * J_MMT[5*7+5]) * J_MpT[1*7+5]);
  J_pp[2*5+2] = (  (J_pMT[3*5+2] * J_MMT[3*7+3] + J_pMT[4*5+2] * J_MMT[3*7+4] + J_pMT[5*5+2] * J_MMT[3*7+5]) * J_MpT[2*7+3]
                 + (J_pMT[3*5+2] * J_MMT[4*7+3] + J_pMT[4*5+2] * J_MMT[4*7+4] + J_pMT[5*5+2] * J_MMT[4*7+5]) * J_MpT[2*7+4]
                 + (J_pMT[3*5+2] * J_MMT[5*7+3] + J_pMT[4*5+2] * J_MMT[5*7+4] + J_pMT[5*5+2] * J_MMT[5*7+5]) * J_MpT[2*7+5]);
  J_pp[2*5+3] = (  (J_pMT[3*5+2] * J_MMT[0*7+3] + J_pMT[4*5+2] * J_MMT[0*7+4] + J_pMT[5*5+2] * J_MMT[0*7+5]) * J_MpT[3*7+0]
                 + (J_pMT[3*5+2] * J_MMT[1*7+3] + J_pMT[4*5+2] * J_MMT[1*7+4] + J_pMT[5*5+2] * J_MMT[1*7+5]) * J_MpT[3*7+1]
                 + (J_pMT[3*5+2] * J_MMT[2*7+3] + J_pMT[4*5+2] * J_MMT[2*7+4] + J_pMT[5*5+2] * J_MMT[2*7+5]) * J_MpT[3*7+2]);
  J_pp[2*5+4] = (  (J_pMT[3*5+2] * J_MMT[0*7+3] + J_pMT[4*5+2] * J_MMT[0*7+4] + J_pMT[5*5+2] * J_MMT[0*7+5]) * J_MpT[4*7+0]
                 + (J_pMT[3*5+2] * J_MMT[1*7+3] + J_pMT[4*5+2] * J_MMT[1*7+4] + J_pMT[5*5+2] * J_MMT[1*7+5]) * J_MpT[4*7+1]
                 + (J_pMT[3*5+2] * J_MMT[2*7+3] + J_pMT[4*5+2] * J_MMT[2*7+4] + J_pMT[5*5+2] * J_MMT[2*7+5]) * J_MpT[4*7+2]);

  J_pp[3*5+0] =  J_pMT[0*5+3] * J_MMT[6*7+0] + J_pMT[1*5+3] * J_MMT[6*7+1] + J_pMT[2*5+3] * J_MMT[6*7+2];
  J_pp[3*5+1] = (  (J_pMT[0*5+3] * J_MMT[3*7+0] + J_pMT[1*5+3] * J_MMT[3*7+1] + J_pMT[2*5+3] * J_MMT[3*7+2]) * J_MpT[1*7+3]
                 + (J_pMT[0*5+3] * J_MMT[4*7+0] + J_pMT[1*5+3] * J_MMT[4*7+1] + J_pMT[2*5+3] * J_MMT[4*7+2]) * J_MpT[1*7+4]
                 + (J_pMT[0*5+3] * J_MMT[5*7+0] + J_pMT[1*5+3] * J_MMT[5*7+1] + J_pMT[2*5+3] * J_MMT[5*7+2]) * J_MpT[1*7+5]);
  J_pp[3*5+2] = (  (J_pMT[0*5+3] * J_MMT[3*7+0] + J_pMT[1*5+3] * J_MMT[3*7+1] + J_pMT[2*5+3] * J_MMT[3*7+2]) * J_MpT[2*7+3]
                 + (J_pMT[0*5+3] * J_MMT[4*7+0] + J_pMT[1*5+3] * J_MMT[4*7+1] + J_pMT[2*5+3] * J_MMT[4*7+2]) * J_MpT[2*7+4]
                 + (J_pMT[0*5+3] * J_MMT[5*7+0] + J_pMT[1*5+3] * J_MMT[5*7+1] + J_pMT[2*5+3] * J_MMT[5*7+2]) * J_MpT[2*7+5]);
  J_pp[3*5+3] = (  (J_pMT[0*5+3] * J_MMT[0*7+0] + J_pMT[1*5+3] * J_MMT[0*7+1] + J_pMT[2*5+3] * J_MMT[0*7+2]) * J_MpT[3*7+0]
                 + (J_pMT[0*5+3] * J_MMT[1*7+0] + J_pMT[1*5+3] * J_MMT[1*7+1] + J_pMT[2*5+3] * J_MMT[1*7+2]) * J_MpT[3*7+1]
                 + (J_pMT[0*5+3] * J_MMT[2*7+0] + J_pMT[1*5+3] * J_MMT[2*7+1] + J_pMT[2*5+3] * J_MMT[2*7+2]) * J_MpT[3*7+2]);
  J_pp[3*5+4] = (  (J_pMT[0*5+3] * J_MMT[0*7+0] + J_pMT[1*5+3] * J_MMT[0*7+1] + J_pMT[2*5+3] * J_MMT[0*7+2]) * J_MpT[4*7+0]
                 + (J_pMT[0*5+3] * J_MMT[1*7+0] + J_pMT[1*5+3] * J_MMT[1*7+1] + J_pMT[2*5+3] * J_MMT[1*7+2]) * J_MpT[4*7+1]
                 + (J_pMT[0*5+3] * J_MMT[2*7+0] + J_pMT[1*5+3] * J_MMT[2*7+1] + J_pMT[2*5+3] * J_MMT[2*7+2]) * J_MpT[4*7+2]);

  J_pp[4*5+0] =  J_pMT[0*5+4] * J_MMT[6*7+0] + J_pMT[1*5+4] * J_MMT[6*7+1] + J_pMT[2*5+4] * J_MMT[6*7+2];
  J_pp[4*5+1] = (  (J_pMT[0*5+4] * J_MMT[3*7+0] + J_pMT[1*5+4] * J_MMT[3*7+1] + J_pMT[2*5+4] * J_MMT[3*7+2]) * J_MpT[1*7+3]
                 + (J_pMT[0*5+4] * J_MMT[4*7+0] + J_pMT[1*5+4] * J_MMT[4*7+1] + J_pMT[2*5+4] * J_MMT[4*7+2]) * J_MpT[1*7+4]
                 + (J_pMT[0*5+4] * J_MMT[5*7+0] + J_pMT[1*5+4] * J_MMT[5*7+1] + J_pMT[2*5+4] * J_MMT[5*7+2]) * J_MpT[1*7+5]);
  J_pp[4*5+2] = (  (J_pMT[0*5+4] * J_MMT[3*7+0] + J_pMT[1*5+4] * J_MMT[3*7+1] + J_pMT[2*5+4] * J_MMT[3*7+2]) * J_MpT[2*7+3]
                 + (J_pMT[0*5+4] * J_MMT[4*7+0] + J_pMT[1*5+4] * J_MMT[4*7+1] + J_pMT[2*5+4] * J_MMT[4*7+2]) * J_MpT[2*7+4]
                 + (J_pMT[0*5+4] * J_MMT[5*7+0] + J_pMT[1*5+4] * J_MMT[5*7+1] + J_pMT[2*5+4] * J_MMT[5*7+2]) * J_MpT[2*7+5]);
  J_pp[4*5+3] = (  (J_pMT[0*5+4] * J_MMT[0*7+0] + J_pMT[1*5+4] * J_MMT[0*7+1] + J_pMT[2*5+4] * J_MMT[0*7+2]) * J_MpT[3*7+0]
                 + (J_pMT[0*5+4] * J_MMT[1*7+0] + J_pMT[1*5+4] * J_MMT[1*7+1] + J_pMT[2*5+4] * J_MMT[1*7+2]) * J_MpT[3*7+1]
                 + (J_pMT[0*5+4] * J_MMT[2*7+0] + J_pMT[1*5+4] * J_MMT[2*7+1] + J_pMT[2*5+4] * J_MMT[2*7+2]) * J_MpT[3*7+2]);
  J_pp[4*5+4] = (  (J_pMT[0*5+4] * J_MMT[0*7+0] + J_pMT[1*5+4] * J_MMT[0*7+1] + J_pMT[2*5+4] * J_MMT[0*7+2]) * J_MpT[4*7+0]
                 + (J_pMT[0*5+4] * J_MMT[1*7+0] + J_pMT[1*5+4] * J_MMT[1*7+1] + J_pMT[2*5+4] * J_MMT[1*7+2]) * J_MpT[4*7+1]
                 + (J_pMT[0*5+4] * J_MMT[2*7+0] + J_pMT[1*5+4] * J_MMT[2*7+1] + J_pMT[2*5+4] * J_MMT[2*7+2]) * J_MpT[4*7+2]);
}


void RKTools::Np_N_NpT(const M7x7& Np, M7x7& N) {

  // N is symmetric

  // Np
  // x x x 0 0 0 0
  // x x x 0 0 0 0
  // x x x 0 0 0 0
  // x x x 1 0 0 0
  // x x x 0 1 0 0
  // x x x 0 0 1 0
  // 0 0 0 0 0 0 1

  // calculate:
  // Np * N * Np^T

  double N00(N[0*7+0]), N11(N[1*7+1]), N22(N[2*7+2]), N33(N[3*7+3]), N44(N[4*7+4]), N55(N[5*7+5]);

  // replace lower left triangle using the original values from upper right triangle plus the cached diagonal values
  N[0*7+0] = (Np[0*7+0] * N00 + Np[0*7+1] * N[0*7+1] + Np[0*7+2] * N[0*7+2]) * Np[0*7+0] + (Np[0*7+0] * N[0*7+1] + Np[0*7+1] * N11 + Np[0*7+2] * N[1*7+2]) * Np[0*7+1] + (Np[0*7+0] * N[0*7+2] + Np[0*7+1] * N[1*7+2] + Np[0*7+2] * N22) * Np[0*7+2];

  N[1*7+0] = (Np[1*7+0] * N00 + Np[1*7+1] * N[0*7+1] + Np[1*7+2] * N[0*7+2]) * Np[0*7+0] + (Np[1*7+0] * N[0*7+1] + Np[1*7+1] * N11 + Np[1*7+2] * N[1*7+2]) * Np[0*7+1] + (Np[1*7+0] * N[0*7+2] + Np[1*7+1] * N[1*7+2] + Np[1*7+2] * N22) * Np[0*7+2];
  N[1*7+1] = (Np[1*7+0] * N00 + Np[1*7+1] * N[0*7+1] + Np[1*7+2] * N[0*7+2]) * Np[1*7+0] + (Np[1*7+0] * N[0*7+1] + Np[1*7+1] * N11 + Np[1*7+2] * N[1*7+2]) * Np[1*7+1] + (Np[1*7+0] * N[0*7+2] + Np[1*7+1] * N[1*7+2] + Np[1*7+2] * N22) * Np[1*7+2];

  N[2*7+0] = (Np[2*7+0] * N00 + Np[2*7+1] * N[0*7+1] + Np[2*7+2] * N[0*7+2]) * Np[0*7+0] + (Np[2*7+0] * N[0*7+1] + Np[2*7+1] * N11 + Np[2*7+2] * N[1*7+2]) * Np[0*7+1] + (Np[2*7+0] * N[0*7+2] + Np[2*7+1] * N[1*7+2] + Np[2*7+2] * N22) * Np[0*7+2];
  N[2*7+1] = (Np[2*7+0] * N00 + Np[2*7+1] * N[0*7+1] + Np[2*7+2] * N[0*7+2]) * Np[1*7+0] + (Np[2*7+0] * N[0*7+1] + Np[2*7+1] * N11 + Np[2*7+2] * N[1*7+2]) * Np[1*7+1] + (Np[2*7+0] * N[0*7+2] + Np[2*7+1] * N[1*7+2] + Np[2*7+2] * N22) * Np[1*7+2];
  N[2*7+2] = (Np[2*7+0] * N00 + Np[2*7+1] * N[0*7+1] + Np[2*7+2] * N[0*7+2]) * Np[2*7+0] + (Np[2*7+0] * N[0*7+1] + Np[2*7+1] * N11 + Np[2*7+2] * N[1*7+2]) * Np[2*7+1] + (Np[2*7+0] * N[0*7+2] + Np[2*7+1] * N[1*7+2] + Np[2*7+2] * N22) * Np[2*7+2];

  N[3*7+0] = (Np[3*7+0] * N00 + Np[3*7+1] * N[0*7+1] + Np[3*7+2] * N[0*7+2] + N[0*7+3]) * Np[0*7+0] + (Np[3*7+0] * N[0*7+1] + Np[3*7+1] * N11 + Np[3*7+2] * N[1*7+2] + N[1*7+3]) * Np[0*7+1] + (Np[3*7+0] * N[0*7+2] + Np[3*7+1] * N[1*7+2] + Np[3*7+2] * N22 + N[2*7+3]) * Np[0*7+2];
  N[3*7+1] = (Np[3*7+0] * N00 + Np[3*7+1] * N[0*7+1] + Np[3*7+2] * N[0*7+2] + N[0*7+3]) * Np[1*7+0] + (Np[3*7+0] * N[0*7+1] + Np[3*7+1] * N11 + Np[3*7+2] * N[1*7+2] + N[1*7+3]) * Np[1*7+1] + (Np[3*7+0] * N[0*7+2] + Np[3*7+1] * N[1*7+2] + Np[3*7+2] * N22 + N[2*7+3]) * Np[1*7+2];
  N[3*7+2] = (Np[3*7+0] * N00 + Np[3*7+1] * N[0*7+1] + Np[3*7+2] * N[0*7+2] + N[0*7+3]) * Np[2*7+0] + (Np[3*7+0] * N[0*7+1] + Np[3*7+1] * N11 + Np[3*7+2] * N[1*7+2] + N[1*7+3]) * Np[2*7+1] + (Np[3*7+0] * N[0*7+2] + Np[3*7+1] * N[1*7+2] + Np[3*7+2] * N22 + N[2*7+3]) * Np[2*7+2];
  N[3*7+3] = (Np[3*7+0] * N00 + Np[3*7+1] * N[0*7+1] + Np[3*7+2] * N[0*7+2] + N[0*7+3]) * Np[3*7+0] + (Np[3*7+0] * N[0*7+1] + Np[3*7+1] * N11 + Np[3*7+2] * N[1*7+2] + N[1*7+3]) * Np[3*7+1] + (Np[3*7+0] * N[0*7+2] + Np[3*7+1] * N[1*7+2] + Np[3*7+2] * N22 + N[2*7+3]) * Np[3*7+2] + Np[3*7+0] * N[0*7+3] + Np[3*7+1] * N[1*7+3] + Np[3*7+2] * N[2*7+3] + N33;

  N[4*7+0] = (Np[4*7+0] * N00 + Np[4*7+1] * N[0*7+1] + Np[4*7+2] * N[0*7+2] + N[0*7+4]) * Np[0*7+0] + (Np[4*7+0] * N[0*7+1] + Np[4*7+1] * N11 + Np[4*7+2] * N[1*7+2] + N[1*7+4]) * Np[0*7+1] + (Np[4*7+0] * N[0*7+2] + Np[4*7+1] * N[1*7+2] + Np[4*7+2] * N22 + N[2*7+4]) * Np[0*7+2];
  N[4*7+1] = (Np[4*7+0] * N00 + Np[4*7+1] * N[0*7+1] + Np[4*7+2] * N[0*7+2] + N[0*7+4]) * Np[1*7+0] + (Np[4*7+0] * N[0*7+1] + Np[4*7+1] * N11 + Np[4*7+2] * N[1*7+2] + N[1*7+4]) * Np[1*7+1] + (Np[4*7+0] * N[0*7+2] + Np[4*7+1] * N[1*7+2] + Np[4*7+2] * N22 + N[2*7+4]) * Np[1*7+2];
  N[4*7+2] = (Np[4*7+0] * N00 + Np[4*7+1] * N[0*7+1] + Np[4*7+2] * N[0*7+2] + N[0*7+4]) * Np[2*7+0] + (Np[4*7+0] * N[0*7+1] + Np[4*7+1] * N11 + Np[4*7+2] * N[1*7+2] + N[1*7+4]) * Np[2*7+1] + (Np[4*7+0] * N[0*7+2] + Np[4*7+1] * N[1*7+2] + Np[4*7+2] * N22 + N[2*7+4]) * Np[2*7+2];
  N[4*7+3] = (Np[4*7+0] * N00 + Np[4*7+1] * N[0*7+1] + Np[4*7+2] * N[0*7+2] + N[0*7+4]) * Np[3*7+0] + (Np[4*7+0] * N[0*7+1] + Np[4*7+1] * N11 + Np[4*7+2] * N[1*7+2] + N[1*7+4]) * Np[3*7+1] + (Np[4*7+0] * N[0*7+2] + Np[4*7+1] * N[1*7+2] + Np[4*7+2] * N22 + N[2*7+4]) * Np[3*7+2] + Np[4*7+0] * N[0*7+3] + Np[4*7+1] * N[1*7+3] + Np[4*7+2] * N[2*7+3] + N[3*7+4];
  N[4*7+4] = (Np[4*7+0] * N00 + Np[4*7+1] * N[0*7+1] + Np[4*7+2] * N[0*7+2] + N[0*7+4]) * Np[4*7+0] + (Np[4*7+0] * N[0*7+1] + Np[4*7+1] * N11 + Np[4*7+2] * N[1*7+2] + N[1*7+4]) * Np[4*7+1] + (Np[4*7+0] * N[0*7+2] + Np[4*7+1] * N[1*7+2] + Np[4*7+2] * N22 + N[2*7+4]) * Np[4*7+2] + Np[4*7+0] * N[0*7+4] + Np[4*7+1] * N[1*7+4] + Np[4*7+2] * N[2*7+4] + N44;

  N[5*7+0] = (Np[5*7+0] * N00 + Np[5*7+1] * N[0*7+1] + Np[5*7+2] * N[0*7+2] + N[0*7+5]) * Np[0*7+0] + (Np[5*7+0] * N[0*7+1] + Np[5*7+1] * N11 + Np[5*7+2] * N[1*7+2] + N[1*7+5]) * Np[0*7+1] + (Np[5*7+0] * N[0*7+2] + Np[5*7+1] * N[1*7+2] + Np[5*7+2] * N22 + N[2*7+5]) * Np[0*7+2];
  N[5*7+1] = (Np[5*7+0] * N00 + Np[5*7+1] * N[0*7+1] + Np[5*7+2] * N[0*7+2] + N[0*7+5]) * Np[1*7+0] + (Np[5*7+0] * N[0*7+1] + Np[5*7+1] * N11 + Np[5*7+2] * N[1*7+2] + N[1*7+5]) * Np[1*7+1] + (Np[5*7+0] * N[0*7+2] + Np[5*7+1] * N[1*7+2] + Np[5*7+2] * N22 + N[2*7+5]) * Np[1*7+2];
  N[5*7+2] = (Np[5*7+0] * N00 + Np[5*7+1] * N[0*7+1] + Np[5*7+2] * N[0*7+2] + N[0*7+5]) * Np[2*7+0] + (Np[5*7+0] * N[0*7+1] + Np[5*7+1] * N11 + Np[5*7+2] * N[1*7+2] + N[1*7+5]) * Np[2*7+1] + (Np[5*7+0] * N[0*7+2] + Np[5*7+1] * N[1*7+2] + Np[5*7+2] * N22 + N[2*7+5]) * Np[2*7+2];
  N[5*7+3] = (Np[5*7+0] * N00 + Np[5*7+1] * N[0*7+1] + Np[5*7+2] * N[0*7+2] + N[0*7+5]) * Np[3*7+0] + (Np[5*7+0] * N[0*7+1] + Np[5*7+1] * N11 + Np[5*7+2] * N[1*7+2] + N[1*7+5]) * Np[3*7+1] + (Np[5*7+0] * N[0*7+2] + Np[5*7+1] * N[1*7+2] + Np[5*7+2] * N22 + N[2*7+5]) * Np[3*7+2] + Np[5*7+0] * N[0*7+3] + Np[5*7+1] * N[1*7+3] + Np[5*7+2] * N[2*7+3] + N[3*7+5];
  N[5*7+4] = (Np[5*7+0] * N00 + Np[5*7+1] * N[0*7+1] + Np[5*7+2] * N[0*7+2] + N[0*7+5]) * Np[4*7+0] + (Np[5*7+0] * N[0*7+1] + Np[5*7+1] * N11 + Np[5*7+2] * N[1*7+2] + N[1*7+5]) * Np[4*7+1] + (Np[5*7+0] * N[0*7+2] + Np[5*7+1] * N[1*7+2] + Np[5*7+2] * N22 + N[2*7+5]) * Np[4*7+2] + Np[5*7+0] * N[0*7+4] + Np[5*7+1] * N[1*7+4] + Np[5*7+2] * N[2*7+4] + N[4*7+5];
  N[5*7+5] = (Np[5*7+0] * N00 + Np[5*7+1] * N[0*7+1] + Np[5*7+2] * N[0*7+2] + N[0*7+5]) * Np[5*7+0] + (Np[5*7+0] * N[0*7+1] + Np[5*7+1] * N11 + Np[5*7+2] * N[1*7+2] + N[1*7+5]) * Np[5*7+1] + (Np[5*7+0] * N[0*7+2] + Np[5*7+1] * N[1*7+2] + Np[5*7+2] * N22 + N[2*7+5]) * Np[5*7+2] + Np[5*7+0] * N[0*7+5] + Np[5*7+1] * N[1*7+5] + Np[5*7+2] * N[2*7+5] + N55;

  N[6*7+0] = Np[0*7+0] * N[0*7+6] + Np[0*7+1] * N[1*7+6] + Np[0*7+2] * N[2*7+6];
  N[6*7+1] = Np[1*7+0] * N[0*7+6] + Np[1*7+1] * N[1*7+6] + Np[1*7+2] * N[2*7+6];
  N[6*7+2] = Np[2*7+0] * N[0*7+6] + Np[2*7+1] * N[1*7+6] + Np[2*7+2] * N[2*7+6];
  N[6*7+3] = Np[3*7+0] * N[0*7+6] + Np[3*7+1] * N[1*7+6] + Np[3*7+2] * N[2*7+6] + N[3*7+6];
  N[6*7+4] = Np[4*7+0] * N[0*7+6] + Np[4*7+1] * N[1*7+6] + Np[4*7+2] * N[2*7+6] + N[4*7+6];
  N[6*7+5] = Np[5*7+0] * N[0*7+6] + Np[5*7+1] * N[1*7+6] + Np[5*7+2] * N[2*7+6] + N[5*7+6];
  //N[6*7+6] = N[6*7+6]; remains unchanged


  // copy the results from the lower left triangle to the upper right triangle as well
  N[5*7+6] = N[6*7+5];
  N[4*7+5] = N[5*7+4];  N[4*7+6] = N[6*7+4];
  N[3*7+4] = N[4*7+3];  N[3*7+5] = N[5*7+3];  N[3*7+6] = N[6*7+3];
  N[2*7+3] = N[3*7+2];  N[2*7+4] = N[4*7+2];  N[2*7+5] = N[5*7+2];  N[2*7+6] = N[6*7+2];
  N[1*7+2] = N[2*7+1];  N[1*7+3] = N[3*7+1];  N[1*7+4] = N[4*7+1];  N[1*7+5] = N[5*7+1];  N[1*7+6] = N[6*7+1];
  N[0*7+1] = N[1*7+0];  N[0*7+2] = N[2*7+0];  N[0*7+3] = N[3*7+0];  N[0*7+4] = N[4*7+0];  N[0*7+5] = N[5*7+0];  N[0*7+6] = N[6*7+0];
}


void RKTools::printDim(const double* mat, unsigned int dimX, unsigned int dimY){

  printOut << dimX << " x " << dimY << " matrix as follows: \n";
  for (unsigned int i=0; i< dimX; ++i){
    for (unsigned int j=0; j< dimY; ++j){
      printf("  %11.5g", mat[i*dimY+j]);
    }
    printOut<<"\n";
  }
  printOut<<std::endl;

}


} /* End of namespace genfit */

