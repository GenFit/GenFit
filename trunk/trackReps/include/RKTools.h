/* Copyright 2008-2009, Technische Universitaet Muenchen,
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

/** @addtogroup RKTrackRep
 * @{
 */

#ifndef genfit_RKTools_h
#define genfit_RKTools_h

namespace genfit {

/**
 * Array Matrix typedefs. They are needed for SSE optimization:
 * gcc can vectorize loops only if the array sizes are known.
 */
typedef double M1x3[1*3];
typedef double M1x4[1*4];
typedef double M1x6[1*6];
typedef double M1x7[1*7];
typedef double M5x5[5*5];
typedef double M6x6[6*6];
typedef double M7x7[7*7];
typedef double M8x7[8*7];
typedef double M6x5[6*5];
typedef double M7x5[7*5];
typedef double M5x6[5*6];
typedef double M5x7[5*7];

/**
 * @brief Array matrix multiplications used in RKTrackRep
 */
namespace RKTools {

  void J_pMTxcov5xJ_pM(const M5x7& J_pM, const M5x5& cov5, M7x7& out7);
  void J_pMTxcov5xJ_pM(const M5x6& J_pM, const M5x5& cov5, M6x6& out6);

  void J_MpTxcov7xJ_Mp(const M7x5& J_Mp, const M7x7& cov7, M5x5& out5);
  void J_MpTxcov6xJ_Mp(const M6x5& J_Mp, const M6x6& cov6, M5x5& out5);

  void J_MMTxcov7xJ_MM(const M7x7& J_MM, M7x7& cov7);

  void J_MMxJ_MM(M7x7& J_MM, const M7x7& J_MM_old);

  void J_pMTTxJ_MMTTxJ_MpTT(const M7x5& J_pMT, const M7x7& J_MMT, const M5x7& J_MpT, M5x5& J_pp);

  void Np_N_NpT(const M7x7& Np, M7x7& N);

  void printDim(const double* mat, unsigned int dimX, unsigned int dimY);

}

} /* End of namespace genfit */
/** @} */

#endif // genfit_RKTools_h

