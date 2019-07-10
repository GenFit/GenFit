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

#pragma once

namespace genfit {

  /* The PDG code used for the magnetic monopole; somewhat arbitrary.
   *
   * PDG allocates all code starting with 99 to be
   * used for user-defined particles. Some folks weave
   * the magnetic charge and mass into the code itself,
   * but here we would like to consider them as floats,
   * and so we have to store them elsewhere as data members.
   */
  const int c_monopolePDGCode = 99666;

}
