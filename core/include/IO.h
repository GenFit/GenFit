/* Copyright 2015, Ludwig-Maximilians-Universität München
   Authors: Tobias Schlüter

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

#ifndef genfit_IO_h
#define genfit_IO_h

/** @brief Defines for I/O streams used for error and debug printing.
 */

#include <ostream>

namespace genfit {

/** Default stream for debug output.  Defaults to std::cout.
   Override destination with debugOut.rdbuf(newStream.rdbuf()).  */
extern std::ostream debugOut;
/** Default stream for error output.  Defaults to std::cerr.
    Override destination with errorOut.rdbuf(newStream.rdbuf()).  */
extern std::ostream errorOut;
/** Default stream for output of Print calls.  Defaults to std::cout.
   Override destination with printOut.rdbuf(newStream.rdbuf()).  */
extern std::ostream printOut;

}

#endif
