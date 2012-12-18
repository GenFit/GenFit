/* Copyright 2008-2010, Technische Universitaet Muenchen,
   Authors: Christian Hoeppner & Sebastian Neubert

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
/** @addtogroup genfit
 * @{ */


#ifndef GFEXCEPTION_H
#define GFEXCEPTION_H

#include <exception>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>

#include "TMatrixD.h"

/** @brief Exception class for error handling in GENFIT (provides storage for diagnostic information)
 *
 *  @author Christian H&ouml;ppner (Technische Universit&auml;t M&uuml;nchen, original author)
 *  @author Sebastian Neubert  (Technische Universit&auml;t M&uuml;nchen, original author)
 *
 * This is the class that is used for all error handling in GENFIT.
 * It is a utility class that allows to store numbers and matrices together
 * with an error string. The exception class can then be thrown when an error
 * is detected and the C++ exception handling facilities can be used to
 * catch and process the exception.
 */
class GFException : public std::exception
{
 private:

  static bool fQuiet;

  std::string fExcString;
  int fLine;
  std::string fFile;

  std::string fErrorMessage;

  std::string fNumbersLabel;
  std::string fMatricesLabel;
  std::vector<double> fNumbers;
  std::vector<TMatrixD> fMatrices;

  bool fFatal;

 public:
  /** @brief Initializing constructor 
   *
   * @param what error message
   * @param line line at which the exception is created. Can be set through
   * __LINE__ macro
   * @param file sorcefile in which the exception is created. 
   * Can be set through __FILE__ macro
   */
  GFException(std::string, int, std::string);
  virtual ~GFException() throw();
  
  /** @brief set fatal flag. if this is true, the fit stops for this current track repr. */
  void setFatal (bool b=true){fFatal=b;}
  /** @brief get fatal flag. */
  bool isFatal (){return fFatal;}
  /** @brief set list of numbers with description */
  void setNumbers (std::string, const std::vector<double>&);
  /** @brief set list of matrices with description */
  void setMatrices(std::string, const std::vector<TMatrixD>&);

  /** @brief print information in the exception object */
  void info();

  //! standard error message handling for exceptions. use like "std::cerr << e.what();"
  virtual const char* what() const throw();

  std::string getExcString(){return fExcString;}

  static void quiet(bool b=true){fQuiet=b;}

};

#endif

/** @} */
