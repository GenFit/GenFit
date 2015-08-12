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
/** @addtogroup genfit
 * @{
 */

#ifndef genfit_Exception_h
#define genfit_Exception_h

#include <exception>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>

#include <TMatrixD.h>


namespace genfit {

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
class Exception : public std::exception {

 public:
  /** @brief Initializing constructor
   *
   * @param excString error message.
   * @param line line at which the exception is created. Can be set through __LINE__ macro.
   * @param file sourcefile in which the exception is created. Can be set through __FILE__ macro.
   */
  Exception(std::string excString, int line, std::string  file);
  virtual ~Exception() throw();

  //! Set fatal flag.
  void setFatal (bool b=true){fatal_=b;}
  //! Get fatal flag.
  bool isFatal (){return fatal_;}
  //! Set list of numbers with description.
  void setNumbers (std::string, const std::vector<double>&);
  //! Set list of matrices with description.
  void setMatrices(std::string, const std::vector<TMatrixD>&);

  //! Print information in the exception object.
  void info();

  //! Standard error message handling for exceptions. use like "std::cerr << e.what();"
  virtual const char* what() const throw();

  std::string getExcString(){return excString_;}

  //! "std::cerr << e.what();" will not write anything.
  static void quiet(bool b=true){quiet_=b;}

 private:

  static bool quiet_;

  std::string excString_;
  int line_;
  std::string file_;

  std::string errorMessage_;

  std::string numbersLabel_;
  std::string matricesLabel_;
  std::vector<double> numbers_;
  std::vector<TMatrixD> matrices_;

  bool fatal_;

  //ClassDef(Exception,1)

};

} /* End of namespace genfit */
/** @} */

#endif // genfit_Exception_h
