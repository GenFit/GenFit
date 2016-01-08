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

#include "Exception.h"
#include "IO.h"

namespace genfit {

bool Exception::quiet_ = false;

Exception::Exception(std::string excString, int line, std::string  file) :
    excString_(excString), line_(line), file_(file), fatal_(false) {
  std::ostringstream ErrMsgStream;
  ErrMsgStream << "genfit::Exception thrown with excString:"
         << std::endl << excString_ << std::endl
         << "in line: " << line_ << " in file: " << file_ << std::endl
         << "with fatal flag " << fatal_ << std::endl;
  errorMessage_ = ErrMsgStream.str();
}

Exception::~Exception() throw() {
}

void Exception::setNumbers(std::string _numbersLabel,
         const std::vector<double>& _numbers) {
  numbersLabel_ = _numbersLabel;
  numbers_ = _numbers;
}

void Exception::setMatrices(std::string _matricesLabel,
          const std::vector<TMatrixD>& _matrices) {
  matricesLabel_ = _matricesLabel;
  matrices_ = _matrices;
}

const char* Exception::what() const throw(){
  return errorMessage_.c_str();
}

void Exception::info() {
  if(quiet_) return;
  if(numbers_.empty() && matrices_.empty()) return; //do nothing
  debugOut << "genfit::Exception Info Output" << std::endl;
  debugOut << "===========================" << std::endl;
  if(numbersLabel_ != "") {
  debugOut << "Numbers Label String:" << std::endl;
  debugOut << numbersLabel_ << std::endl;
  }
  if(!numbers_.empty()) {
  debugOut << "---------------------------" << std::endl;
  debugOut << "Numbers:" << std::endl;
  for(unsigned int i=0;i<numbers_.size(); ++i ) debugOut << numbers_[i] << std::endl;
  }
  if(matricesLabel_ != "") {
  debugOut << "---------------------------" << std::endl;
  debugOut << "Matrices Label String:" << std::endl;
  debugOut << matricesLabel_ << std::endl;
  }
  if(!matrices_.empty()) {
  debugOut << "---------------------------" << std::endl;
  debugOut << "Matrices:" << std::endl;
  for(unsigned int i=0;i<matrices_.size(); ++i ) matrices_[i].Print();
  }
  debugOut << "===========================" << std::endl;
}

} /* End of namespace genfit */

