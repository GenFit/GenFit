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
#include "GFException.h"

bool GFException::fQuiet = false;

GFException::GFException(std::string _excString, int _line, std::string  _file) : fExcString(_excString), fLine(_line), fFile(_file),fFatal(false) {
  std::ostringstream ErrMsgStream;
  ErrMsgStream << "GFException thrown with excString:"
	       << std::endl << fExcString << std::endl
	       << "in line: " << fLine << " in file: " << fFile << std::endl
	       << "with fatal flag " << fFatal << std::endl;
  fErrorMessage = ErrMsgStream.str();
}

GFException::~GFException() throw() {
}

void GFException::setNumbers(std::string _numbersLabel,
				 const std::vector<double>& _numbers) {
  fNumbersLabel = _numbersLabel;
  fNumbers = _numbers;
}

void GFException::setMatrices(std::string _matricesLabel,
				  const std::vector<TMatrixD>& _matrices) {
  fMatricesLabel = _matricesLabel;
  fMatrices = _matrices;
}

const char* GFException::what() const throw(){
  if(fQuiet) return "";
  return fErrorMessage.c_str();
}

void GFException::info() {
  if(fQuiet) return;
  if(fNumbers.empty() && fMatrices.empty()) return; //do nothing
  std::cout << "GFException Info Output" << std::endl;
  std::cout << "===========================" << std::endl;
  if(fNumbersLabel != "") {
	std::cout << "Numbers Label String:" << std::endl;
	std::cout << fNumbersLabel << std::endl;
  }
  if(!fNumbers.empty()) {
	std::cout << "---------------------------" << std::endl;
	std::cout << "Numbers:" << std::endl;
	for(unsigned int i=0;i<fNumbers.size(); ++i ) std::cout << fNumbers[i] << std::endl;
  }
  if(fMatricesLabel != "") {
	std::cout << "---------------------------" << std::endl;
	std::cout << "Matrices Label String:" << std::endl;
	std::cout << fMatricesLabel << std::endl;
  }
  if(!fMatrices.empty()) {
	std::cout << "---------------------------" << std::endl;
	std::cout << "Matrices:" << std::endl;
	for(unsigned int i=0;i<fMatrices.size(); ++i ) fMatrices[i].Print();
  }
  std::cout << "===========================" << std::endl;  
}
