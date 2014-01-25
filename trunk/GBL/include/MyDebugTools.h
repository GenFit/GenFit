#ifndef MYDEBUGTOOLS_H
#define MYDEBUGTOOLS_H


#include <map>
#include <iostream>
#include <iomanip>
#include <TMatrixD.h>
#include <assert.h>
#include <sstream>
#include <TMath.h>
#include <TVector3.h>

/*
using namespace std;

ofstream file("debug.txt");



void outputMatrix(TMatrixDSym matrix, std::string caption = "")
{
	if (caption != "") file << caption << endl;
	for (int i=0; i<matrix.GetNrows(); i++)
	{
		for (int j=0; j<matrix.GetNcols(); j++)
		{
			file << setw(12) << matrix[i][j] << "  ";
		}
		file << endl;
	}
}

void outputMatrix(TMatrixD matrix, std::string caption = "")
{
	if (caption != "") file << caption << endl;
	for (int i=0; i<matrix.GetNrows(); i++)
	{
		for (int j=0; j<matrix.GetNcols(); j++)
		{
			file << setw(12) << matrix[i][j] << "  ";
		}
		file << endl;
	}
}

void outputVector(TVectorD vector, std::string caption)
{
	if (caption != "") file << caption << endl;
	for (int i=0; i<vector.GetNoElements(); i++)
	{
		file << setw(12) << vector[i] << "  ";
	}
	file << endl;
}

void outputVector(TVector3 vector, std::string caption)
{
	if (caption != "") file << caption << endl;
	for (int i=0; i<3; i++)
	{
		file << setw(12) << vector[i] << "  ";
	}
	file << endl;
}
*/
#endif
