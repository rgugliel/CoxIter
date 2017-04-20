/*
Copyright (C) 2013, 2014, 2015, 2016
Rafael Guglielmetti, rafael.guglielmetti@unifr.ch
*/

/*
This file is part of CoxIter.

CoxIter is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

CoxIter is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with CoxIter. If not, see <http://www.gnu.org/licenses/>.
*/

/*!
 * \file signature.h
 * \author Rafael Guglielmetti
 * 
 * \class GrowthRate
 * \brief To compute the signature of a matrix (PARI format)
*/

#ifndef SIGNATURE_H
#define SIGNATURE_H

#include <pari/pari.h>
#include <string>
#include <array>

using namespace std;

class Signature
{
	private:
		GEN gEpsilon; ///< Some small number (typically 10^-50)
		
	public:
		Signature();
		~Signature();
		
		array< unsigned int, 3 > iComputeSignature(string strMatrix);
};

#endif // SIGNATURE_H
