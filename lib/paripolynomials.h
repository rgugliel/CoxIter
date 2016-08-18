/*
Copyright (C) 2013, 2014, 2015, 2016
Rafael Guglielmetti, rafael.guglielmetti@unifr.ch
*/

/*
This file is part of CoxIter and AlVin.

These are free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

CIVA is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with CIVA. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef __PARI_POLYNOMIALS_H__
#define __PARI_POLYNOMIALS_H__

#include <vector>
#include <pari/pari.h>
#ifdef _USE_LOCAL_GMP_
#include "gmpxx.h"
#else
#include <gmpxx.h>
#endif
using namespace std;

namespace PariPolynomials
{
	/*!	\fn vector2t_POL
	* 	\brief Create a polynomial in PARI form from a C++ vector
	* 
	* 	\param iCoefficients( const vector< int >& ) Coefficients of the vector
	* 
	* 	\return Polynomial (GEN > t_POL)
	*/
	GEN vector2t_POL( const vector< int >& iCoefficients );
	
	/*!	\fn vector2t_POL
	* 	\brief Create a polynomial in PARI form from a C++ vector
	* 
	* 	\param iCoefficients( const vector< mpz_class >& ) Coefficients of the vector
	* 
	* 	\return Polynomial (GEN > t_POL)
	*/
	GEN vector2t_POL( const vector< mpz_class >& iCoefficients );
	
	/*!	\fn t_POL2vector
	* 	\brief Create a C++ vector from a t_POL (PARI)
	* 
	* 	\param poly( GEN ) Polynomial
	* 
	* 	\return Vector (vector< long int >)
	*/
	vector< long int > t_POL2vector( GEN poly );
}

#endif