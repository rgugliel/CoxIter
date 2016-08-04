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
 * \file growthrate.h
 * \author Rafael Guglielmetti
 * 
 * \class GrowthRate
 * \brief To compute the growth rate
*/

#ifndef GROWTHRATE_H
#define GROWTHRATE_H

#include <vector>
#include <pari/pari.h>
#include <gmpxx.h>
#include <iostream>

#include "../tools/paripolynomials.h"

using namespace std;
using namespace PariPolynomials;

struct GrowthRate_Result
{
	int iPerron;
	int iPisot;
	int iSalem;
	string strGrowthRate;
	bool bComputed;
};

class GrowthRate
{
	private:
		vector< GEN > t_POLfactors; ///< Irreducible factors of the denominator of the growth series (only those having a root between 0 and 1)
		
		GEN gGrowthRate; ///< Maximal positive root of the polynomial
		GEN gMaximalRoots; ///< Roots of the polynomial which has the maximal root
		long int iIndexMaximalRoot; ///< Factor which contains the minimal root // TODO: utile?
		
		GEN gEpsilon; ///< Some small number (typically 10^-50)
		long int iPariPrecision; ///< Given as prec (typically 8)
		
	public:
		GrowthRate();
		~GrowthRate();
		
		GrowthRate_Result grrComputations( vector< mpz_class > iPolynomial, const bool& bOnlyGrowthRate = false );
		
	private:
		/*!	\fn irreducibleFactors( const vector< mpz_class >& iPolynomial )
		 * 	Factor the polynomial iPolynomial and store the factors into t_POLfactors
		 * 	\param iPolynomial( vector< mpz_class > ) The polynomial (coefficients in GMPlib)
		 */
		void irreducibleFactors( const vector< mpz_class >& iPolynomial );
		
		/*!	\fn minimalRoot()
		 * 	Find which fact has the smallest (positive, <1) real root and store the index into iIndexMinimalRoot
		 */
		void minimalRoot();

		/*!	\fn iNumberRootsUnitCircle( GEN gPol );
		 * 	For a palindromic polynomial, try to compute the number of zeros on the unit circle
		 * 	\param gPol (GEN, PARI polynomial)
		 * 	\return Number of roots on the unit circle of -1 if we cannot decide
		 */
		long int iNumberRootsUnitCircle( GEN gPol );
};

#endif // GROWTHRATE_H
