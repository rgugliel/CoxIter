/*
Copyright (C) 2013, 2014
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

#ifndef ARITHMETICITY_ADMISSIBLE_NUMBER_H
#define ARITHMETICITY_ADMISSIBLE_NUMBER_H

#include <iostream>

using namespace std;

/*!
 * \file arithmeticity.number.h
 * \author Rafael Guglielmetti
 * 
 * \class ArithmeticityNumber
 * \brief A number of the form: a * a2 * a3, with a an integer, a2 = 1 or a2 = sqrt(2), a3 = 1 or a3 = sqrt(3)
*/

class ArithmeticityNumber
{
	private:
		void update( );
		
	public:
		ArithmeticityNumber( );
		ArithmeticityNumber( const int& );
		ArithmeticityNumber( const ArithmeticityNumber & );
		ArithmeticityNumber( const int&, const bool&, const bool&, const bool& );
		
		int parti;
		bool part2;
		bool part3;
		bool part5;
		
		ArithmeticityNumber operator*( ArithmeticityNumber const &n ) const;
		ArithmeticityNumber& operator*=( ArithmeticityNumber const &n );
		ArithmeticityNumber& operator=( const int &n );
		
		void reduceIntergerPart( );
		
		bool operator==( const int& i );
		bool operator==( const ArithmeticityNumber& an );
};

ostream& operator<<( ostream& , ArithmeticityNumber const& );

#endif // ARITHMETICITY_ADMISSIBLE_NUMBER_H
