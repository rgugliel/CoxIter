/*
Copyright (C) 2013, 2014, 2016, 2016
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
 * \file rational.h
 * \brief Rational numbers (template)
 * \author Rafael Guglielmetti
 * \class Rational rational.h
*/

#ifndef __RATIONAL_GENERIC_H__
#define __RATIONAL_GENERIC_H__

#include <iostream>

using namespace std;

template <typename T>
class Rational 
{
	public:
		T a;
		T b;
		
	private:
		// TODO init
		bool isZero;
		bool hasDenominatorOne;
		bool isOne;
		bool isMinusOne;

	private:
		/*! \fn update
		 * 	\brief Met à jour les attributs (gcd, isInt, ...)
		 */
		void update( );

	public:
		Rational( );
		
		Rational( T a, T b );
		Rational( T a );
		
		Rational( const int& i );

		bool operator>( Rational const& ) const;
		bool operator>=( int const& ) const;
		bool operator==( int const& ) const;
		bool operator==( Rational const& ) const;
		bool operator!=( Rational const& ) const;
		
		Rational& operator=( long int );
		
		Rational operator+( Rational const &n ) const;
		Rational& operator+=( Rational const &n );

		Rational operator-( Rational const &n ) const;
		Rational operator-( ) const;
		void opp( Rational* &_c ) const;
		Rational& operator-=( Rational const &n );
		
		Rational operator*( Rational const &n ) const;
		Rational& operator*=( Rational const &n );

		Rational operator/( Rational const &n ) const;
		Rational& operator/=( Rational const &n );

		void print( ostream & ) const;	
		
		string to_string( ) const;
		
		bool get_hasDenominatorOne() const;
};

template <typename T>
ostream& operator<<( ostream& , Rational<T> const & );

template <typename T>
T abs( const T& r )
{
	T rabs( r );
	if( rabs.a.bIsLessThan( 0 ) )
		rabs.a.multiplyBy( -1 );
	
	return rabs;
}

#include "rational.tpp"

#endif