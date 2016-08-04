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

#include "rational.h"

template <typename T>   
Rational<T>::Rational( ) 
	: a( 0 ), b( 1 )
{
	update( );
}

template <typename T>
Rational<T>::Rational( const int& i )
	: a( i ), b( 1 )
{
	update( );
}

template <typename T>
Rational<T>::Rational( T a ) 
	: a( a ), b( 1 )
{
	update( );
}

template <typename T>
Rational<T>::Rational( T a, T b ) 
	: a( a ), b( b )
{
	update( );
}

template <typename T>
void Rational<T>::update( )
{
	if( a.divideByIfDivisible( &b ) )
	{
		b = 1;
	}
	else
	{
		T iGCD( a );
		iGCD.gcd( &b );
		
		a /= iGCD;
		b /= iGCD;
	}

	if( b.bIsLessThan( 0 ) )
	{
		a.multiplyBy( -1 );
		b.multiplyBy( -1 );
	}

	isZero = a.bIsEqualTo( 0 );
	hasDenominatorOne = b.bIsEqualTo( 1 );
	isOne = hasDenominatorOne && a.bIsEqualTo( 1 );
	isMinusOne = hasDenominatorOne && a.bIsEqualTo( -1 );
}

template <typename T>
Rational<T> Rational<T>::operator+( Rational<T> const &n ) const
{
	return Rational( a * n.b + b * n.a, b * n.b );
}

template <typename T>
Rational<T>& Rational<T>::operator+=( Rational<T> const &n )
{
	a = a * n.b + b * n.a ;
	b = b * n.b ;
	update( );

	return *this;
}

template <typename T>
Rational<T> Rational<T>::operator-( Rational<T> const &n ) const
{
	return Rational( a * n.b - b * n.a, b * n.b );
}

template <typename T>
Rational<T> Rational<T>::operator-( ) const
{
	return Rational( -a, b );
}

template <typename T>
void Rational<T>::opp( Rational<T>* &_c ) const
{
	_c = new Rational( -a, -b );
}

template <typename T>
Rational<T>& Rational<T>::operator-=( Rational<T> const &n )
{
	a = a * n.b - b * n.a ;
	b = b * n.b ;
	update( );

	return *this;
}

template <typename T>
Rational<T> Rational<T>::operator*( Rational<T> const &n ) const
{
	return Rational( a * n.a, b * n.b );
}

template <typename T>
Rational<T>& Rational<T>::operator*=( Rational<T> const &n )
{
	a *= n.a;
	b *= n.b;
	update( );

	return *this;
}

template <typename T>
Rational<T> Rational<T>::operator/( Rational<T> const &n ) const
{
	if( n.isZero )
		throw( 0 );

	return Rational( a * n.b, b * n.a );
}

template <typename T>
Rational<T>& Rational<T>::operator/=( Rational<T> const &n )
{
	if( n.isZero )
		throw( 0 );
	
	a *= n.b;
	b *= n.a;
	update( );

	return *this;
}

template <typename T>
Rational<T>& Rational<T>::operator=( long int i )
{
	a = i;
	b = 1;
	
	isZero = a.bIsEqualTo( 0 );
	hasDenominatorOne = b.bIsEqualTo( 1 );
	isOne = hasDenominatorOne && a.bIsEqualTo( 1 );
	isMinusOne = hasDenominatorOne && a.bIsEqualTo( -1 );
	
	return *this;
}

template <typename T>
bool Rational<T>::operator>(const Rational<T>& r) const
{
	Rational<T> n( *this - r );
	return ( n.a > 0 );
}

template <typename T>
bool Rational<T>::operator>=(const int& ni) const
{
	Rational<T> n( *this - ni );
	return ( n.a.bIsGreaterOEThan(0) );
}

template <typename T>
bool Rational<T>::operator==( const int& i ) const
{
	return ( a == i && b == 1 );
}

template <typename T>
bool Rational<T>::operator==( Rational<T> const& n ) const
{
	return ( a == n.a && b == n.b );
}

template <typename T>
bool Rational<T>::operator!=( Rational<T> const& n ) const
{
	return ( a != n.a || b != n.b );
}

template <typename T>
void Rational<T>::print( ostream &o ) const
{
	if( hasDenominatorOne )
		o << a;
	else
		o << a << "/" << b;
}

template <typename T>
string Rational<T>::to_string() const
{
	if( hasDenominatorOne )
		return a.get_str( );
	else
		return ( a.get_str() + "/" + b.get_str() );
}

template <typename T>
ostream& operator<<( ostream &o, Rational<T> const &n )
{
	n.print( o );

	return o;
}

template <typename T>
bool Rational<T>::get_hasDenominatorOne( ) const
{
	return hasDenominatorOne;
}

/*
BigInteger gcd( const BigInteger &n, const BigInteger &m )
{
	return ( m == 0 ? n : gcd( m, n % m ) );
}*/
