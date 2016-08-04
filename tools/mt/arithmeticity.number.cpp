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

#include "arithmeticity.number.h"

ArithmeticityNumber::ArithmeticityNumber( )
{
}

ArithmeticityNumber::ArithmeticityNumber( const int& iVal )
: parti( iVal ), part2( false ), part3( false ), part5( false )
{
}

ArithmeticityNumber::ArithmeticityNumber( const ArithmeticityNumber& n )
: parti( n.parti ), part2( n.part2 ), part3( n.part3 ), part5( n.part5 )
{
	update( );
}

ArithmeticityNumber::ArithmeticityNumber( const int& parti, const bool& part2, const bool& part3, const bool& part5 )
: parti( parti ), part2( part2 ), part3( part3 ), part5( part5 )
{
	update( );
}

ArithmeticityNumber ArithmeticityNumber::operator*( ArithmeticityNumber const &n ) const
{
	int n_parti( parti * n.parti );
	bool n_part2( part2 ^ n.part2 ), n_part3( part3 ^ n.part3 ), n_part5( part5 ^ n.part5 );
	
	if( part2 && n.part2 )
		n_parti *= 2;
	
	if( part3 && n.part3 )
		n_parti *= 3;
	
	if( part5 && n.part5 )
		n_parti *= 5;
	
	return ArithmeticityNumber( n_parti, n_part2, n_part3, n_part5 );
}

void ArithmeticityNumber::update( )
{
	if( !parti )
		part2 = part3 = part5 = false;
}

void ArithmeticityNumber::reduceIntergerPart( )
{
	if( parti != 0 )
		parti = 1;
}

ArithmeticityNumber& ArithmeticityNumber::operator=(const int& n)
{
	part2 = part3 = part5 = false;
	parti = n;
	
	return *this;
}

ArithmeticityNumber& ArithmeticityNumber::operator*=( ArithmeticityNumber const &n )
{
	parti *= n.parti;
	
	if( part2 && n.part2 )
	{
		part2 = false;
		parti *= 2;
	}
	else if( part2 || n.part2 )
		part2 = true;
	
	if( part3 && n.part3 )
	{
		part3 = false;
		parti *= 3;
	}
	else if( part3 || n.part3 )
		part3 = true;
	
	if( part5 && n.part5 )
	{
		part5 = false;
		parti *= 5;
	}
	else if( part5 || n.part5 )
		part5 = true;
	
	update( );

	return *this;
}

bool ArithmeticityNumber::operator==( const int& i )
{
	return ( !part2 && !part3 && !part5 && parti == i );
}

bool ArithmeticityNumber::operator==( const ArithmeticityNumber& an )
{
	return ( part2 == an.part2 && part3 == an.part3 && part5 == an.part5 && parti == an.parti );
}

ostream& operator<<( ostream &o, ArithmeticityNumber const &n )
{
	o << "(" << n.parti << "," << ( n.part2 ? "1," : "0," ) << ( n.part3 ? "1," : "0," ) << ( n.part5 ? "1)" : "0)" );

	return o;
}