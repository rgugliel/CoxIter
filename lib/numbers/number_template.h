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

/*!
 * \file number_template.h
 * \brief Permet d'avoir un type générique (regroupant Number_rational, ...)
 * \author Rafael Guglielmetti
 * \class Number_template number_template.h
*/

#ifndef __NUMBER_TEMPLATE_H__
#define __NUMBER_TEMPLATE_H__

#include <iostream>

using namespace std;

class Number_template
{
	public:
		bool isZero;
		bool isMinusOne;
		bool isOne;
		bool isInt;

	public:
		virtual bool operator==(Number_template const &) const;
		virtual bool operator!=(Number_template const &) const;
		
		virtual Number_template operator+( Number_template const &n ) const;
		virtual Number_template& operator+=( Number_template const &n );

		virtual Number_template operator-( Number_template const &n ) const;
		virtual Number_template& operator-=( Number_template const &n );
		
		virtual Number_template operator*( Number_template const &n ) const;
		virtual Number_template& operator*=( Number_template const &n );

		virtual Number_template operator/( Number_template const &n ) const;
		virtual Number_template& operator/=( Number_template const &n );

		virtual void print( ostream & ) const;
};

ostream& operator<<( ostream& , Number_template const & );

#endif
