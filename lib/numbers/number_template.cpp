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
 * \file number_rational.cpp
 * \brief Tout ceci ne sert pas vraiment puisque tout est virtuel (mais on peut pas faire du virtuel pur)
 * \author Rafael Guglielmetti
*/

#include "number_template.h"

bool Number_template::operator==(Number_template const &) const
{
	return false;
}

bool Number_template::operator!=(Number_template const &) const
{
	return false;
}

Number_template Number_template::operator+( Number_template const &n ) const
{
	return *this;
}

Number_template& Number_template::operator+=( Number_template const &n )
{
	return *this;
}

Number_template Number_template::operator-( Number_template const &n ) const
{
	return *this;
}

Number_template& Number_template::operator-=( Number_template const &n )
{
	return *this;
}

Number_template Number_template::operator*( Number_template const &n ) const
{
	return *this;
}

Number_template& Number_template::operator*=( Number_template const &n )
{
	return *this;
}

Number_template Number_template::operator/( Number_template const &n ) const
{
	return *this;
}

Number_template&Number_template:: operator/=( Number_template const &n )
{
	return *this;
}

void Number_template::print( ostream & ) const
{
}