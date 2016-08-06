/*
Copyright (C) 2013, 2014
Rafael Guglielmetti, rafael.guglielmetti@unifr.ch
*/

/*
This file is part of CoxIter and AlVin.

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
 * \file regexp.h
 * \brief To simplify the use of regexp
 * \author Rafael Guglielmetti
 * \class PCRERegexp regexp.h
*/

#include <string>
#include <iostream>
#include <vector>

#include <pcre.h>

using namespace std;

#ifndef __REGEXP_H__
#define __REGEXP_H__

#include <stdexcept>      // std::invalid_argument // TODO: remove

typedef vector< vector< string > > PCREResult;

class PCRERegexp
{
	private:
		string strError; ///< Eventually, error code
		bool bClassUsed; ///< True if the class was used

		pcre *regexp;
		const char *regPattern;
		const char *regSubject;
		const char *regError;
		
		int regErrorOffset;
		int regSubjectLength;

		int *ovector; ///< Information about the result
		int iOvectorSize; ///< Size of the ovector array

	public:
		PCRERegexp( int iOvectorSize = 30 );
		~PCRERegexp( );

		/*!	\fn preg_match_all
		* 	\brief As the usual PHP preg_match_all: executes the regexp and fetch all the occurrences
		* 	\param pattern( const string & ) The pattern to search for
		* 	\param subject( const string & ) The input string
		* 	\param results( PCREResult & ) Array of all matches in multi-dimensional array
		* 	\param optionsCompile( const int& ) Options for the regexp. For example: 0, PCRE_CASELESS
		* 	\return Number of results (unsigned int) or -1 if an error occurred
		*/
		int preg_match_all( const string &pattern, const string &subject, PCREResult &results, const int& optionsCompile = 0 );

		/*! \fn get_strError
		* 	\brief Return the error ccode
		* 	\return Error code (string)
		*/
		string get_strError( );
};

#endif
