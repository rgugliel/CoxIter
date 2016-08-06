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

#include "regexp.h"

PCRERegexp::PCRERegexp( int _iOvectorSize )
: 	strError( "" ),
	bClassUsed( false ),
	regexp( 0 ),
	regPattern( 0 ),
	regSubject( 0 ),
	regError( 0 ),
	ovector( 0 ),
	iOvectorSize( _iOvectorSize )
{
	// ----------------------------------------------
	// ovector
	if( iOvectorSize <= 3 )
		iOvectorSize = 30;

	ovector = new int[ iOvectorSize ];
}

PCRERegexp::~PCRERegexp( )
{
	delete [ ] ovector;
	ovector = NULL;

	if( bClassUsed )
		pcre_free( regexp ); 
}

string PCRERegexp::get_strError( )
{
	return strError;
}

int PCRERegexp::preg_match_all( const string &pattern, const string &subject, PCREResult &results, const int& optionsCompile )
{
	int rc;
	int options( 0 );
	
	regPattern = pattern.c_str( );
	regSubject = subject.c_str( );
	regSubjectLength = subject.size( );

	// ----------------------------------------------------
	// compilation of the regexp
	regexp = pcre_compile( regPattern, optionsCompile, &regError, &regErrorOffset, NULL );
	if( regexp == NULL )
	{
		throw( string( "Error: PCRE_COMPILATION_ERROR" ) );
		strError = "PCRE_COMPILATION_ERROR";
		return -1;
	}
	
	bClassUsed = true;
	results.clear( );

	// ----------------------------------------------------
	// execution of the regexp
	options = 0;
	
	rc = pcre_exec( regexp, NULL, regSubject, regSubjectLength, 0, options, ovector, iOvectorSize );
	if( rc < 0 ) // no result
		return 0;
	else if( rc == 0 ) // if the vector is too small
	{
		strError = "PCRE_OVECTOR_TOO_SMALL";
		return -1;
	}
	
	// ----------------------------------------------------
	// first result
	for( int i = 0; i < rc; i++ )
	{
		results.push_back( vector<string>( ) );
		results[i].push_back( subject.substr( ovector[2*i], ovector[2*i+1] - ovector[2*i] ) );
	}

	// ----------------------------------------------------
	// next results
	int start_offset;
	unsigned int match_count( 1 );
	while( true )
	{
		options = 0;
		start_offset = ovector[1]; // end of the last occurrence

		if( ovector[0] == ovector[1] )
		{
			if( ovector[0] == regSubjectLength )
				break;
			
			options = PCRE_NOTEMPTY | PCRE_ANCHORED;
		}

		rc = pcre_exec( regexp, NULL, regSubject, regSubjectLength, start_offset, options, ovector, iOvectorSize );

		if( rc == PCRE_ERROR_NOMATCH )
		{
			if( options == 0 )
				break;

			ovector[1] = start_offset + 1;
			continue;
		}

		if( rc < 0 )
		{
			strError = "UNKNOWN_ERROR";
			pcre_free( regexp );
			return -1;
		}
		else if(rc == 0)
		{
			strError = "PCRE_OVECTOR_TOO_SMALL";
			return -1;
		}

		match_count++;
		for( int i = 0; i < rc; i++ )
			results[i].push_back( subject.substr( ovector[2*i], ovector[2*i+1] - ovector[2*i] ) );
	}
	
	/*
	cout << "Match count: " << match_count << " / size: " << results.size() << endl;
	for( auto row : results )
	{
		cout << "New: " << endl;
		for( auto str : row )
			cout << "\t" << str << endl;
	}*/
	
	return match_count;
}
