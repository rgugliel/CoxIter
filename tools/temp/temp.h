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

#ifndef TEMP_H
#define TEMP_H

#include <string>
#include <fstream>
#include <vector>
#include <map>
#include <iostream>

#include "../../../tools/regexp.h"
#include "../../coxiter.h"
#include "../../arithmeticity.h"
#include "../../growthrate.h"
#include "../../../tools/string.h"

using namespace std;

class Temp
{
	private:
		string strError;
		vector< string > strGraphsFilename;
		
		string strPARItoMathematica( string str );
	
	public:
		bool bReadGraphsList( const string& strList  );
		vector< string > get_strGraphsFilename( ) const;
		
		void exportGraphs( );
		void dropEnds( const string& strFolderOutput );
		
		void compareWithTerragni( );
		
		void testArithmeticity( );
};

#endif // TEMP_H
