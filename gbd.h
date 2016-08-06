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
 * \file gbd.h
 * \author Rafael Guglielmetti
 * 
 * \class GBD
 * \brief Gal, Bonnafé-Dyer
 * This part is still experimental!
*/

#ifndef GBD_H
#define GBD_H

#include <string>

using namespace std;

#include "coxiter.h"

struct NewVertex
{
	string strLabel;
	unsigned int iIndex;
	unsigned int iOriginVertex;
};

class GBD
{
	private:
		CoxIter* ci;
		vector< vector< unsigned int > > iCox; ///< Coxeter matrix
		unsigned int iVerticesCount; ///< Number of vertices of the starting graph
		
		unsigned int iVertex; ///< Index of the vertex
		
		vector< NewVertex > newVertices;
		vector< vector< unsigned int > > iNewCox; ///< New Coxeter matrix
		unsigned int iNewVerticesCount;
		
		string strError;
		
	public:
		GBD( CoxIter* ci );
		
		bool removeVertex( const string& strVertexName );
		void printMatrix( vector< vector< unsigned int > >* iMatrix );
		
		string get_strError( ) const;
};

#endif // GBD_H
