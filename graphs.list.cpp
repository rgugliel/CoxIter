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

#include "graphs.list.h"

GraphsList::GraphsList( size_t iVerticesCount, vector< string > *ptr_map_vertices_indexToLabel ) 
: iVerticesCount( iVerticesCount ), graphsCount( vector< size_t >( iVerticesCount + 1, 0 ) )
{
	for( unsigned int i(0); i <= iVerticesCount; ++i )
		graphs.push_back( GraphsListN( i, ptr_map_vertices_indexToLabel ) );
	
	iGraphsCount = 0;
}

void GraphsList::addGraph( const std::vector< short unsigned int >& iVertices, const std::vector< bool >& bVerticesLinkable, const unsigned int& iType, bool bSpherical, const unsigned int &iVertexSupp1, const unsigned int &iVertexSupp2, const unsigned int &iDataSupp )
{
	size_t sizeTemp, iVCount( iVertices.size( ) );
	
	if( bSpherical )
	{
		if( iType == 1 || iType == 3 || iType == 4 || iType == 5 || iType == 7 )
			iVCount++;
	}
	else
	{
		if( ( iType == 1 && iVCount >= 3 && iDataSupp ) || iType == 2 || ( iType == 3 && iVCount >= 3 ) || ( iType == 4 && iVCount == 5 )) // TBn, TCn, TE6
			iVCount += 2;
		else if( ( iType == 0 && iDataSupp ) || iType == 1 || iType == 3 || iType == 4 || iType == 5 || ( iType == 6 && iDataSupp ) ) // TAn ou TB3 ou TD4 ou TEn ou TG_2
			iVCount++;
	}
	
	sizeTemp = graphs[ iVCount ].size( );
	graphs[ iVCount ].addGraph( iVertices, bVerticesLinkable, iType, bSpherical, iVertexSupp1, iVertexSupp2, iDataSupp );
	if( sizeTemp != graphs[ iVCount ].size( ) )
	{
		iGraphsCount++;
		graphsCount[ iVCount ]++;
	}
}

Graph* GraphsList::begin( )
{
	if( !iGraphsCount )
		return 0;
	
	for( unsigned int i(1); i <= iVerticesCount; i++ )
	{
		if( graphsCount[i] )
			return graphs[i].begin( );
	}
	
	return 0;
}

Graph* GraphsList::next( size_t &iVCount, size_t &iGraphIndex )
{
	if( ( ++iGraphIndex ) < graphsCount[ iVCount ] )
		return graphs[ iVCount ].next( iGraphIndex );
	else
	{
		for( size_t i( iVCount + 1 ); i <= iVerticesCount; i++ )
		{
			if( graphsCount[i] )
			{
				iVCount = i;
				iGraphIndex = 0;
				return graphs[ iVCount ].next( iGraphIndex );
			}
		}
	}
	
	return 0;
}

ostream& operator<<( ostream &o, const GraphsList &g )
{
	unsigned int iMax( g.graphs.size( ) );
	for( unsigned int i(0); i < iMax; ++i )
		o << g.graphs[i];
	
	return o;
}