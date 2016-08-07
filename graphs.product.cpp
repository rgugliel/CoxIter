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

#include "graphs.product.h"

GraphsProduct::GraphsProduct( )
: iRank( 0 )
{
}

vector< vector< short unsigned int > > GraphsProduct::createFootPrint( )
{
	vector< vector< short unsigned int > > iGraphsCountByType( 8 );
	
	size_t i, j, iMax, iIndex, iGraphsCount( graphs.size( ) );

	// on compte les graphes qui apparaissent
	for( i = 0; i < iGraphsCount; i++ )
	{
		iIndex = graphs[i]->iGraphType;
		
		// Si graphe de type G, ce qui compte c'est le poids
		iMax = ( iIndex == 6 && graphs[i]->bSpherical ) ? graphs[i]->iDataSupp : graphs[i]->iVertices.size( );
		
		for( j = iGraphsCountByType[ iIndex ].size( ); j < iMax; j++ )
			iGraphsCountByType[ iIndex ].push_back( 0 );
		
		iGraphsCountByType[ iIndex ][ iMax - 1 ]++;
	}

	return iGraphsCountByType;
}

ostream& operator<<(ostream &o, const GraphsProduct& gp)
{
	for( vector< Graph* >::const_iterator it( gp.graphs.begin( ) ); it != gp.graphs.end( ); ++it )
		o << **it;
	
	return o;
}