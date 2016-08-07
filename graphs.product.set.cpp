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

#include "graphs.product.set.h"

GraphsProductSet::GraphsProductSet( )
: iRank( 0 )
{
}

GraphsProductSet::GraphsProductSet( const GraphsProduct& gp )
{
	iRank = gp.iRank;
	
	for( vector< Graph* >::const_iterator it( gp.graphs.begin( ) ); it != gp.graphs.end( ); ++it )
		graphs.insert( *it );
}

vector< short unsigned int > GraphsProductSet::get_iVertices( ) const
{
	vector< short unsigned int > iVertices;
	
	for( set< Graph* >::const_iterator itProd( graphs.begin( ) ); itProd != graphs.end( ); ++itProd )
	{
		for( vector< short unsigned int >::const_iterator itV( (*itProd)->iVertices.begin( ) ); itV != (*itProd)->iVertices.end( ); ++itV )
		{
			iVertices.push_back( *itV );
		}
	}
	
	return iVertices;
}


bool GraphsProductSet::b_areVerticesSubsetOf( const GraphsProductSet& gp ) const
{
	/* NOTE:
		* We could sort de vectors and use the standard function includes but this is slower (tested with 17-vinb85).
	*/
	vector< short unsigned int > iVerticesSmall( get_iVertices( ) ), iVerticesBig( gp.get_iVertices( ) );
	
	for( vector< short unsigned int >::const_iterator it( iVerticesSmall.begin( ) ); it != iVerticesSmall.end( ); ++it )
	{
		if( find( iVerticesBig.begin( ), iVerticesBig.end( ), *it ) == iVerticesBig.end( ) )
			return false;
	}
	
	return true;
}

ostream& operator<<(ostream &o, const GraphsProductSet& gp)
{
	for( set< Graph* >::const_iterator it( gp.graphs.begin( ) ); it != gp.graphs.end( ); ++it )
		o << **it;
	
	return o;
}