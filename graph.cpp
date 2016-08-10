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

#include "graph.h"

Graph::Graph( const vector< short unsigned int >& iVertices, vector< string > *ptr_map_vertices_indexToLabel, const vector< bool >& bVerticesLinkable, const unsigned int &iType, const bool& bSpherical, const unsigned int &iDataSupp )
: iGraphType( iType ), iVertices( iVertices ), bVerticesLinkable( bVerticesLinkable ), iDataSupp( iDataSupp ), bSpherical( bSpherical ), ptr_map_vertices_indexToLabel( ptr_map_vertices_indexToLabel ), b_map_vertices_indexToLabelIsEmpty( !ptr_map_vertices_indexToLabel || ptr_map_vertices_indexToLabel->size( ) == 0 )
{
}

ostream& operator<<( ostream &o, const Graph &g )
{
	unsigned int iMax( g.iVertices.size( ) ), i;
	
	// -----------------------------------------------------------------------
	// Nom du graphe
	o << "\t\t" << ( g.bSpherical ? "" : "T" ) << (char)( g.iGraphType + 65 ) << ( g.bSpherical ? iMax : iMax - 1 ) << " ; ";

	// -----------------------------------------------------------------------
	// sommets qui constituent le graphe
	for( i = 0; i < iMax; i++ )
	{
		if( g.b_map_vertices_indexToLabelIsEmpty )
			o << ( g.iVertices[i] + 1 ) << " ";
		else
			o << ((*g.ptr_map_vertices_indexToLabel)[g.iVertices[i]]) << " ";

		if( ( ( g.iGraphType == 3 || g.iGraphType == 4 ) && i == ( iMax - 2 ) ) || ( ( g.iGraphType == 1 && !g.bSpherical ) && ( i == ( iMax - 2 ) || i == 0 ) ) ) // ( Dn ou En ) ou ( \tilde Bn )
			o << "| ";
	}
	
	// -----------------------------------------------------------------------
	// Poids (pour G_2^k)
	if( g.iDataSupp && g.bSpherical )
		o << " (" << g.iDataSupp << ")";
	
	o << endl;
	
	return o;
}

bool Graph::bIsSubgraphOf( const Graph* grBig ) const
{
	if( *this == *grBig )
		return true;
	
	if( bSpherical && grBig->bSpherical )
		return bIsSubgraphOf_spherical_spherical( grBig );
	else if( bSpherical && !grBig->bSpherical )
		return bIsSubgraphOf_spherical_euclidean( grBig );
	else // We do not test other cases
		throw( 0 );
}

bool Graph::bIsSubgraphOf_spherical_euclidean( const Graph* grBig ) const
{
	size_t iVerticesCount( iVertices.size( ) ), iVerticesBigCount( grBig->iVertices.size( ) );
	
	if( iVerticesCount == 1 )
		return ( find( grBig->iVertices.begin( ), grBig->iVertices.end( ), iVertices[0] ) != grBig->iVertices.end( ) );
	
	// ---------------------------------------------------------------------
	// A_n < TA_m
	if( iGraphType == 0 && grBig->iGraphType == 0 && iVerticesBigCount > 2 )
	{
		vector< short unsigned int >::const_iterator itBig( find( grBig->iVertices.begin( ), grBig->iVertices.end( ), iVertices[0] ) );
		if( itBig == grBig->iVertices.end( ) )
			return false;
		
		int iDir;
		
		if( ( itBig + 1 ) != grBig->iVertices.end( ) && *( itBig + 1 ) == iVertices[1] )
			iDir = 1;
		else if( ( itBig + 1 ) == grBig->iVertices.end( ) && grBig->iVertices[0] == iVertices[1] )
			iDir = 1;
		else if( ( itBig - 1 ) >= grBig->iVertices.begin( ) && *( itBig - 1 ) == iVertices[1] )
			iDir = -1;
		else if( ( itBig - 1 ) < grBig->iVertices.begin( ) && grBig->iVertices[ iVerticesBigCount - 1 ] == iVertices[1] )
			iDir = -1;
		else
			return false;
		
		for( vector< short unsigned int >::const_iterator itSub( iVertices.begin( ) ); itSub != iVertices.end( ); ++itSub )
		{
			if( *itSub != *itBig )
				return false;
			
			itBig += iDir;
			if( itBig == grBig->iVertices.end( ) )
				itBig = grBig->iVertices.begin( ); // move to the beginning
			else if( itBig < grBig->iVertices.begin( ) )
				itBig = grBig->iVertices.end( ) - 1; // move to the end
		}
		
		return true;
	}
	
	// ---------------------------------------------------------------------
	// A_n < TB_m
	if( iGraphType == 0 && grBig->iGraphType == 1 )
	{
		if( iVerticesBigCount == 4 )
		{
			/* In this case, if we remove the first edge, the obtained A3 is not correct */
			vector<short unsigned int> iVTemp;
			iVTemp.push_back( grBig->iVertices[2] );
			iVTemp.push_back( grBig->iVertices[1] );
			iVTemp.push_back( grBig->iVertices[3] );
			
			return bAnSubAm( iVertices, iVTemp );
		}
		else
		{
			Graph g( vector<short unsigned int>( grBig->iVertices.begin( ) + 1, grBig->iVertices.end( ) ), 0, vector< bool >( 0 ), 3, true, 0 );
			
			return bIsSubgraphOf_spherical_spherical( &g );
		}
	}
	
	// ---------------------------------------------------------------------
	// A_n < TC_m
	if( iGraphType == 0 && grBig->iGraphType == 2 )
	{
		return bAnSubAm( iVertices, vector<short unsigned int>( grBig->iVertices.begin( ) + 1, grBig->iVertices.end( ) - 1 ) );
	}
	
	// ---------------------------------------------------------------------
	// A_n < TD_m
	if( iGraphType == 0 && grBig->iGraphType == 3 )
	{
		if( iVerticesCount <= 3 ) // could be in one of the ends
		{
			vector<short unsigned int> iVTemp;
			
			// left hand side
			iVTemp.push_back( grBig->iVertices[0] );
			iVTemp.push_back( grBig->iVertices[2] );
			iVTemp.push_back( grBig->iVertices[1] );
			if( bAnSubAm( iVertices, iVTemp ) )
				return true;
			iVTemp.clear( );
			
			// right hand side
			iVTemp.push_back( grBig->iVertices[ iVerticesBigCount - 2 ] );
			iVTemp.push_back( grBig->iVertices[ iVerticesBigCount - 3 ] );
			iVTemp.push_back( grBig->iVertices[ iVerticesBigCount - 1 ] );
			if( bAnSubAm( iVertices, iVTemp ) )
				return true;
		}
		
		vector<short unsigned int>::const_iterator itStart( find( grBig->iVertices.begin( ), grBig->iVertices.end( ), iVertices[0] ) );
		vector<short unsigned int>::const_iterator itEnd( find( grBig->iVertices.begin( ), grBig->iVertices.end( ), iVertices[ iVerticesCount - 1 ] ) );
		
		vector<short unsigned int> iVTemp( grBig->iVertices.begin( ) + 2, grBig->iVertices.end( ) - 2 ); // basis (middle)
		
		if( itStart == grBig->iVertices.begin( ) )
			iVTemp.insert( iVTemp.begin( ), grBig->iVertices[0] );
		else if( itStart == ( grBig->iVertices.begin( ) + 1 ) )
			iVTemp.insert( iVTemp.begin( ), grBig->iVertices[1] );
		
		if( itEnd == grBig->iVertices.end( ) )
			iVTemp.push_back( grBig->iVertices[ iVerticesBigCount - 1 ] );
		else if( itEnd == ( grBig->iVertices.end( ) - 1 ) )
			iVTemp.push_back( grBig->iVertices[ iVerticesBigCount - 2 ] );
		
		return bAnSubAm( iVertices, iVTemp );
	}
	
	// ---------------------------------------------------------------------
	// A_n < TE_6
	if( iGraphType == 0 && grBig->iGraphType == 4 && iVerticesBigCount == 7 )
	{
		vector<short unsigned int> iVTemp( grBig->iVertices.begin( ), grBig->iVertices.begin( ) + 2 );
		iVTemp.push_back( grBig->iVertices[6] );
		iVTemp.push_back( grBig->iVertices[4] );
		iVTemp.push_back( grBig->iVertices[5] );
		if( bAnSubAm( iVertices, iVTemp ) )
			return true;
		
		iVTemp = vector<short unsigned int>( grBig->iVertices.begin( ), grBig->iVertices.begin( ) + 2 );
		iVTemp.push_back( grBig->iVertices[6] );
		iVTemp.push_back( grBig->iVertices[2] );
		iVTemp.push_back( grBig->iVertices[3] );
		if( bAnSubAm( iVertices, iVTemp ) )
			return true;
		
		iVTemp.clear( );
		iVTemp.push_back( grBig->iVertices[5] );
		iVTemp.push_back( grBig->iVertices[4] );
		iVTemp.push_back( grBig->iVertices[6] );
		iVTemp.push_back( grBig->iVertices[2] );
		iVTemp.push_back( grBig->iVertices[3] );
		if( bAnSubAm( iVertices, iVTemp ) )
			return true;
	}
	
	// ---------------------------------------------------------------------
	// A_n < TE_m
	if( iGraphType == 0 && grBig->iGraphType == 4 )
	{
		vector<short unsigned int> iVTemp;
		
		//-------------------------------------
		// basis
		iVTemp = vector<short unsigned int>( grBig->iVertices.begin( ), grBig->iVertices.end( ) - 1 );
		
		if( bAnSubAm( iVertices, iVTemp ) )
			return true;
		
		iVTemp.clear( );
		
		// -------------------------------------
		// left & queue or right & queue
		if( iVerticesBigCount == 8 ) // TE7
		{
			iVTemp = vector<short unsigned int>( grBig->iVertices.begin( ), grBig->iVertices.begin( ) + 4 );
			iVTemp.push_back( grBig->iVertices[7] );
			if( bAnSubAm( iVertices, iVTemp ) )
				return true;
			
			iVTemp.clear( );
			iVTemp = vector<short unsigned int>( grBig->iVertices.begin( ) + 3, grBig->iVertices.end( ) - 1 );
			iVTemp.insert( iVTemp.begin( ), grBig->iVertices[7] );
			
			
			
			return bAnSubAm( iVertices, iVTemp );
		}
		else if( iVerticesBigCount == 9 ) // TE8
		{
			iVTemp = vector<short unsigned int>( grBig->iVertices.begin( ), grBig->iVertices.begin( ) + 3 );
			iVTemp.push_back( grBig->iVertices[8] );
			if( bAnSubAm( iVertices, iVTemp ) )
				return true;
			
			iVTemp.clear( );
			iVTemp = vector<short unsigned int>( grBig->iVertices.begin( ) + 2, grBig->iVertices.end( ) - 1 );
			iVTemp.insert( iVTemp.begin( ), grBig->iVertices[8] );
			return bAnSubAm( iVertices, iVTemp );
		}
	}
	
	// ---------------------------------------------------------------------
	// A_n < TF_4
	if( iGraphType == 0 && grBig->iGraphType == 5 )
	{
		if( iVerticesCount > 3 )
			return false;
		
		if( bAnSubAm( iVertices, vector<short unsigned int>( grBig->iVertices.begin( ), grBig->iVertices.begin( ) + 3 ) ) )
			return true;
		
		if( iVerticesCount == 2 )
			return bAnSubAm( iVertices, vector<short unsigned int>( grBig->iVertices.end( ) - 2, grBig->iVertices.end( ) ) );
		
		return false;
	}
	
	// ---------------------------------------------------------------------
	// A_n < TG_2
	if( iGraphType == 0 && grBig->iGraphType == 6 )
	{
		if( iVerticesCount > 2 )
			return false;
		
		return bAnSubAm( iVertices, vector<short unsigned int>( grBig->iVertices.begin( ), grBig->iVertices.begin( ) + 2 ) );
	}
	
	// ---------------------------------------------------------------------
	// B_n < TB_m
	if( iGraphType == 1 && grBig->iGraphType == 1 )
	{
		vector<short unsigned int> iVTemp( iVertices );
		reverse( iVTemp.begin( ), iVTemp.end( ) );
		
		vector<short unsigned int>::const_iterator itSearch( find( grBig->iVertices.begin( ), grBig->iVertices.end( ), iVertices[0] ) );
		
		if( itSearch == grBig->iVertices.end( ) )
			return false;
		
		// if we don't use any of the two ends (right ends)
		if( itSearch + 1 != grBig->iVertices.end( ) && itSearch + 2 != grBig->iVertices.end( ) )
			return ( iVTemp == vector<short unsigned int>( grBig->iVertices.begin( ), itSearch + 1 ) );
		else if( itSearch + 1 == grBig->iVertices.end( ) )
		{
			
			vector<short unsigned int> iVTemp2( vector<short unsigned int>( grBig->iVertices.begin( ), itSearch - 1 ) );
			iVTemp2.push_back( grBig->iVertices[ iVerticesBigCount - 1 ] );
			
			return ( iVTemp == iVTemp2 );
			
		}
		else // itSearch + 2 == grBig->iVertices.end( )
		{
			vector<short unsigned int> iVTemp2( vector<short unsigned int>( grBig->iVertices.begin( ), itSearch ) );
			iVTemp2.push_back( grBig->iVertices[ iVerticesBigCount - 2 ] );
			
			return ( iVTemp == iVTemp2 );
		}
		
		return false;
	}
	
	// ---------------------------------------------------------------------
	// B_n < TC_m
	if( iGraphType == 1 && grBig->iGraphType == 2 )
	{
		vector<short unsigned int>::const_iterator itSearch( find( grBig->iVertices.begin( ), grBig->iVertices.end( ), iVertices[0] ) );
		
		if( itSearch == grBig->iVertices.end( ) )
			return false;
		
		// [3, ..., 3, 4] < [3, ..., 3, 4]
		if( iVertices[ iVerticesCount - 1 ] == grBig->iVertices[ iVerticesBigCount - 1 ] )
			return ( iVertices == vector<short unsigned int>( itSearch, grBig->iVertices.end( ) ) );
		
		// [3, ..., 3, 4] < [4, 3, ..., 3]
		if( iVertices[ iVerticesCount - 1 ] == grBig->iVertices[0] )
		{
			vector<short unsigned int> iVTemp( grBig->iVertices.begin( ), itSearch + 1 );
			reverse( iVTemp.begin( ), iVTemp.end( ) );
			
			return ( iVTemp == iVertices );
		}
	}
	
	// ---------------------------------------------------------------------
	// B_n < TF_4
	if( iGraphType == 1 && grBig->iGraphType == 5 )
	{
		if( iVerticesCount > 4 )
			return false;
		
		// ----------------------------------------
		// ( [4,3] OR [4,3,3] ) < ( [3,3,4] )
		vector<short unsigned int> iVTemp( grBig->iVertices.begin( ) + ( iVerticesCount == 4 ? 0 : 1 ), grBig->iVertices.end( ) - 1 );
		if( iVertices == iVTemp )
			return true;
		
		// ----------------------------------------
		// [4,3] < ( [4,3] )
		iVTemp = vector<short unsigned int>( grBig->iVertices.begin( ) + 2, grBig->iVertices.end( ) );
		reverse( iVTemp.begin( ), iVTemp.end( ) );
		return ( iVertices == iVTemp );
	}
	
	// ---------------------------------------------------------------------
	// D_n < TB_m
	if( iGraphType == 3 && grBig->iGraphType == 1 )
	{
		if( iVerticesBigCount < 5 )
			return false;
		
		Graph g( vector<short unsigned int>( grBig->iVertices.begin( ) + 1, grBig->iVertices.end( ) ), 0, vector< bool >( false ), 3, true, 0 );
		
		return bIsSubgraphOf_spherical_spherical( &g );
	}
	
	// ---------------------------------------------------------------------
	// D_4 < TD_4
	if( iGraphType == 3 && iVerticesCount == 4 && grBig->iGraphType == 3 && iVerticesBigCount == 5 )
	{
		if( iVertices[1] != grBig->iVertices[4] ) // comparison of the central node
			return false;
		
		/*
		 * Remark:
		 * We should use find( grBig->iVertices.begin( ), grBig->iVertices.end( ) - 1, - ), instead of
		 * find( grBig->iVertices.begin( ), grBig->iVertices.end( ), - ) but it doesn't really matter
		 * since the elements of grBig->iVertices are all distincts.
		 */ 
		if( find( grBig->iVertices.begin( ), grBig->iVertices.end( ), iVertices[0] ) == grBig->iVertices.end( ) )
			return false;
		if( find( grBig->iVertices.begin( ), grBig->iVertices.end( ), iVertices[2] ) == grBig->iVertices.end( ) )
			return false;
		if( find( grBig->iVertices.begin( ), grBig->iVertices.end( ), iVertices[3] ) == grBig->iVertices.end( ) )
			return false;
		
		return true;
	}
	
	// ---------------------------------------------------------------------
	// D_n < TD_m
	if( iGraphType == 3 && grBig->iGraphType == 3 )
	{
		// Test 1
		Graph g( vector<short unsigned int>( grBig->iVertices.begin( ) + 1, grBig->iVertices.end( ) ), 0, vector< bool >( false ), 3, true, 0 );
		if( bIsSubgraphOf_spherical_spherical( &g ) )
			return true;
		
		// Test 2
		vector<short unsigned int> iVTemp( vector<short unsigned int>( grBig->iVertices.begin( ) + 2, grBig->iVertices.end( ) ) );
		iVTemp.insert( iVTemp.begin( ), grBig->iVertices[0] );
		g = Graph( iVTemp, 0, vector< bool >( false ), 3, true, 0 );
		if( bIsSubgraphOf_spherical_spherical( &g ) )
			return true;
		
		// Test 3
		iVTemp = vector<short unsigned int>( vector<short unsigned int>( grBig->iVertices.begin( ), grBig->iVertices.end( ) - 1 ) );
		iVTemp[0] = max( grBig->iVertices[0], grBig->iVertices[1] );
		iVTemp[1] = min( grBig->iVertices[0], grBig->iVertices[1] );
		reverse( iVTemp.begin( ), iVTemp.end( ) );
		g = Graph( iVTemp, 0, vector< bool >( false ), 3, true, 0 );
		if( bIsSubgraphOf_spherical_spherical( &g ) )
			return true;
		
		// Test 4
		iVTemp = vector<short unsigned int>( vector<short unsigned int>( grBig->iVertices.begin( ), grBig->iVertices.end( ) - 2 ) );
		iVTemp[0] = max( grBig->iVertices[0], grBig->iVertices[1] );
		iVTemp[1] = min( grBig->iVertices[0], grBig->iVertices[1] );
		iVTemp.push_back( grBig->iVertices[ iVerticesBigCount - 1 ] );
		reverse( iVTemp.begin( ), iVTemp.end( ) );
		g = Graph( iVTemp, 0, vector< bool >( false ), 3, true, 0 );
		if( bIsSubgraphOf_spherical_spherical( &g ) )
			return true;
	}
	
	// ---------------------------------------------------------------------
	// D_n < TE_6
	if( iGraphType == 3 && grBig->iGraphType == 4 && iVerticesBigCount == 7 )
	{
		// -------------------------------------------------
		// first test
		vector<short unsigned int> iVTemp( grBig->iVertices.begin( ), grBig->iVertices.begin( ) + 2 );
		iVTemp.push_back( grBig->iVertices[6] );
		iVTemp.push_back( min( grBig->iVertices[2], grBig->iVertices[4] ) );
		iVTemp.push_back( max( grBig->iVertices[2], grBig->iVertices[4] ) );
		Graph g( iVTemp, 0, vector< bool >( false ), 3, true, 0 );
		
		if( bIsSubgraphOf_spherical_spherical( &g ) )
			return true;
		
		// -------------------------------------------------
		// second test
		iVTemp.clear( );
		iVTemp.push_back( grBig->iVertices[5] );
		iVTemp.push_back( grBig->iVertices[4] );
		iVTemp.push_back( grBig->iVertices[6] );
		iVTemp.push_back( min( grBig->iVertices[2], grBig->iVertices[1] ) );
		iVTemp.push_back( max( grBig->iVertices[2], grBig->iVertices[1] ) );
		
		g = Graph( iVTemp, 0, vector< bool >( false ), 3, true, 0 );
		
		if( bIsSubgraphOf_spherical_spherical( &g ) )
			return true;
		
		// -------------------------------------------------
		// third test
		iVTemp.clear( );
		iVTemp.push_back( grBig->iVertices[3] );
		iVTemp.push_back( grBig->iVertices[2] );
		iVTemp.push_back( grBig->iVertices[6] );
		iVTemp.push_back( min( grBig->iVertices[4], grBig->iVertices[1] ) );
		iVTemp.push_back( max( grBig->iVertices[4], grBig->iVertices[1] ) );
		
		g = Graph( iVTemp, 0, vector< bool >( false ), 3, true, 0 );
		
		if( bIsSubgraphOf_spherical_spherical( &g ) )
			return true;
		
		return false;
	}
	
	// ---------------------------------------------------------------------
	// D_n < TE_m
	if( iGraphType == 3 && grBig->iGraphType == 4 )
	{
		unsigned int iBigQueueIndex( iVerticesBigCount == 8 ? 3 : 2 );
		
		// -------------------------------------------------
		// first test
		vector<short unsigned int> iVTemp( grBig->iVertices.begin( ), grBig->iVertices.begin( ) + iBigQueueIndex + 1 );
		iVTemp.push_back( min( grBig->iVertices[ iBigQueueIndex + 1 ], grBig->iVertices[ iVerticesBigCount - 1 ] ) );
		iVTemp.push_back( max( grBig->iVertices[ iBigQueueIndex + 1 ], grBig->iVertices[ iVerticesBigCount - 1 ] ) );
		Graph g( iVTemp, 0, vector< bool >( false ), 3, true, 0 );
		
		if( bIsSubgraphOf_spherical_spherical( &g ) )
			return true;
		
		// -------------------------------------------------
		// second test
		iVTemp = vector<short unsigned int>( grBig->iVertices.begin( ) + iBigQueueIndex, grBig->iVertices.end( ) - 1 );
		reverse( iVTemp.begin( ), iVTemp.end( ) );
		iVTemp.push_back( min( grBig->iVertices[ iBigQueueIndex - 1 ], grBig->iVertices[ iVerticesBigCount - 1 ] ) );
		iVTemp.push_back( max( grBig->iVertices[ iBigQueueIndex - 1 ], grBig->iVertices[ iVerticesBigCount - 1 ] ) );
		g = Graph( iVTemp, 0, vector< bool >( false ), 3, true, 0 );
		
		if( bIsSubgraphOf_spherical_spherical( &g ) )
			return true;
		
		return false;
	}
	
	// ---------------------------------------------------------------------
	// E_6 < TE_6
	if( iGraphType == 4 && iVerticesCount == 6 && grBig->iGraphType == 4 && iVerticesBigCount == 7 )
	{
		vector<short unsigned int> iVBigBasis;
		
		// We must find out what the base is
		if( iVertices[5] == grBig->iVertices[4] )
		{
			iVBigBasis.push_back( grBig->iVertices[0] );
			iVBigBasis.push_back( grBig->iVertices[1] );
			iVBigBasis.push_back( grBig->iVertices[6] );
			iVBigBasis.push_back( grBig->iVertices[2] );
			iVBigBasis.push_back( grBig->iVertices[3] );
		}
		else if( iVertices[5] == grBig->iVertices[1] )
		{
			iVBigBasis.push_back( grBig->iVertices[5] );
			iVBigBasis.push_back( grBig->iVertices[4] );
			iVBigBasis.push_back( grBig->iVertices[6] );
			iVBigBasis.push_back( grBig->iVertices[2] );
			iVBigBasis.push_back( grBig->iVertices[3] );
		}
		else if( iVertices[5] == grBig->iVertices[2] )
		{
			iVBigBasis.push_back( grBig->iVertices[0] );
			iVBigBasis.push_back( grBig->iVertices[1] );
			iVBigBasis.push_back( grBig->iVertices[6] );
			iVBigBasis.push_back( grBig->iVertices[4] );
			iVBigBasis.push_back( grBig->iVertices[5] );
		}
		else
			return false;
		
		if( iVBigBasis[0] > iVBigBasis[4] )
			reverse( iVBigBasis.begin( ), iVBigBasis.end( ) );
		
		return ( iVBigBasis == vector<short unsigned int>( iVertices.begin( ), iVertices.end( ) - 1 ) );
	}
	
	// ---------------------------------------------------------------------
	// E_n < TE_7, TE_8
	if( iGraphType == 4 && grBig->iGraphType == 4 )
	{	
		if( iVertices[ iVerticesCount - 1 ] != grBig->iVertices[ iVerticesBigCount - 1 ] )
			return false;
		
		// E_6, E_7 < TE7
		if( iVerticesCount <= 7 && iVerticesBigCount == 8 )
		{
			vector<short unsigned int> iVBigBasis( grBig->iVertices.begin( ) + 1, grBig->iVertices.begin( ) + iVerticesCount );
			if( bAnSubAm( vector<short unsigned int>( iVertices.begin( ), iVertices.end( ) - 1 ), iVBigBasis ) )
				return true;
			
			if( iVerticesCount == 7 && iVerticesBigCount == 8 ) // E_7 < TE_7
			{
				iVBigBasis = vector<short unsigned int>( grBig->iVertices.begin( ), grBig->iVertices.begin( ) + 6 );
				reverse( iVBigBasis.begin( ), iVBigBasis.end( ) );
				
				return (  vector<short unsigned int>( iVertices.begin( ), iVertices.end( ) - 1 ) == iVBigBasis );
			}
			
			return false;
		}
		else if( iVerticesBigCount == 9 ) // < TE_8
		{
			vector<short unsigned int> iVBigBasis( grBig->iVertices.begin( ), grBig->iVertices.begin( ) + iVerticesCount - 1 );
			
			if( vector<short unsigned int>( iVertices.begin( ), iVertices.end( ) - 1 ) == iVBigBasis )
				return true;
			
			if( iVerticesCount == 6 ) // E_5 is symmetric
			{
				reverse( iVBigBasis.begin( ), iVBigBasis.end( ) );
				
				return ( ( vector<short unsigned int>( iVertices.begin( ), iVertices.end( ) - 1 ) == iVBigBasis ) );
			}
			
			return false;
		}
		else
			return false;
	}
	
	// ---------------------------------------------------------------------
	// F_4 < TF_4
	if( iGraphType == 5 && grBig->iGraphType == 5 )
	{
		return bAnSubAm( iVertices, vector<short unsigned int>( grBig->iVertices.begin( ) + 1, grBig->iVertices.end( ) ) );
	}
	
	// ---------------------------------------------------------------------
	// G_4 < TB_m
	if( iGraphType == 6 && iDataSupp == 4 && grBig->iGraphType == 1 )
	{
		return ( ( iVertices[0] == grBig->iVertices[0] && iVertices[1] == grBig->iVertices[1] ) || ( iVertices[0] == grBig->iVertices[1] && iVertices[1] == grBig->iVertices[0] ) );
	}
	
	// ---------------------------------------------------------------------
	// G_4 < TC_m
	if( iGraphType == 6 && iDataSupp == 4 && grBig->iGraphType == 2 )
	{
		if( ( iVertices[0] == grBig->iVertices[0] && iVertices[1] == grBig->iVertices[1] ) || ( iVertices[0] == grBig->iVertices[1] && iVertices[1] == grBig->iVertices[0] ) )
			return true;
		
		return ( ( iVertices[0] == grBig->iVertices[ iVerticesBigCount - 1 ] && iVertices[1] == grBig->iVertices[ iVerticesBigCount - 2 ] ) || ( iVertices[0] == grBig->iVertices[ iVerticesBigCount - 2 ] && iVertices[1] == grBig->iVertices[ iVerticesBigCount - 1 ] ) );
	}
	
	// ---------------------------------------------------------------------
	// G_4 < TF_4
	if( iGraphType == 6 && iDataSupp == 4 && grBig->iGraphType == 5 )
	{
		return ( ( iVertices[0] == grBig->iVertices[2] && iVertices[1] == grBig->iVertices[3] ) || ( iVertices[0] == grBig->iVertices[3] && iVertices[1] == grBig->iVertices[2] ) );
	}
	
	// ---------------------------------------------------------------------
	// G_6 < TG_2
	if( iGraphType == 6 && iDataSupp == 6 && grBig->iGraphType == 6 )
	{
		return ( ( iVertices[0] == grBig->iVertices[1] && iVertices[1] == grBig->iVertices[2] ) || ( iVertices[1] == grBig->iVertices[1] && iVertices[0] == grBig->iVertices[2] ) );
	}
	
	return false;
}

bool Graph::bIsSubgraphOf_spherical_spherical( const Graph* grBig ) const
{
	size_t iVerticesCount( iVertices.size( ) ), iVerticesBigCount( grBig->iVertices.size( ) );
	
	if( iVerticesBigCount < iVerticesCount )
		return false;
	
	if( iVerticesCount == 1 )
		return ( find( grBig->iVertices.begin( ), grBig->iVertices.end( ), iVertices[0] ) != grBig->iVertices.end( ) );
	
	// ---------------------------------------------------------------------
	// A_n < A_m
	if( iGraphType == 0 && grBig->iGraphType == 0 )
	{
		return bAnSubAm( iVertices, grBig->iVertices );
	}
	
	// ---------------------------------------------------------------------
	// A_n < B_m
	if( iGraphType == 0 && grBig->iGraphType == 1 )
	{
		return bAnSubAm( iVertices, vector<short unsigned int>( grBig->iVertices.begin( ), grBig->iVertices.end( ) - 1 ) );
	}
	
	// ---------------------------------------------------------------------
	// A_n < D_m
	if( iGraphType == 0 && grBig->iGraphType == 3 )
	{
		// base plus une des deux extrémités
		vector<short unsigned int> iVTemp( grBig->iVertices.begin( ), grBig->iVertices.end( ) - 1 );
		if( bAnSubAm( iVertices, iVTemp ) )
			return true;
		
		// base plus l'autre extrémité
		iVTemp = vector<short unsigned int>( grBig->iVertices.begin( ), grBig->iVertices.end( ) - 2 );
		iVTemp.push_back( grBig->iVertices[ iVerticesBigCount - 1 ] );
		if( bAnSubAm( iVertices, iVTemp ) )
			return true;
		
		// A_3 dans le "bout du Y"
		if( iVerticesCount <= 3 && iVerticesBigCount > 3 )
		{
			iVTemp.clear( );
			iVTemp.push_back( min( grBig->iVertices[ iVerticesBigCount - 1 ], grBig->iVertices[ iVerticesBigCount - 2 ] ) );
			iVTemp.push_back(  grBig->iVertices[ iVerticesBigCount - 3 ] );
			iVTemp.push_back( max( grBig->iVertices[ iVerticesBigCount - 1 ], grBig->iVertices[ iVerticesBigCount - 2 ] ) );
		
			return bAnSubAm( iVertices, iVTemp );
		}
		
		return false;
	}
	
	// ---------------------------------------------------------------------
	// A_n < E_m
	if( iGraphType == 0 && grBig->iGraphType == 4 )
	{	
		// A_n dans la base du E_n
		vector<short unsigned int> iVTemp( grBig->iVertices.begin( ), grBig->iVertices.end( ) - 1 );
		if( bAnSubAm( iVertices, iVTemp ) )
			return true;
		
		// A_n contient la queue du E_m (début)
		iVTemp = vector<short unsigned int>( grBig->iVertices.begin( ), grBig->iVertices.begin( ) + 3 );
		iVTemp.push_back( grBig->iVertices[ iVerticesBigCount - 1 ] );
		if( bAnSubAm( iVertices, iVTemp ) )
			return true;
		
		// A_n contient la queue du E_m (fin)
		iVTemp = vector<short unsigned int>( grBig->iVertices.begin( ) + 2, grBig->iVertices.end( ) - 1 );
		iVTemp.insert( iVTemp.begin( ), grBig->iVertices[ iVerticesBigCount - 1 ] );
		
		return bAnSubAm( iVertices, iVTemp );
	}
	
	// ---------------------------------------------------------------------
	// A_2 < F_4
	if( iGraphType == 0 && iVerticesCount == 2 && grBig->iGraphType == 5 )
	{
		return ( ( iVertices[0] == grBig->iVertices[0] && iVertices[1] == grBig->iVertices[1] )
			|| ( iVertices[1] == grBig->iVertices[0] && iVertices[0] == grBig->iVertices[1] )
			|| ( iVertices[0] == grBig->iVertices[2] && iVertices[1] == grBig->iVertices[3] )
			|| ( iVertices[1] == grBig->iVertices[2] && iVertices[0] == grBig->iVertices[3] )
		);
	}
	
	// ---------------------------------------------------------------------
	// A_n < H_m
	if( iGraphType == 0 && grBig->iGraphType == 7 )
	{
		return bAnSubAm( iVertices, vector<short unsigned int>( grBig->iVertices.begin( ), grBig->iVertices.end( ) - 1 ) );
	}
	
	// ---------------------------------------------------------------------
	// B_n < B_m
	if( iGraphType == 1 && grBig->iGraphType == 1 )
	{
		return bAnSubAm( iVertices, grBig->iVertices );
	}
	
	// ---------------------------------------------------------------------
	// B_3 < F_4
	if( iGraphType == 1 && grBig->iGraphType == 5 && iVerticesCount == 3 )
	{
		if( iVertices == vector<short unsigned int>( grBig->iVertices.begin( ), grBig->iVertices.begin( ) + 3 ) )
			return true;
		
		vector<short unsigned int> iVTemp( grBig->iVertices.begin( ) + 1, grBig->iVertices.begin( ) + 4 );
		reverse( iVTemp.begin( ), iVTemp.end( ) );
		
		return ( iVertices == iVTemp );	
	}
	
	// ---------------------------------------------------------------------
	// D_n < D_m
	if( iGraphType == 3 && grBig->iGraphType == 3 )
	{
		if( iVerticesCount == 4 ) // D_4 est très symétrique
		{	
			vector<short unsigned int> iVTemp1;
			iVTemp1.push_back( iVertices[0] );
			iVTemp1.push_back( iVertices[2] );
			iVTemp1.push_back( iVertices[3] );
			
			vector<short unsigned int> iVTemp2;
			iVTemp2.push_back( grBig->iVertices[ iVerticesBigCount - 1 ] );
			iVTemp2.push_back( grBig->iVertices[ iVerticesBigCount - 2 ] );
			iVTemp2.push_back( grBig->iVertices[ iVerticesBigCount - 4 ] );
			
			sort( iVTemp1.begin( ), iVTemp1.end( ) );
			sort( iVTemp2.begin( ), iVTemp2.end( ) );
			
			return ( iVTemp1 == iVTemp2 );
		}
		
		// Autres D_n
		if( iVertices[ iVerticesCount - 2 ] != grBig->iVertices[ iVerticesBigCount - 2 ] || iVertices[ iVerticesCount - 1 ] != grBig->iVertices[ iVerticesBigCount - 1 ] )
			return false;
		
		vector<short unsigned int> iVTemp1( iVertices.begin( ), iVertices.begin( ) + iVerticesCount - 2 );
		vector<short unsigned int> iVTemp2( grBig->iVertices.begin( ) + iVerticesBigCount - iVerticesCount, grBig->iVertices.begin( ) + iVerticesBigCount - 2  );
		
		return ( iVTemp1 == iVTemp2 );
	}
	
	// ---------------------------------------------------------------------
	// D_n < E_m
	if( iGraphType == 3 && grBig->iGraphType == 4 )
	{
		// -------------------------------------------------
		// first test
		vector<short unsigned int> iVTemp( grBig->iVertices.begin( ), grBig->iVertices.begin( ) + 3 );
		iVTemp.push_back( min( grBig->iVertices[3], grBig->iVertices[ iVerticesBigCount - 1 ] ) );
		iVTemp.push_back( max( grBig->iVertices[3], grBig->iVertices[ iVerticesBigCount - 1 ] ) );
		Graph g( iVTemp, 0, vector< bool >( false ), 3, true, 0 );
		
		if( bIsSubgraphOf_spherical_spherical( &g ) )
			return true;
		
		// -------------------------------------------------
		// second test
		iVTemp = vector<short unsigned int>( grBig->iVertices.begin( ) + 2, grBig->iVertices.end( ) - 1 );
		reverse( iVTemp.begin( ), iVTemp.end( ) );
		iVTemp.push_back( min( grBig->iVertices[1], grBig->iVertices[ iVerticesBigCount - 1 ] ) );
		iVTemp.push_back( max( grBig->iVertices[1], grBig->iVertices[ iVerticesBigCount - 1 ] ) );
		g = Graph( iVTemp, 0, vector< bool >( false ), 3, true, 0 );
		
		if( bIsSubgraphOf_spherical_spherical( &g ) )
			return true;
		
		return false;
	}
	
	// ---------------------------------------------------------------------
	// E_n < E_m
	if( iGraphType == 4 && grBig->iGraphType == 4 )
	{
		if( iVertices[ iVerticesCount - 1 ] != grBig->iVertices[ iVerticesBigCount - 1 ] )
			return false;
		
		vector<short unsigned int> iVTemp1( iVertices.begin( ), iVertices.end( ) - 1 );
		vector<short unsigned int> iVTemp2( grBig->iVertices.begin( ), grBig->iVertices.begin( ) + iVerticesCount - 1 );
		
		if( iVTemp1 == iVTemp2 )
			return true;
		
		if( iVerticesCount == 6 ) // E_5 a une base symétrique
		{
			reverse( iVTemp1.begin( ), iVTemp1.end( ) );
			if( iVTemp1 == iVTemp2 )
				return true;
		}
		
		return false;
	}
	
	// ---------------------------------------------------------------------
	// G_4 < B_n
	if( iGraphType == 6 && iDataSupp == 4 && grBig->iGraphType == 1 )
	{
		return ( ( iVertices[0] == grBig->iVertices[iVerticesBigCount-2] && iVertices[1] == grBig->iVertices[iVerticesBigCount-1] ) || ( iVertices[0] == grBig->iVertices[iVerticesBigCount-1] && iVertices[1] == grBig->iVertices[iVerticesBigCount-2] ) );
	}
	
	// ---------------------------------------------------------------------
	// G_m < G_m
	if( iGraphType == 6 && grBig->iGraphType == 6 && iDataSupp == grBig->iDataSupp && iVertices == grBig->iVertices )
	{
		return ( ( iVertices[0] == grBig->iVertices[iVerticesBigCount-2] && iVertices[1] == grBig->iVertices[iVerticesBigCount-1] ) || ( iVertices[0] == grBig->iVertices[iVerticesBigCount-1] && iVertices[1] == grBig->iVertices[iVerticesBigCount-2] ) );
	}
	
	// ---------------------------------------------------------------------
	// G_5 < H_m
	if( ( iGraphType == 6 && iDataSupp == 5 ) && grBig->iGraphType == 7 )
	{
		return ( ( iVertices[0] == grBig->iVertices[ iVerticesBigCount - 2 ] && iVertices[1] == grBig->iVertices[ iVerticesBigCount - 1 ] ) || ( iVertices[1] == grBig->iVertices[ iVerticesBigCount - 2 ] && iVertices[0] == grBig->iVertices[ iVerticesBigCount - 1 ] ) );
	}
	
	// ---------------------------------------------------------------------
	// H_3 < H_4
	if( iGraphType == 7 && grBig->iGraphType == 7 )
	{
		return ( vector<short unsigned int>( grBig->iVertices.begin( ) + 1, grBig->iVertices.begin( ) + iVerticesCount + 1 ) == iVertices );
	}
	
	return false;
}

bool Graph::bAnSubAm( const vector<short unsigned int>& iSubV, const vector<short unsigned int>& iBigV )
{
	vector<short unsigned int>::const_iterator itBig, itSub, it( find( iBigV.begin( ), iBigV.end( ), iSubV[0] ) );
	
	if( it == iBigV.end( ) )
		return false;
	
	if( iSubV.size( ) == 1 )
		return true;
	
	if( iSubV.size( ) > iBigV.size( ) )
		return false;
	
	if( ( it + 1 ) == iBigV.end( ) && *( it - 1 ) != iSubV[1] )
		return false;
	
	// ----------------------------------------------------------------------
	// recherche en avant depuis it
	if( ( it + 1 ) != iBigV.end( ) && *( it + 1 ) == iSubV[1] )
	{
		itBig = it;
		for( itSub = iSubV.begin( ); itSub != iSubV.end( ) && itBig != iBigV.end( ); ++itSub )
		{
			if( *itSub != *itBig )
				break;
			++itBig;
		}
		
		return ( itSub == iSubV.end( ) );
	}
	
	// ----------------------------------------------------------------------
	// recherche en arrière depuis it
	itBig = it;
	for( itSub = iSubV.begin( ); itSub != iSubV.end( ) && itBig >= iBigV.begin( ); ++itSub )
	{
		if( *itSub != *itBig )
			return false;
		
		--itBig;
	}
	
	return (  itSub == iSubV.end( ) );
}

bool operator==( const Graph &g1, const Graph &g2 )
{
	return ( g1.iGraphType == g2.iGraphType && g1.iDataSupp == g2.iDataSupp && g1.iVertices == g2.iVertices && g1.bSpherical == g2.bSpherical );
}

bool operator<( const Graph &g1, const Graph &g2 )
{
	if( g1 == g2 )
		return false;
	
	if( g1.bSpherical && !g2.bSpherical )
		return true;
	if( !g1.bSpherical && g2.bSpherical )
		return false;
	
	if( g1.iVertices.size() < g2.iVertices.size() )
		return true;
	if( g1.iVertices.size() > g2.iVertices.size() )
		return false;
	
	if( g1.iGraphType < g2.iGraphType )
		return true;
	if( g1.iGraphType > g2.iGraphType )
		return false;
	
	if( g1.iGraphType == 6 && g2.iGraphType == 6 )
	{
		if( g1.iDataSupp < g2.iDataSupp )
			return true;
		if( g1.iDataSupp > g2.iDataSupp )
			return false;
	}
	
	if( g1.iVertices < g2.iVertices )
		return true;
	if( g1.iVertices > g2.iVertices )
		return false;
	
	throw( string( "Graph::operator<: One missed case" ) );
	
	return false;
}
