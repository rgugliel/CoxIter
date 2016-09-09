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

#include "gbd.h"

GBD::GBD( CoxIter* ci )
: ci( ci ), iVerticesCount( ci->get_iVerticesCount( ) ), iNewVerticesCount( 0 ), iCox( ci->get_iCoxeterMatrix( ) )
{
}

bool GBD::bIsVertexAdmissible( const string& strVertexName )
{
	
	if( !ci->bIsVertexValid( strVertexName ) )
	{
		strError = "Invalid vertex name: " + strVertexName;
		return false;
	}
	
	unsigned int iV( ci->get_iVertexIndex( strVertexName ) );
	
	for( unsigned int i( 0 ); i < iVerticesCount; i++ )
	{
		if( iCox[ iV ][ i ] != 0 && iCox[ iV ][ i ] != 1 && ( iCox[ iV ][ i ] % 2 ) )
		{
			strError = "m(" + strVertexName + "," + ci->get_strVertexLabel( i ) + ") is not even";
			return false;
		}
	}
	
	return true;
}

bool GBD::removeVertex( const string& strVertexName )
{
	unsigned int iWeight;
	
	// -------------------------------------------------------
	// some verifications, firsts constructions
	if( !bIsVertexAdmissible( strVertexName ) )
		return false;
	
	iVertex = ci->get_iVertexIndex( strVertexName );
	
	// Remark: we don't merge the previous and the next loop. Thus, if there is a problem, we do not change the CoxIter object.
	
	// -------------------------------------------------------
	// new vertices
	for( unsigned int i( 0 ); i < iVerticesCount; i++ )
	{
		// If they don't commute, this will add a new vertex
		if( iCox[ iVertex ][ i ] != 2 )
		{
			NewVertex nv;
			nv.iIndex = iVerticesCount + iNewVerticesCount - 1;

			nv.iOriginVertex = i;
			nv.strLabel = ci->get_strVertexLabel( iVertex ) + "_" + ci->get_strVertexLabel( i );
			
			newVertices.push_back( nv );
			ci->map_vertices_labels_addReference( nv.strLabel );
			
			iNewVerticesCount++;
		}
	}

	ci->map_vertices_labels_removeReference( iVertex );
	
	// -------------------------------------------------------
	// new adjacency matrix
	for( unsigned int i( 0 ); i <= ( iVerticesCount + iNewVerticesCount - 1 ); i++ ) // the -1 is for the vertex we removed
	{
		if( i == iVertex )
			continue;
		
		iNewCox.push_back( vector< unsigned int >( iVerticesCount + iNewVerticesCount - 1, 2 ) ); // TODO: tout construire en une fois
	}
	
	// we fill the old values (upper left square)
	for( unsigned int i( 0 ); i < iVerticesCount; i++ )
	{
		if( i == iVertex )
			continue;
		
		for( unsigned int j( 0 ); j <= i; j++ )
			iNewCox[ i > iVertex ? i-1 : i ][ j > iVertex ? j-1 : j ] = iNewCox[ j > iVertex ? j-1 : j ][ i > iVertex ? i-1 : i ] = iCox[i][j];
	}
	
	// -------------------------------------------------------
	// weights of the new edges
	for( vector< NewVertex >::const_iterator itNew( newVertices.begin( ) ); itNew != newVertices.end( ); ++itNew ) // Going through the set of new vertices
	{
		// A new vertex is of the form: t0 * sj * t0, where sj = NewVertex.iOriginVertex
		
		// for every old vertex
		for( unsigned int i( 0 ); i < iVerticesCount; i++ ) // Goal: find m( t0 * sj * t0, si )
		{
			if( i == iVertex ) // we removed this vertex
				continue;
			
			if( (*itNew).iOriginVertex == i ) // m( t0 * si * t0, si ) = m( t0, si ) / 2
				iWeight = iCox[ iVertex ][ i ] == 0 || iCox[ iVertex ][ i ] == 1 ? iCox[ iVertex ][ i ] : iCox[ iVertex ][ i ] / 2;
			else
			{
				if( iCox[i][ iVertex ] == 2 ) // If si commutes with t0, then m( t0 * sj * t0, si ) = m( sj, si )
					iWeight = iCox[ (*itNew).iOriginVertex ][ i ];
				else if( iCox[ (*itNew).iOriginVertex ][iVertex] != 2 )
				{
					if( iCox[iVertex][i] == 4 && iCox[iVertex][ (*itNew).iOriginVertex ] == 4 && iCox[i][ (*itNew).iOriginVertex ] == 2 )
						iWeight = 0;
					else
						iWeight = 1;
				}
				else
					throw( string( "GBD::removeVertex: Error" ) ); // TODO: idem
			}
			
			iNewCox[ i > iVertex ? i - 1 : i ][ (*itNew).iIndex ] = iWeight;
			iNewCox[ (*itNew).iIndex ][ i > iVertex ? i - 1 : i ] = iWeight;
		}

		// for every new vertex t0 * si * t0
		for( vector< NewVertex >::const_iterator itNewSub( newVertices.begin( ) ); itNewSub != itNew; ++itNewSub )
		{
			iNewCox[ (*itNew).iIndex ][ (*itNewSub).iIndex ] = iCox[ (*itNew).iOriginVertex ][ (*itNewSub).iOriginVertex ];
			iNewCox[ (*itNewSub).iIndex ][ (*itNew).iIndex ] = iCox[ (*itNew).iOriginVertex ][ (*itNewSub).iOriginVertex ];
		}
	}
	
	ci->set_iCoxeterMatrix( iNewCox );
	
	return true;
}

void GBD::printMatrix( vector< vector< unsigned int > >* iMatrix )
{
	for( vector< vector< unsigned int > >::const_iterator itRow( iMatrix->begin( ) ); itRow != iMatrix->end( ); ++itRow )
	{
		for( vector< unsigned int >::const_iterator itCol( (*itRow).begin( ) ); itCol != (*itRow).end( ); ++itCol )
			cout << *itCol << " ";
		
		cout << endl;
	}
}

string GBD::get_strError( ) const
{
	return strError;
}