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

#include "graphs.list.n.h"

GraphsListN::GraphsListN( unsigned int iVerticesCount, vector< string > *ptr_map_vertices_indexToLabel )
: iVerticesCount( iVerticesCount ), ptr_map_vertices_indexToLabel( ptr_map_vertices_indexToLabel )
{
}

size_t GraphsListN::size( ) const
{
	return graphs.size( );
}

Graph* GraphsListN::next( const size_t &iGraphIndex )
{
	return &graphs[ iGraphIndex ];
}

Graph* GraphsListN::begin( )
{
	if( graphs.empty( ) )
		return 0;
	
	return &graphs[0];
}

void GraphsListN::addGraph( vector< short unsigned int > iVertices, const vector< bool > &bVerticesLinkable, const unsigned int &iType, bool bSpherical, const short unsigned int &iVertexSupp1, const short unsigned int &iVertexSupp2, const unsigned int &iDataSupp )
{
	unsigned int short iTemp, i;
	
	if( bSpherical )
	{
		if( iType == 0 ) // An
		{
			// on veut une fois les A_n (et pas une fois pour chaque sens de lecture)
			if( iVertices.front( ) > iVertices.back( ) )
				reverse( iVertices.begin( ), iVertices.end( ) );
		}
		else if( iType == 1 ) // Bn
		{
			iVertices.push_back( iVertexSupp1 );
		}
		else if( iType == 3 ) // Dn
		{
			if( iVerticesCount > 4 )
			{
				iTemp = iVertices[ iVerticesCount - 2 ];
				iVertices[ iVerticesCount - 2 ] = min( iTemp, iVertexSupp1 );
				iVertices.push_back( max( iTemp, iVertexSupp1 ) );
			}
			else
			{
				if( iVertices[0] > iVertices[2] )
				{
					iTemp = iVertices[0];
					iVertices[0] = iVertices[2];
					iVertices[2] = iTemp;
				}
				
				if( iVertexSupp1 < iVertices[0] )
				{
					iVertices.push_back( iVertices[2] );
					iVertices[2] = iVertices[0];
					iVertices[0] = iVertexSupp1;
				}
				else if( iVertexSupp1 < iVertices[2] )
				{
					iVertices.push_back( iVertices[2] );
					iVertices[2] = iVertexSupp1;
				}
				else
					iVertices.push_back( iVertexSupp1 );
				
			}
		}
		else if( iType == 4 ) // En
		{
			if( iVerticesCount == 6 )
			{
				if( iVertices.front( ) > iVertices.back( ) )
					reverse( iVertices.begin( ), iVertices.end( ) );
			}
			iVertices.push_back( iVertexSupp1 );
		}
		else if( iType == 5 ) // Fn
		{
			iVertices.push_back( iVertexSupp1 );
			if( iVertices.front( ) > iVertices.back( ) )
				reverse( iVertices.begin( ), iVertices.end( ) );
		}
		else if( iType == 7 ) // Hn
		{
			iVertices.push_back( iVertexSupp1 );
		}
	}
	else
	{
		if( iType == 0 && iDataSupp )
		{
			int iMinIndex(0);
			unsigned int iMinValue( iVertices[0] );
			
			iVertices.push_back( iVertexSupp1 );
			vector< short unsigned int > iVerticesTemp( iVertices );
			
			// on répère le sommet avec l'indice le plus petit
			for( i = 1; i < iVerticesCount; i++ )
			{
				if( iVertices[i] < iMinValue )
				{
					iMinValue = iVertices[i];
					iMinIndex = i;
				}
			}
				
			// et le sens de parcours
			int iDirection( iVertices[ !iMinIndex ? (iVerticesCount - 1) : (iMinIndex - 1) ] > iVertices[ iMinIndex == (int)( iVerticesCount - 1 ) ? 0 : ( iMinIndex + 1 ) ] ? 1 : -1 );
			
			// on réordonne
			for( unsigned int j(0); j < iVerticesCount; j++ )
			{
				iVertices[j] = iVerticesTemp[iMinIndex];
				
				iMinIndex += iDirection;
				if( iMinIndex == -1 )
					iMinIndex = iVerticesCount - 1;
				else if( iMinIndex == (int)iVerticesCount )
					iMinIndex = 0;
			}
		}
		else if( iType == 1 ) // TBn
		{
			if( iVerticesCount == 4 ) // TB3
			{
				i = iVertices[1];
				iTemp = max( iVertices[0], iVertices[2] );
				iVertices[1] = min( iVertices[0], iVertices[2] ); 
				iVertices[2] = iTemp;
				iVertices[0] = i;
				
				iVertices.insert( iVertices.begin( ), iVertexSupp1 );
			}
			else // autres \tilde Bn
			{	
				iTemp = min( iVertices[ iVerticesCount - 3 ], iVertexSupp1 );
				iVertices.push_back( max( iVertices[ iVerticesCount - 3 ], iVertexSupp1 ) );
				iVertices[ iVerticesCount - 3 ] = iTemp;
				
				iVertices.insert( iVertices.begin( ), iVertexSupp2 );
			}
		}
		else if( iType == 2 ) // \tilde Cn
		{
			iVertices.insert( iVertices.begin( ), iVertexSupp1 );
			iVertices.push_back( iVertexSupp2 );
				
			if( iVertices[0] > iVertices[ iVerticesCount - 1 ] )
				reverse( iVertices.begin( ), iVertices.end( ) );
		}
		else if( iType == 3 ) // \tilde Dn
		{
			if( iVerticesCount >= 6 )
			{
				// le vecteur iVerticesBase contient la base (i.e. sans les 4 extrémités)
				#ifdef WIN32 // TODO: utile?
					vector< short unsigned int > iVerticesBase;
					iTemp = iVerticesCount - 3;
					for( i = 1; i < iTemp; i++ )
						iVerticesBase.push_back( iVertices[i] );
				#else
					vector< short unsigned int > iVerticesBase( iVertices.begin( ) + 1, iVertices.begin( ) + iVerticesCount - 3 ); // ne marche pas sous Visual Studio 2010 & 2012 malheureusement
				#endif

				// on regarde si on doit modifier l'ordre
				if( iVerticesBase[0] > iVerticesBase[ iVerticesCount - 5 ] )
				{
					reverse( iVerticesBase.begin( ), iVerticesBase.end( ) );
					
					// ajout des 4 extrémités
					iVerticesBase.push_back( min( iVertices[0], iVertexSupp1 ) );
					iVerticesBase.push_back( max( iVertices[0], iVertexSupp1 ) );
					iVerticesBase.insert( iVerticesBase.begin( ), max( iVertices[ iVerticesCount - 3 ], iVertexSupp2 ) );
					iVerticesBase.insert( iVerticesBase.begin( ), min( iVertices[ iVerticesCount - 3 ], iVertexSupp2 ) );
				}
				else
				{
					// ajout des 4 extrémités
					iVerticesBase.push_back( min( iVertices[ iVerticesCount - 3 ], iVertexSupp2 ) );
					iVerticesBase.push_back( max( iVertices[ iVerticesCount - 3 ], iVertexSupp2 ) );
					iVerticesBase.insert( iVerticesBase.begin( ), max( iVertices[0], iVertexSupp1 ) );
					iVerticesBase.insert( iVerticesBase.begin( ), min( iVertices[0], iVertexSupp1 ) );
				}

				iVertices = iVerticesBase;
			}
			else // le TD4, "+" est traité différemment
			{
				vector<short unsigned int> iVerticesTemp( 4, 0 );
				iVerticesTemp[0] = iVertices[0];
				iVerticesTemp[1] = iVertices[2];
				iVerticesTemp[2] = iVertexSupp1;
				iVerticesTemp[3] = iVertexSupp2;
				sort( iVerticesTemp.begin( ), iVerticesTemp.end( ) );
				iVerticesTemp.push_back( iVertices[1] );
				iVertices = iVerticesTemp;
			}
		}
		else if( iType == 4 ) // \tilde En
		{
			if( iVerticesCount == 7 ) // TE6
			{
				// TODO: refaire l'encodage de ce graphe et modifer la fonction bIsSubgraphOf_spherical_euclidean?
				vector< short unsigned int > iVerticesTemp;
				unsigned int iMin( min( min( iVertices[1], iVertices[3] ), iVertexSupp1 ) );
				
				if( iMin == iVertices[1] )
				{
					iVerticesTemp.push_back( iVertices[0] );
					iVerticesTemp.push_back( iVertices[1] );
					
					if( min( iVertices[3], iVertexSupp1 ) == iVertices[3] )
					{
						iVerticesTemp.push_back( iVertices[3] );
						iVerticesTemp.push_back( iVertices[4] );
						iVerticesTemp.push_back( iVertexSupp1 );
						iVerticesTemp.push_back( iVertexSupp2 );
					}
					else
					{
						iVerticesTemp.push_back( iVertexSupp1 );
						iVerticesTemp.push_back( iVertexSupp2 );
						iVerticesTemp.push_back( iVertices[3] );
						iVerticesTemp.push_back( iVertices[4] );
					}
				}
				else if( iMin == iVertices[3] )
				{
					iVerticesTemp.push_back( iVertices[4] );
					iVerticesTemp.push_back( iVertices[3] );
					
					if( min( iVertices[1], iVertexSupp1 ) == iVertices[1] )
					{
						iVerticesTemp.push_back( iVertices[1] );
						iVerticesTemp.push_back( iVertices[0] );
						iVerticesTemp.push_back( iVertexSupp1 );
						iVerticesTemp.push_back( iVertexSupp2 );
					}
					else
					{
						iVerticesTemp.push_back( iVertexSupp1 );
						iVerticesTemp.push_back( iVertexSupp2 );
						iVerticesTemp.push_back( iVertices[1] );
						iVerticesTemp.push_back( iVertices[0] );
					}
				}
				else
				{
					iVerticesTemp.push_back( iVertexSupp2 );
					iVerticesTemp.push_back( iVertexSupp1 );
					
					if( min( iVertices[1], iVertices[3] ) == iVertices[1] )
					{
						iVerticesTemp.push_back( iVertices[1] );
						iVerticesTemp.push_back( iVertices[0] );
						iVerticesTemp.push_back( iVertices[3] );
						iVerticesTemp.push_back( iVertices[4] );
					}
					else
					{
						iVerticesTemp.push_back( iVertices[3] );
						iVerticesTemp.push_back( iVertices[4] );
						iVerticesTemp.push_back( iVertices[1] );
						iVerticesTemp.push_back( iVertices[0] );
					}	
				}
				
				iVerticesTemp.push_back( iVertices[2] );
				iVertices = iVerticesTemp;
			}
			else if( iVerticesCount == 8  ) // TE7
			{
				if( iVertices.front( ) > iVertices.back( ) ) // la base du \tilde E7 est symétrique
					reverse( iVertices.begin( ), iVertices.end( ) );
				
				iVertices.push_back( iVertexSupp1 );
			}
			else // TE8
			{
				iVertices.push_back( iVertexSupp1 );
			}
		}
		else if( iType == 5 )
		{
			iVertices.push_back( iVertexSupp1 );
		}
		else if( iType == 6 && iDataSupp ) // \tilde G_2
		{
			iVertices.push_back( iVertexSupp1 );
		}
	}
	
	Graph g( iVertices, ptr_map_vertices_indexToLabel, bVerticesLinkable, iType, bSpherical, iDataSupp );
	
	/*
	unsigned int iMax( graphs.size( ) );
	for( unsigned int i(0); i < iMax; i++ )
	{
		if( g == graphs[i] )
			return;
	} TODO*/

	auto it( lower_bound( graphs.begin(), graphs.end(), g ) );
	if( it == graphs.end() || !(*it == g) )
		graphs.insert( it, g );
}

bool GraphsListN::addGraphsList( const GraphsListN& gln )
{
	if( iVerticesCount != gln.get_iVerticesCount( ) )
		return false;
	
	vector< Graph > gr( gln.get_graphs( ) );
	graphs.insert( graphs.end( ), gr.begin( ), gr.end( ) );
	
	return true;
}

unsigned int GraphsListN::get_iVerticesCount( ) const
{
	return iVerticesCount;
}

vector< Graph > GraphsListN::get_graphs( ) const
{
	return graphs;
}

ostream& operator<<( ostream &o, const GraphsListN &g )
{
	o << "\tGraphs of rank " << g.iVerticesCount << endl;
	unsigned int iMax( g.graphs.size( ) );
	
	for( unsigned int i(0); i < iMax; ++i )
		o << g.graphs[i];
	
	return o;
}
