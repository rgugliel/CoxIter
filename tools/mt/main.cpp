#include <iostream>
#include <chrono>

#include "../../gbd.h"
#include "bruteforce.h"
#include "multitools.h"
#include "../../coxiter.h"

using namespace std;

int main_gen( )
{
	MultiTools m;
	m.createCycles_main( );
	return 0;
}

int main_vinb77( )
{
	// Vinberg, 2014, dimension 18, partie 1
	CoxIter ci( "../../../graphs/temp/18-vinb77_base.in", false, "", false, true, false, false );
	if( !ci.readGraph( ) )
	{
		cout << "Erreur lecture graphe" << endl;
		return false;
	}
	vector< vector< unsigned int > > iCox( ci.get_iCoxeterMatrix( ) );
	
	/*
	vector< unsigned int > iVerticesStar( { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20, 21, 22, 23, 25, 26, 27 } );
	vector< unsigned int > iFirstType( { 19, 28, 29, 30, 31 } );
	vector< unsigned int > iSecondType( { 24, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50 } );*/
	
	
	vector< unsigned int > iVerticesStar( { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 19, 20, 21, 22, 24 } );
	vector< unsigned int > iFirstType( { 18, 25, 26 } );
	vector< unsigned int > iSecondType( { 23, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37 } );
	
	for( vector< unsigned int >::const_iterator itF( iFirstType.begin( ) ); itF != iFirstType.end( ); itF++ )
	{
		for( vector< unsigned int >::const_iterator itS( iSecondType.begin( ) ); itS != iSecondType.end( ); itS++ )
		{
			bool bBrokenLine( true );
			
			// premier test
			for( vector< unsigned int >::const_iterator itStar( iVerticesStar.begin( ) ); itStar != iVerticesStar.end( ); ++itStar )
			{
				if( iCox[ *itStar - 1 ][ *itF - 1 ] != 2 && iCox[ *itStar - 1 ][ *itS - 1 ] != 2 )
				{
					bBrokenLine = false;
					break;
				}
	
				/*
				for( vector< unsigned int >::const_iterator it1( iFirstType.begin( ) ); it1 != iFirstType.end( ); it1++ )
				{
					if( iCox[ *itStar - 1 ][ *it1 - 1 ] != 2 )
					{
						for( vector< unsigned int >::const_iterator it2( iSecondType.begin( ) ); it2 != iSecondType.end( ); it2++ )
						{
							if( iCox[ *itStar - 1 ][ *it2 - 1 ] != 2 )
							{
								bBrokenLine = false;
								break;
							}
						}
					}
					
					if( !bBrokenLine )
						break;
				}*/
				
				if( !bBrokenLine )
					break;
			}
			
			// second test
			if( bBrokenLine )
			{
				for( vector< unsigned int >::const_iterator it1( iVerticesStar.begin( ) ); it1 != iVerticesStar.end( ); ++it1 )
				{
					for( vector< unsigned int >::const_iterator it2( iVerticesStar.begin( ) ); it2 != iVerticesStar.end( ); ++it2 )
					{
						if( *it1 == *it2 || iCox[ *it1 - 1 ][ *it2 - 1 ] == 2 )
							continue;
						
						if( iCox[ *it1 - 1 ][ *itF - 1 ] != 2 && iCox[ *it2 - 1 ][ *itS - 1 ] != 2 )
						{
							bBrokenLine = false;
							break;
						}
					}
					if( !bBrokenLine )
						break;
				}
			}

			iCox[ *itF - 1 ][ *itS - 1 ] = iCox[ *itS - 1 ][ *itF - 1 ] = ( bBrokenLine ? 1 : 0 );
		}
	}
	
	ci.set_iCoxeterMatrix( iCox );
	ci.writeGraph( "../../../graphs/temp/18-vinb77_final" );
	
	cout << "Fin" << endl;
	return 0;
}

int main_( )
{
	chrono::time_point< std::chrono::system_clock > timeStart, timeEnd;
	
	CoxIter ci( "../../../graphs/17-vinb85.in", false, "", false, true, false, false );
	
	if( !ci.readGraph( ) )
		return false;
	
	ci.exploreGraph( );
	ci.computeGraphsProducts( );
	
	timeStart = chrono::system_clock::now();
	for( unsigned int i( 0 ); i < 1000000000; i++ )
	{
		ci.isFiniteCovolume( );
		ci.iIsGraphCocompact( );
	}
	timeEnd = chrono::system_clock::now();
	cout << "\tComputation time: " << chrono::duration <double, milli>(timeEnd-timeStart).count( ) / 1000 << "\n" << endl;
	
	return 0;
}

int main__( ) 
{
	unsigned int iDimension( 3 );
	try{
		BruteForce bf( iDimension );
		
		bf.basisMatrix_create( iDimension + 1 );
		
		for( unsigned int i( 0 ); i <= iDimension; i++ )
		{
			for( unsigned int j( 0 ); j < i; j++ )
				bf.variableEdge_add( i, j, vector< unsigned int >{ 2, 3, 4, 5 } );
		}
		
		bf.main( );
	}
	catch( string strE )
	{
		cout << "Erreur: " << strE << endl;
	}
	
	return 0;
}

int main( ) 
{
	string strDropVertex( "6" ), strDropVertex_connected( "7" ), strFilename( "10-tum04_12_01-k=4-l=4" ), strFolder( "../../../graphs/gbd/" );
	unsigned int iDropVertex( 10 ), iDropVertex_connected( 11 ), iTemp;
	
	for( unsigned int i( 0 ); i < 10; i++ )
	{
		CoxIter ci( strFolder + strFilename + "_" + to_string( i ) + ".coxiter", false, strFolder + strFilename + "_" + to_string( i + 1 ), false, true, false, false);
		if( !ci.readGraph( ) )
		{
			cout << "Erreur lors de la lecture du fichier de l'étape " << ( i + 1 ) << ": " + strFolder + strFilename + "_" + to_string( i ) + ".coxiter\n\t" << ci.get_strError( ) << endl;
			return 0;
		}
		
		/*
		GBD gbd( &ci );
		if( !gbd.removeVertex( strDropVertex ) )
		{
			cout << "Erreur GBD: \n\tEtape: " << i << "\n\tErreur: " << gbd.get_strError() << "\n\tFichier: " + strFolder + strFilename + "_" + to_string( i ) + ".coxiter" << endl;
			return 0;
		}
		else
			cout << "GBD " << ( i + 1 ) << ": OK" << endl;
		
		ci.writeGraph( strFolder + strFilename + "_" + to_string( i + 1 ) );
		strDropVertex += "_" + strDropVertex_connected; // s0
		
		cout << "Drop: " << strDropVertex << ", Connected: " << strDropVertex_connected << endl;*/
		
		cout << "Avant: Drop: " << iDropVertex << ", Connected: " << iDropVertex_connected << endl;
		
		GBD gbd( &ci );
		if( !gbd.removeVertex( to_string( iDropVertex ) ) )
		{
			cout << "Erreur GBD: \n\tEtape: " << i << "\n\tErreur: " << gbd.get_strError() << "\n\tFichier: " + strFolder + strFilename + "_" + to_string( i ) + ".coxiter" << endl;
			return 0;
		}
		else
			cout << "GBD " << ( i + 1 ) << ": OK" << endl;
		
		iTemp = iDropVertex;
		
		iDropVertex =  ci.get_iVertexIndex( to_string( iDropVertex ) + "_" + to_string( iDropVertex_connected ) );
		if( iTemp < iDropVertex_connected )
			iDropVertex_connected--;
		
		//iDropVertex_connected = ci.get_iVertexIndex( to_string( iDropVertex_connected ) );
		ci.map_vertices_labels_reinitialize();
		
		ci.writeGraph( strFolder + strFilename + "_" + to_string( i + 1 ) );
		
		cout << "Next: Drop: " << iDropVertex << ", Connected: " << iDropVertex_connected << endl;
		
		cout << endl;
	}
	
	return 0;
}
