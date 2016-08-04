#include "bruteforce.h"

Edge::Edge(unsigned int i, unsigned int j, std::vector< unsigned int > iPossibleValues)
: i(i), j(j), iPossibleValues( iPossibleValues )
{
	iPossibleValuesCount = iPossibleValues.size( );
}

BruteForce::BruteForce( unsigned int iDimension )
: iDimension( iDimension )
{
}

void BruteForce::basisMatrix_create( unsigned int iVerticesNumber )
{
	iBasisMatrix = vector< vector< unsigned int > >( iVerticesNumber,  vector< unsigned int >( iVerticesNumber, 2 ) );
	iMatrixSize = iVerticesNumber;
}

void BruteForce::basisMatrix_addEdge( unsigned int i, unsigned int j, unsigned int iWeight )
{
	if( i >= iBasisMatrix.size( ) || j >= iBasisMatrix.size( ) )
		throw( (string) ("INVALID COORDINATES: (" + to_string( i ) + "," + to_string( j ) + ")" ) );
	
	if( i == j )
		throw( (string) ("INVALID COORDINATES: SAME" ) );
	
	iBasisMatrix[i][j] = iWeight;
	iBasisMatrix[j][i] = iWeight;
}

void BruteForce::variableEdge_add(unsigned int i, unsigned int j, vector< unsigned int > iWeights )
{
	if( i >= iBasisMatrix.size( ) || j >= iBasisMatrix.size( ) )
		throw( (string) ("INVALID COORDINATE: (" + to_string( i ) + "," + to_string( j ) + ")") );
	
	edges.push_back( Edge( i, j, iWeights ) );
}


void BruteForce::main( )
{
	iEdgesCount = edges.size( );
	iWeights.push_back( vector< unsigned int >( iEdgesCount, 0 ) );
	
	for( unsigned int i( 0 ); i < iEdgesCount; i++ )
	{
		iBasisMatrix[ edges[i].i ][ edges[i].j ] = edges[i].iPossibleValues[0];
		iBasisMatrix[ edges[i].j ][ edges[i].i ] = edges[i].iPossibleValues[0];
		
		iWeights[0][i] = edges[i].iPossibleValues[0];
	}
	
	for( unsigned int i( 1 ); i < (unsigned int) omp_get_max_threads( ); i++ )
		iWeights.push_back( iWeights[0] );
	
	for( unsigned int i( 0 ); i < (unsigned int) omp_get_max_threads( ); i++ )
		iMatrices.push_back( iBasisMatrix );
	
	enumerate_main( );
}

void BruteForce::enumerate_main( )
{
	unsigned int i;
	
	#pragma omp parallel for private(i) schedule(static,1)
	for( i = 0; i < edges[0].iPossibleValuesCount; i++ )
	{
		iMatrices[ omp_get_thread_num() ][ edges[0].i ][ edges[0].j ] = edges[0].iPossibleValues[i];
		iMatrices[ omp_get_thread_num() ][ edges[0].j ][ edges[0].i ] = edges[0].iPossibleValues[i];
		iWeights[ omp_get_thread_num() ][0] = edges[0].iPossibleValues[i];
		
		enumerate( 1, omp_get_thread_num() );
	}
}

void BruteForce::enumerate( const unsigned int& iStartingIndex, const unsigned int& iWorkingIndex )
{
	if( iEdgesCount == iStartingIndex )
		return;
	
	for( unsigned int i( 0 ); i < edges[iStartingIndex].iPossibleValuesCount; i++ )
	{
		iWeights[ iWorkingIndex ][ iStartingIndex ] = edges[iStartingIndex].iPossibleValues[i];
		
		iMatrices[ iWorkingIndex ][ edges[iStartingIndex].i ][ edges[iStartingIndex].j ] = edges[iStartingIndex].iPossibleValues[i];
		iMatrices[ iWorkingIndex ][ edges[iStartingIndex].j ][ edges[iStartingIndex].i ] = edges[iStartingIndex].iPossibleValues[i];
		
		enumerate( iStartingIndex + 1, iWorkingIndex );
		
		if( iEdgesCount == iStartingIndex + 1 )
		{
			CoxIter ci( iMatrices[ iWorkingIndex ], iDimension, true, true );
			ci.exploreGraph( );
			ci.computeGraphsProducts( );
			
			if( !ci.euler( )  )
				continue;
			
			if( !isConnected( iWorkingIndex ) )
				continue;
			
			ci.writeGraph( "bf/" + implode( "", iWeights[ iWorkingIndex ] ) );
		}
	}
}

bool BruteForce::isConnected( unsigned int iIndex )
{
	vector< bool > bVisitedVertices( iMatrixSize, false );
	
	bVisitedVertices[0] = true;
	for( unsigned int i( 1 ); i < iMatrixSize; i++ )
	{
		if( iMatrices[iIndex][0][i] != 2 )
			isConnected_sub( iIndex, i, bVisitedVertices );
	}
	
	for( unsigned int i( 0 ); i < iMatrixSize; i++ )
	{
		if( !bVisitedVertices[i] )
			return false;
	}
	
	return true;
}

void BruteForce::isConnected_sub( const unsigned int& iIndex, const unsigned int& iStartFrom, vector< bool >& bVisitedVertices )
{
	bVisitedVertices[iStartFrom] = true;
	
	for( unsigned int i( 0 ); i < iMatrixSize; i++ )
	{
		if( i != iStartFrom && !bVisitedVertices[i] && iMatrices[iIndex][0][i] != 2 )
			isConnected_sub( iIndex, i, bVisitedVertices );
	}
}