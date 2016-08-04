#ifndef BRUTEFORCE_H
#define BRUTEFORCE_H

using namespace std;

#include <vector>
#include <iostream>
#include <iterator>
#include <omp.h>

#include "../../coxiter.h"

class Edge
{
	public:
		unsigned int i;
		unsigned int j;
		
		vector< unsigned int > iPossibleValues;
		unsigned int iPossibleValuesCount;
		
	public:
		Edge( unsigned int i, unsigned int j, vector< unsigned int > iPossibleValues );
};

typedef vector< vector< unsigned int > > Matrix;

class BruteForce
{
	private:
		vector< vector< unsigned int > > iWeights;
		
		vector< Edge > edges;
		unsigned int iEdgesCount;
		
		unsigned int iDimension;
		unsigned int iMatrixSize;
		
		Matrix iBasisMatrix;
		vector< Matrix > iMatrices;
		
	public:
		BruteForce( unsigned int iDimension );
		
		void basisMatrix_create( unsigned int iVerticesNumber );
		void basisMatrix_addEdge( unsigned int i, unsigned int j, unsigned int iWeight );
		void variableEdge_add( unsigned int i, unsigned int j, vector< unsigned int > iWeights );
		
		void main( );
		void enumerate( const unsigned int& iStartingIndex, const unsigned int& iWorkingIndex );
		void enumerate_main( );
		
		bool isConnected( unsigned int iIndex );
		void isConnected_sub( const unsigned int& iIndex, const unsigned int& iStartFrom, vector< bool >& bVisitedVertices );
};

#endif // BRUTEFORCE_H
