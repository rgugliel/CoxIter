#include <iostream>
#include <chrono>
#include <stdio.h>

#include "../../arithmeticity.h"
#include "../../gbd.h"
#include "../../../tools/polynomials.h"
#include "../../../Vino/quadraticinteger.h"
#include "temp.h"

#include <pari/pari.h>

#include <omp.h>

using namespace std;

GEN vector2t_VEC( vector< int > iCoefficients ) // [i] = x^i
{
	unsigned int iSize( iCoefficients.size() );
	const long iVar = 0; // Polynomial variable
	
	GEN pol( cgetg(iSize+2, t_POL) ); // iSize coefficients, sign and variable number
	pol[1] = evalsigne(1) | evalvarn(iVar) | evallgef(iSize+2); // not equal to zero, variable, coefficients
	
	for( unsigned int i(0); i < iSize; i++ )
		pol[i+2] = (long)stoi( iCoefficients[i] );
	
	return pol;
}

vector< long int > t_VEC2vector( GEN poly )
{
	unsigned int iSize( lg(poly) );
	vector< long int > iV;
	
	for( unsigned int i(2); i < iSize; i++ )
		iV.push_back( gtolong(gel(poly,i)) );
	
	return iV;
}

void factorInFp( vector< int > iCoefficients, unsigned int iPrime )
{
	pari_sp av = avma;
	
	GEN pol( vector2t_VEC( { 1, 1, 2, 1, 1 } ) ), gPrime( stoi( iPrime ) );
	
	pol = RgX_to_FpX( pol, gPrime ); // Reduction mod gPrime
	GEN gFactors( FpX_factor( pol, gPrime ) ), gExponents;
	
	gExponents = gel(gFactors,2); // Exponents of the irreducible factors
	gFactors = gel(gFactors,1); // Irreducible factors
	long k( lg(gFactors) - 1 );
	
	for( unsigned int i(1); i <= k; i++ )
		vector< long int > vec( t_VEC2vector( gel(gFactors,i) ) );
	
	avma = av;
}

void factorInFp_( vector< int > iCoefficients, unsigned int iPrime )
{
	GEN pol( vector2t_VEC( iCoefficients ) ), gPrime( stoi( iPrime ) );
	auto test( t_VEC2vector( pol ) );
	
	cout << "Factorisation de ";
	Polynomials::polynomialDisplay( test );
	cout << " sur F_" << gtolong( gPrime ) << endl;
	
	pol = RgX_to_FpX( pol, gPrime ); // Reduction mod gPrime
	GEN gFactors( FpX_factor( pol, gPrime ) ), gExponents;
	
	gExponents = gel(gFactors,2); // Exponents of the irreducible factors
	gFactors = gel(gFactors,1); // Irreducible factors
	long k( lg(gFactors) - 1 );
	
	cout << "Nombre de facteurs: " << k << endl;
	for( unsigned int i(1); i <= k; i++ )
	{
		vector< long int > vec( t_VEC2vector( gel(gFactors,i) ) );
		
		cout << "(";
		Polynomials::polynomialDisplay( vec );
		cout << ")^" << gExponents[i] << endl;
	}
}

int main( )
{
	string strCmd( "echo \"2+2\" | sage -q" ), strResult;
	FILE *FileStream;
	char stdbuffer[1024];
	FileStream = popen( strCmd.c_str(), "r" );
	while( fgets( stdbuffer, 1024, FileStream ) != NULL )
		strResult.append(stdbuffer);
	pclose(FileStream);
	
	cout << "Res: \n" << strResult << endl;
}

int main_terragni( )
{
	Temp t;
	
	try{
	t.compareWithTerragni();
	}
	catch( string strE )
	{
		cout << "Error: " << strE << endl;
	}
	return 0;
}

int main__( )
{
	/*
	QuadraticInteger::set_d(5);
	vector< QuadraticInteger > qiS;
	
	QuadraticInteger qiA( 2, 1 ), qiB( -1, 0 ), qiP( 2, 0 );
	
	cout << "P= " << qiP << ", Norm=" << qiP.iNorm() << endl;
	cout << "a=" << qiA << ", Norm=" << qiA.iNorm() << ", val=" << qiA.iValuation( qiP ) << endl;
	cout << "b=" << qiB << ", Norm=" << qiB.iNorm() << ", val=" << qiB.iValuation( qiP ) << endl;

	unsigned int m( 3 ), p( 5 );
	
	unsigned int i, j, k;
	int iVal;
	
	for( i = 0; i < 25; i++ )
	{
		for( j = 0; j < 5; j++ )
			qiS.push_back( QuadraticInteger( i - j, j ) );
	}
	
	unsigned int iQISSize( qiS.size() );
	
	unsigned int iDCountTot( 0 ), iDCountVal( 0 ), iDCountMax( 0 ); 
	
	for( i = 0; i < iQISSize; i++ )
	{
		QuadraticInteger X( qiS[i] );
		for( j = 0; j < iQISSize; j++ )
		{
			QuadraticInteger Y( qiS[j] );
			for( k = 0; k < iQISSize; k++ )
			{
				iDCountTot++;
				QuadraticInteger Z( qiS[k] );
				
				QuadraticInteger qiTemp( 0 ), qiSum( 0 );
				
				qiSum = qiA;
				qiSum.multiplyBy( &X );
				qiSum.multiplyBy( &X );
				
				qiTemp = qiB;
				qiTemp.multiplyBy( &Y );
				qiTemp.multiplyBy( &Y );
				qiSum.add( &qiTemp );
				
				qiTemp = Z;
				qiTemp.multiplyBy( &Z );
				
				qiSum.substract( &qiTemp );
				
				iVal = qiSum.iValuation( qiP );
				if( iVal >= 0 && iVal < m )
					continue;
				
				iDCountVal++;
				
				if( X.iValuation( qiP ) != 0 && Y.iValuation( qiP ) != 0 && Z.iValuation( qiP ) != 0 )
					continue;
				
				iDCountMax++;
				
				cout << "X=" << X << ", val=" << X.iValuation( qiP ) << endl;
				cout << "Y=" << Y << ", val=" << Y.iValuation( qiP ) << endl;
				cout << "Z=" << Z << ", val=" << Z.iValuation( qiP ) << endl;
				cout << "Sum= " << qiSum << ", val=" << qiSum.iValuation( qiP ) << endl;
				
				exit(0);
			}
		}
	}
	
	cout << "Total: " << iDCountTot << endl;
	cout << "Valuation: " << iDCountVal << endl;
	cout << "Max: " << iDCountMax << endl;
	
	return 0;*/
}

int main_( ) 
{
	chrono::time_point< std::chrono::system_clock > timeStart, timeEnd;
	unsigned int iMax( 100 );
	
	/*
	timeStart = chrono::system_clock::now();
	for( unsigned int i(0); i < iMax; i++ )
	{
		string cmd = "echo \"G:={{1,2},{2,1}}; pos := 0; neg := 0; eig = Eigenvalues[G]; For[i = 1, i <= Dimensions[G][[1]], i++, If[eig[[i]] > 0, pos = pos + 1,  If[eig[[i]] < 0, neg = neg + 1, Blank[]]]]; Print[pos,\\\",\\\",neg];\" | math -noprompt";
		string OutString;
		FILE *FileStream;
		char stdbuffer[1024];
		FileStream = popen(cmd.c_str(), "r");
		while (fgets(stdbuffer, 1024, FileStream) != NULL)
			OutString.append(stdbuffer);
		pclose(FileStream);
	}*/
	
	
	string cmd = "echo \"G:={{ 1, -Sqrt[2]/2, 0, 0, 0, 0}, {-Sqrt[2]/2, 1, -1/2, 0, 0, 0}, {0, -1/2, 1, -Sqrt[2]/2, 0, l25}, {0, 0, -Sqrt[2]/2, 1, l34, 0}, {0, 0, 0, l34, 1, 0}, {0, 0, l25, 0, 0, 1}}; f[x_] = CharacteristicPolynomial[G, x]; Reduce[Coefficient[f[x], x, 0] == 0 && l34 < -1 && l25 < -1, {l34, l25}] \" | math -noprompt";
	string OutString;
	FILE *FileStream;
	char stdbuffer[1024];
	FileStream = popen(cmd.c_str(), "r");
	while (fgets(stdbuffer, 1024, FileStream) != NULL)
		OutString.append(stdbuffer);
	pclose(FileStream);
	cout << OutString << endl;
	
	
	timeEnd = chrono::system_clock::now();
	cout << "\tComputation time: " << chrono::duration <double, milli>(timeEnd-timeStart).count( ) / ( 1000 * iMax ) << "\n" << endl;
	return 0;
}