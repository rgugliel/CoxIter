#include <iostream>
#include <chrono>
#include <stdio.h>

#include "temp.h"
#include "../../../tools/regexp.h"
#include "../../../Vino/rationalinteger_invariantsqf.h"
#include "../../../Vino/quadraticinteger.h"

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

vector< int > iDiagQF( string str )
{
	PCRERegexp regexp;
	PCREResult regexp_res;
	
	int iResCount( regexp.preg_match_all( "\\[[[:space:]\\*]* ([[:digit:]\\-\\/]+)", str, regexp_res ) );
	vector< int > iRes;
	
	for( int i(0); i < iResCount; i++ )
	{
		vector< string > strFrac( explode( "/", regexp_res[1][i] ) );
		int iVal( strFrac.size() == 1 ? stoi( strFrac[0] ) : stoi( strFrac[0] ) * stoi( strFrac[1] ) );
		
		iRes.push_back( iRemoveSquareFactors( iVal ) );
	}
	
	sort( iRes.begin(), iRes.end() );
	
	return iRes;
}

int main( )
{
	QuadraticInteger::set_d(5);
	QuadraticInteger qi( -2, 5 );
	auto qiF( qi.qiPrimeDecomposition() );

	cout << qi << endl;
	
	for( auto f : qiF )
		cout << f.first << ", " << f.second << endl;
	
	return 0;
	
	
	/*
	string strCmd( "echo \"Q=QuadraticForm(QQ,6,[ 1,-2,0,-10,0,0,1,-10,0,0,0,10,-10,-10,0,10,0,0,10,-10,10]);\nQ.rational_diagonal_form()\" | sage -q" ), strResult;
	FILE *FileStream;
	char stdbuffer[1024];
	FileStream = popen( strCmd.c_str(), "r" );
	while( fgets( stdbuffer, 1024, FileStream ) != NULL )
		strResult.append(stdbuffer);
	pclose(FileStream);
	
	auto iQF( iDiagQF( strResult ) );
	
	return 0;
	*/
	
	ifstream fileIn( "graphs roberts.csv" );
	string strLine, strFileOut;
	if( fileIn.fail( ) )
		throw( string( "Error opening file" ) );
	
	while( getline( fileIn, strLine ) ) // tant que....
	{
		auto strExpLine( explode( ";", strLine ) );
		if( strExpLine[2] != "A" )
		{
			strFileOut += strLine + "\n";
			continue;
		}
		
		string strCmd( "echo \"" + strExpLine[5] + ";\nQ.rational_diagonal_form()\" | sage -q" ), strResult;
		FILE *FileStream;
		char stdbuffer[1024];
		FileStream = popen( strCmd.c_str(), "r" );
		while( fgets( stdbuffer, 1024, FileStream ) != NULL )
			strResult.append(stdbuffer);
		pclose(FileStream);
		
		auto iQF2( iDiagQF( strResult ) );
		vector< int > iQF1;
		
		strExpLine[4] = implode( ",", iQF2 );
		
		explode( ",", strExpLine[3], iQF1 );
		
		
		InvariantsQF iqf1( iQF1 );
		InvariantsQF iqf2( iQF2 );
		
		strFileOut += implode( ";", strExpLine ) + ";" + iqf1.get_strInvariant() + ";" + iqf2.get_strInvariant() + "\n";
	}
	
	cout << strFileOut << endl;
	
	fileIn.close();
	
	
	return 0;
}