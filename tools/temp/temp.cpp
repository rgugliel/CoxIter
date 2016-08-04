/*
Copyright (C) 2013, 2014
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

#include "temp.h"

void Temp::compareWithTerragni()
{
	PCRERegexp regexp;
	PCREResult regexpRes;
	string strLine;
	
	map< string, string > strCorrespondance;
	
	string strGroup, strGrowthSeries, strGrowthRate;
	
	// ----------------------------------------------------------------------------------
	// correspondence
	ifstream fileIn( "../../../documents/growth/correspondence terragni" );
	if( fileIn.fail( ) )
		throw( string( "Error opening file" ) );
	
	while( getline( fileIn, strLine ) ) // tant que....
	{
		regexpRes.clear( );
		if( regexp.preg_match_all( "(hc|hnc)([0-9]+)[[:space:]]+([[:alnum:]_-]+)", strLine, regexpRes ) == 1 )
		{
			strCorrespondance[ regexpRes[1][0] + regexpRes[2][0] ] = regexpRes[3][0];
		}
	}
	
	fileIn.close();
	
	// ----------------------------------------------------------------------------------
	// Data
	fileIn.open( "../../../documents/growth/output.magma" );
	if( fileIn.fail( ) )
		throw( string( "Error opening file" ) );

	while( getline( fileIn, strLine ) ) // tant que....
	{
		regexpRes.clear( );
		if( regexp.preg_match_all( "(hc|hnc)([0-9]+)", strLine, regexpRes ) == 1 )
		{
			strGroup = regexpRes[0][0];
		}
		
		//M(x)=(-x^11 - 2*x^10 - 2*x^9 - 2*x^8 - 2*x^7 - 3*x^6 - 3*x^5 - 2*x^4 - 2*x^3 - 2*x^2- 2*x - 1)/(x^11 - 2*x^10 + x^9 - x^7 + 2*x^6 - 2*x^5 + x^4 - x^2 + 2*x - 1)
		
		regexpRes.clear( );
		if( regexp.preg_match_all( "M\\[x_\\]\\:=([^\\n]+)", strLine, regexpRes ) == 1 )
			strGrowthSeries = regexpRes[0][0] + ";";
		
		regexpRes.clear( );
		if( regexp.preg_match_all( "tau\\=([[:digit:].]+)", strLine, regexpRes ) == 1 )
		{
			strGrowthRate = regexpRes[1][0];
			
			CoxIter ci( "../../../graphs/simplices/" + strCorrespondance[ strGroup ] + ".coxiter", false, "", false, true, true, true, true, false, vector< string >( 0 ), "mathematica" );
			if( !ci.runAllComputations() )
				throw( string( ci.get_strError() ) );
			
			GrowthRate_Result grr;
			grr.bComputed = false;
			grr.iPerron = -1;
			grr.iPisot = -1;
			grr.iSalem = -1;
			
			GrowthRate gr;
			grr = gr.grrComputations( ci.get_iGrowthSeries_denominator() );
			
			//ci.printGrowthSeries();
			//cout << "\n" << strGrowthSeries << endl;
			//cout << "Print[\"" << strGroup << "=" << strCorrespondance[ strGroup ] << ": \",Exponent[Together[f[x]-M[x]],x]]" << endl;
			
			
			cout << "Print[\"" << strGroup << "=" << strCorrespondance[ strGroup ] << ": \"," << grr.strGrowthRate << "-" << strGrowthRate << "]" << endl;
			
			//exit(0);
		}
	}
	
	fileIn.close();
}

bool Temp::bReadGraphsList( const string& strList )
{
	// ------------------------------------------------------
	// ouverture fichier
	
	
	ifstream fileIn( strList.c_str( ) );
	if( fileIn.fail( ) )
	{
		strError = "READ_INPUT_FILE";
		return false;
	}
	
	// ------------------------------------------------------
	// lecture
	PCRERegexp regexp;
	PCREResult regexpRes;
	string szLine, szLineWork, szNumber;

	while( getline( fileIn, szLine ) ) // tant que....
	{
		// --------------------------------------------------------
		// nom du fichier à tester
		regexpRes.clear( );
		if( regexp.preg_match_all( "\"([^\"]+)\"", szLine, regexpRes ) == 0 ) // si on ne peut pas lire le nom du fichier
		{
			cout << "Ligne non lue: " << szLine << endl;
			continue;
		}
		
		strGraphsFilename.push_back( regexpRes[1][0] );
	}
	
	fileIn.close( );
	
	return true;
}

void Temp::testArithmeticity()
{
	for( auto strF : strGraphsFilename )
	{
		CoxIter ci( string( "../../../graphs/" ) + strF, false, "", false, true, true, true, true );
		if( !ci.runAllComputations() )
		{
			cout << "Error" << endl;
			return;
		}
		unsigned int iVerticesCount( ci.get_iVerticesCount() );
		
		Arithmeticity arithmeticity;
		arithmeticity.test( ci, true );
		auto strCycles( arithmeticity.get_strListCycles() );
		
		// We remove duplicates
		sort( strCycles.begin( ), strCycles.end( ) );
		strCycles = vector< string >( strCycles.begin( ), unique( strCycles.begin( ), strCycles.end( ) ) );
		
		// Weights
		auto strWeights( ci.get_strWeights() );
		string strTestWeights;
		for( auto it : strWeights )
		{
			unsigned int i( iLinearizationMatrix_row( it.first, iVerticesCount ) ), j( iLinearizationMatrix_col( it.first, iVerticesCount ) );
			
			strTestWeights += ( strTestWeights == "" ? "l" : " l" ) + to_string( i ) + "m" + to_string( j ) + " := " + strPARItoMathematica( it.second ) + ";";
		}
		
		// Tests
		string strTestA( "IntegerQ[" + implode( "] && IntegerQ[", strCycles ) + "]" );
		string strTestQA( "Element[Expand[" + implode( "], Rationals] && Element[Expand[", strCycles ) + "], Rationals]" );
		string strTest( strTestWeights + " Print[If[" + strTestA + ", \"A\", If[" + strTestQA + ", \"QA\", \"NQA\"]]];" );
		
		// Preparing the output
		str_replace( strF, "non-cocompact n+3/", "" ); // cleaning
		str_replace( strF, ".coxiter", "" ); // cleaning
		
		PCRERegexp regexp;
		PCREResult regexp_res;
		if( regexp.preg_match_all( "([[:digit:]]{1,1})-roberts15_", strF, regexp_res, PCRE_CASELESS ) > 0 )
			str_replace( strF, regexp_res[0][0], "" );
		
		if( ci.get_iIsArithmetic() == -1 )
		{
			// Calling Mathematica
			string strCmd( "echo \"2+2\" | sage -q" ), strResult;
			FILE *FileStream;
			char stdbuffer[1024];
			FileStream = popen( strCmd.c_str(), "r" );
			while( fgets( stdbuffer, 1024, FileStream ) != NULL )
				strResult.append(stdbuffer);
			pclose(FileStream);

			cout << ci.get_iDimension() << ";" << strF << ";" << strResult;
		}
		else
		{
			cout << ci.get_iDimension() << ";" << strF << ";" << ( ci.get_iIsArithmetic() == 1 ? "A" : "NA" ) << endl;
		}
	}
}

/*
 * This is a prototype, incomplete and too simple (non recursive among other things) function!
 */
string Temp::strPARItoMathematica(string str)
{
	PCRERegexp regexp;
	PCREResult regexp_res;
	int regexp_iResultCount;
	
	regexp_iResultCount = regexp.preg_match_all( "sqrt\\(([[:alnum:] /]+)\\)", str, regexp_res, PCRE_CASELESS );
	for( int i(0); i < regexp_iResultCount; i++ )
		str_replace( str, regexp_res[0][i], "Sqrt[" + regexp_res[1][i] + "]" );
	
	regexp_res.clear();
	regexp_iResultCount = regexp.preg_match_all( "cos\\(([[:alnum:] /]+)\\)", str, regexp_res, PCRE_CASELESS );
	for( int i(0); i < regexp_iResultCount; i++ )
		str_replace( str, regexp_res[0][i], "Cos[" + regexp_res[1][i] + "]" );
	
	return str;
}


vector< string > Temp::get_strGraphsFilename( ) const
{
	return strGraphsFilename;
}

void Temp::exportGraphs( )
{
	cout << "Nombre de graphes à lire: " << strGraphsFilename.size( ) << endl;
	
	for( vector< string >::const_iterator it( strGraphsFilename.begin( ) ); it != strGraphsFilename.end( ); ++it )
	{
		/*
		CoxIter ci( *it, false, true, ( *it ) + ".out", false, true );
		
		if( !ci.inputRead( ) )
		{
			cout << "Erreur lors de la lecture du fichier: " << ci.get_strError( ) << endl;
			continue;
		}
		
		if( !ci.writeGraph( ) )
			cout << "Erreur lors de l'écriture du fichier: " << ci.get_strError( ) << endl;*/
	}
}

/*
 * But: A partir du graphe créer des graphes en enlevant une des extrémité
 */
void Temp::dropEnds( const string& strFolderOutput )
{
	for( vector< string >::const_iterator it( strGraphsFilename.begin( ) ); it != strGraphsFilename.end( ); ++it )
	{
		string strFilename( *it );
		str_replace( strFilename, ".in", "" );
		str_replace( strFilename, "../", "" );
		str_replace( strFilename, "graphs/", "" );
		
		CoxIter ci( *it, false, "", false, true, false, false );
		vector< vector< unsigned int > > iAdj;
		unsigned int i, j, k, iVerticesCount, iConnectedCount;
		
		if( !ci.readGraph( ) )
		{
			cout << "Erreur: " << *it << " / " << ci.get_strError() << endl;
			continue;
		}
		
		iAdj = ci.get_iCoxeterMatrix( );
		iVerticesCount = ci.get_iVerticesCount( );
		
		for( i = 0; i < iVerticesCount; i++ )
		{
			iConnectedCount = 0;
			for( j = 0; j < iVerticesCount; j++ )
			{
				if( i != j && iAdj[i][j] != 2 )
					iConnectedCount++;
			}
			
			if( iConnectedCount == 1 )
			{
				string strFilenameOut( strFolderOutput + strFilename + "-drop_" + to_string( i + 1 ) + ".in" );
				ofstream out( strFilenameOut.c_str( ) );
				if( !out.is_open( ) )
				{
					cout << "Cannot open: " << strFilenameOut << endl;
					continue;
				}
				
				out << ( iVerticesCount - 1 ) << " " << ci.get_iDimension( ) << endl;
				
				for( j = 0; j < iVerticesCount; j++ )
				{
					if( j == i )
						continue;
					
					for( k = 0; k < j; k++ )
					{
						if( k == i )
							continue;
						
						if( iAdj[j][k] != 2 )
							out << ( j > i ? j : j + 1 ) << " " << ( k > i ? k : k + 1 ) << " " << iAdj[j][k] << endl;
					}
				}
				
				cout << "../../../graphs/tumarkin_simplices/" << strFilename + "-drop_" + to_string( i + 1 ) + ".in\t\t\tnon-cocompact\tnon-fv" << endl;
				
				out.close( );
			}
		}
	}
}