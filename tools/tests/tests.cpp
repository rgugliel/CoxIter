#include "tests.h"

Tests::~Tests( )
{
	// si l'affichage est redirigé dans un fichier
	if( bCoutFile )
	{
		outCout->close( );
		cout.rdbuf( sBufOld ); // restauration du cout vers la console
	}
}

bool Tests::readGraphsFile( string strInputFilename )
{
	// ------------------------------------------------------
	// ouverture fichier
	ifstream fileIn( strInputFilename.c_str( ) );
	if( fileIn.fail( ) )
	{
		szError = "READ_INPUT_FILE";
		return false;
	}
	
	// ------------------------------------------------------
	// lecture
	PCRERegexp regexp;
	PCREResult regexpRes;
	string szLine, szLineWork, szNumber, strEmpty;
	
	size_t iFindPos;
	
	Test test;
	
	while( getline( fileIn, szLine ) ) // tant que....
	{
		strEmpty = szLine;
		str_replace( strEmpty, " ", "" );
		if( strEmpty == "" )
			continue;
		
		test.bTestEuler = false;
		test.bTestFVector = false;
		test.bTestCompacity = false;
		test.bIsFiniteVolume = true;
		test.bTestArithmeticity = false;
		test.bTestGrowthSeries = false;
		test.strGrowthRate = "";
		
		if( ( iFindPos = szLine.find( "#" ) ) != string::npos ) // Forget what is after #
			szLine = szLine.substr( 0, iFindPos );
		
		// --------------------------------------------------------
		// nom du fichier à tester
		regexpRes.clear( );
		if( regexp.preg_match_all( "\"([^\"]+)\"", szLine, regexpRes ) == 0 ) // si on ne peut pas lire le nom du fichier
		{
			cout << "Ligne non lue: " << szLine << endl;
			continue;
		}
		
		test.szFile = regexpRes[1][0];
		szLineWork = szLine.substr( test.szFile.length() + 2 ); // on enlève le nom du fichier
		str_replace( szLineWork, " ", "" );
		
		// --------------------------------------------------------
		// Growth rate TODO: tester
		regexpRes.clear( );
		if( regexp.preg_match_all( "tau\\=([[:digit:].]+)[;]?", szLineWork, regexpRes ) != 0 )
		{
			str_replace( szLineWork, regexpRes[0][0], " " );
			test.strGrowthRate = regexpRes[1][0];
		}
		
		// --------------------------------------------------------
		// Growth series
		regexpRes.clear( );
		if( regexp.preg_match_all( "f\\(x\\)=[C\\(]*([0-9,]*)[\\)]?/\\(([0-9\\-\\+\\*\\^x]+)\\)[;|\n]{1,1}", szLineWork, regexpRes ) )
		{
			test.bTestGrowthSeries = true;
			str_replace( szLineWork, regexpRes[0][0], "" );
			
			unsigned int iPower;
			test.growthSeries_iPolynomialDenominator = vector< mpz_class >( 1, 0 );
			test.growthSeries_iCyclotomicNumerator.clear( );
			
			vector< unsigned int >  iTemp;
			
			if( regexpRes[1][0] != "" ) // Cyclotomic factors
				explode( ",", regexpRes[1][0], test.growthSeries_iCyclotomicNumerator );
			
			string strPol( regexpRes[2][0] );
			
			regexpRes.clear( );
			int iResCount( regexp.preg_match_all( "(\\+|\\-?)([0-9]*)[\\*]?x[\\^]?([0-9]*)", strPol, regexpRes ) );
			
			for( int i(iResCount - 1); i >= 0; i-- )
			{
				iPower = regexpRes[3][i] == "" ? 1 : abs( stoi( regexpRes[3][i] ) ); // power
				
				if( iPower + 1 > test.growthSeries_iPolynomialDenominator.size( ) )
					test.growthSeries_iPolynomialDenominator.insert( test.growthSeries_iPolynomialDenominator.end( ), iPower - test.growthSeries_iPolynomialDenominator.size() + 1, mpz_class( 0 ) );
				
				test.growthSeries_iPolynomialDenominator[iPower] = ( regexpRes[1][i] == "-" ? -1 : 1 ) * ( regexpRes[2][i] == "" ? 1 : stoi( regexpRes[2][i] ) );
				
				string strRep( regexpRes[1][i] + regexpRes[2][i] + ( regexpRes[2][i] != "" && stoi( regexpRes[2][i] ) > 1 ? "*" : "" ) + "x" + ( iPower > 1 ? "^" : "" ) + regexpRes[3][i] );
				
				str_replace( strPol, strRep, "" );
			}
			
			if( strPol != "" )
				test.growthSeries_iPolynomialDenominator[0] = stoi( strPol );
		}

		// --------------------------------------------------------
		// f-vecteur
		regexpRes.clear( );
		if( regexp.preg_match_all( "\\(([[:digit:] ,]+)\\)", szLineWork, regexpRes ) != 0 )
		{
			str_replace( szLineWork, "(" + regexpRes[1][0] + ")", " " );
			vector< string > strFVector( explode( ",", regexpRes[1][0] ) );
			test.iFVector = vector< unsigned int >( 0 );
			test.bTestFVector = true;
			
			for( vector< string >::const_iterator strIt( strFVector.begin( ) ); strIt != strFVector.end( ); ++strIt )
				test.iFVector.push_back( stoi( *strIt ) );
		}
		
		// --------------------------------------------------------
		// caractéristique d'Euler
		regexpRes.clear( );
		if( regexp.preg_match_all( "([-]{0,1})([[:digit:]/]+)", szLineWork, regexpRes ) != 0 )
		{
			szNumber = regexpRes[1][0] + regexpRes[2][0];
			test.bTestEuler = true;
			
			if( regexp.preg_match_all( "([[:digit:]]+)", szNumber, regexpRes ) != 0 ) // si l'expression contient au moins un chiffre
			{
				try{
					test.brResult = Number_rational( szNumber );
				}
				catch( int iCode )
				{
					cout << "Error\t " << szLine << endl;
					cout << "\t\tErrreur lors de la lecture de la caractéristique" << endl;
					cout << "\t\t" << szNumber << " n'est pas reconnu comme un nombre valide" << endl;
					test.bTestEuler = false;
				}
			}
			else
				test.bTestEuler = false;
			
		}
		
		// --------------------------------------------------------
		// cocompacity
		if( szLineWork.find( "compact" ) != string::npos || szLineWork.find( "cocompact" ) != string::npos )
		{
			test.bTestCompacity = true;
			test.bIsCompact = true;
		}
		if( szLineWork.find( "non-compact" ) != string::npos || szLineWork.find( "non-cocompact" ) != string::npos )
		{
			test.bTestCompacity = true;
			test.bIsCompact = false;
		}
		
		// --------------------------------------------------------
		// arithmeticity
		if( szLineWork.find( "arithmetic" ) != string::npos )
		{
			test.bTestCompacity = true;
			test.bTestArithmeticity = true;
			test.bIsArithmetic = true;
		}
		if( szLineWork.find( "non-arithmetic" ) != string::npos )
		{
			test.bTestCompacity = true;
			test.bTestArithmeticity = true;
			test.bIsArithmetic = false;
		}
		
		// --------------------------------------------------------
		// covolume fini
		if( szLineWork.find( "non-fv" ) != string::npos )
		{
			test.bIsFiniteVolume = false;
		}
		
		tests_EulerCarac.push_back( test );
	}
	
	fileIn.close( );
	
	// ------------------------------------------------------
	// redirection du cout
	bCoutFile = true;
	string str_tempFilename( strInputFilename + ".output" );
	outCout = new ofstream( str_tempFilename.c_str( ) );
	
	if( !outCout->is_open( ) )
		this->bCoutFile = false;
	else
		sBufOld = cout.rdbuf( outCout->rdbuf( ) );
	
	return true;
}

void Tests::runTests()
{
	CoxIter ci;

	bool bCanBeFiniteCovolume;
	int i, iMax( tests_EulerCarac.size( ) ), iCompacity, iFiniteVolume, iArithmeticity, iSignatureComputed;
	unsigned int iDim;
	string strGrowthRate;
	array< unsigned int, 3 > iSignature;
	
	//#pragma omp parallel for private(ci, i, iCompacity, iFiniteVolume, iArithmeticity, bCanBeFiniteCovolume, iDim, strGrowthRate, iSignatureComputed, iSignature) shared(iMax) schedule(static,1)
	for( i = 0; i < iMax; i++ )
	{
		ci = CoxIter( "../../../graphs/" + tests_EulerCarac[i].szFile, false, "", false, true, tests_EulerCarac[i].bTestCompacity, true, true );
		if( !ci.readGraph( ) )
		{
			cout << "Error\t " << tests_EulerCarac[i].szFile << endl;
			cout << "\t\tErreur lors de l'ouverture du fichier (" << ci.get_strError( ) << ")" << endl;
			continue;
		}
		
		iDim = ci.get_iDimension();
		ci.set_iDimension(0); // TODO: remove if doesn't want to check the guessing of dimension
		
		ci.exploreGraph( );
		ci.computeGraphsProducts( );

		if( !ci.euler( ) )
		{
			#pragma omp critical
			{
				cout << "Error\t " << tests_EulerCarac[i].szFile << endl;
				cout << "\t\tProblème avec l'encodage du graphe ou son covolume n'est pas fini" << endl;
			}
		}
		else
		{
			if( tests_EulerCarac[i].bTestCompacity )
				ci.iIsGraphCocompact( );
			
			// Finite volume
			iFiniteVolume = ci.isFiniteCovolume( );
			bCanBeFiniteCovolume = ci.bCanBeFiniteCovolume( );
			ci.growthSeries( );
			
			if( tests_EulerCarac[i].bTestArithmeticity )
			{
				Arithmeticity arithmeticity;
				arithmeticity.test( ci, false );
			}
			
			#pragma omp critical
			{
				if( tests_EulerCarac[i].strGrowthRate != "" ) // TODO: move outside
				{
					GrowthRate* gr( new GrowthRate() );
					GrowthRate_Result grr( gr->grrComputations( ci.get_iGrowthSeries_denominator() ) );
					delete gr;
					strGrowthRate = grr.strGrowthRate;
				}
			
				try{ // TODO: move outside
					if( ci.get_iHasDottedLineWithoutWeight() == 0 )
					{
						//cout << "Test: " << tests_EulerCarac[i].szFile << endl;
						Signature s;
						iSignature = s.iComputeSignature( ci.get_strGramMatrix_PARI() );
						iSignatureComputed = 1;
						//cout << "Fin signature" << endl;
					}
					else
						iSignatureComputed = 0;
				}
				catch( string strE )
				{
					iSignatureComputed = -1;
				}
				
				
				if( tests_EulerCarac[i].bTestEuler )
				{
					if( ci.get_brEulerCaracteristic( ) != tests_EulerCarac[i].brResult )
					{
						cout << "Error Euler\t " << tests_EulerCarac[i].szFile << endl;
						cout << "\t\tExpected: " << tests_EulerCarac[i].brResult << endl;
						cout << "\t\tComputed: " << ci.get_brEulerCaracteristic( ) << endl;
					}
					else
						cout << "OK\tEuler\t\t\t" << tests_EulerCarac[i].szFile << endl;
				}
				
				if( ci.get_iIsFiniteCovolume( ) > 0 )
				{
					if( ci.get_iFVectorAlternateSum( ) == ( ci.get_iDimension( ) % 2 ? 2 : 0 ) )
						cout << "OK\tSomme alt\t\t" << tests_EulerCarac[i].szFile << endl;
					else
					{
						cout << "Error\t\t " << tests_EulerCarac[i].szFile << endl;
						cout << "\tAlt sum: " << ci.get_iFVectorAlternateSum( ) << endl;
					}
				}
				
				if( tests_EulerCarac[i].strGrowthRate != "" )
				{
					if( tests_EulerCarac[i].strGrowthRate == strGrowthRate )
						cout << "OK\tGrowth rate\t\t" << tests_EulerCarac[i].szFile << endl;
					else
					{
						cout << "Error\tGrowth rate\t\t" << tests_EulerCarac[i].szFile << endl;
						cout << "\t\tComputed: " << strGrowthRate << endl;
						cout << "\t\tExpected: " << tests_EulerCarac[i].strGrowthRate << endl;
					}
				}
				
				if( ci.get_bDimensionGuessed() )
				{
					if( iDim == ci.get_iDimension() )
						cout << "OK\tDimension guessed\t" << tests_EulerCarac[i].szFile << endl;
					else
						cout << "Error\tDimension guessed\t" << tests_EulerCarac[i].szFile << endl;
				}
				
				if( iSignatureComputed == 1 && ci.get_iDimension() )
				{
					if( iSignature[0] != ci.get_iDimension() || iSignature[1] != 1 )
					{
						cout << "Error\t " << tests_EulerCarac[i].szFile << endl;
						cout << "\t\tSignature computed: (" << iSignature[0] << "," << iSignature[1] << "," << iSignature[2] << ")" << endl;
					}
					else
						cout << "OK\tSignature\t\t" << tests_EulerCarac[i].szFile << endl;
				}
				else if( iSignatureComputed == -1 )
				{
					cout << "Error\t " << tests_EulerCarac[i].szFile << endl;
					cout << "\t\tSignature: format" << endl;
				}
				else if( iSignatureComputed == 0 && tests_EulerCarac[i].szFile.find( "roberts", 0 ) != string::npos )
				{
					cout << "Error\t " << tests_EulerCarac[i].szFile << endl;
					cout << "\t\tWeights unkown" << endl;
				}
				
				// A small test of the growth series
				if( ci.get_iIsFiniteCovolume( ) > 0  )
				{
					vector< mpz_class > iDenom;
					vector< unsigned int > iCyclotomic;
					bool bReduced;
					
					ci.get_iGrowthSeries( iCyclotomic, iDenom, bReduced );
					
					if( tests_EulerCarac[i].bTestGrowthSeries )
					{
						if( tests_EulerCarac[i].growthSeries_iCyclotomicNumerator == iCyclotomic && tests_EulerCarac[i].growthSeries_iPolynomialDenominator == iDenom )
						{
							cout << "OK\tGrowth series\t\t" << tests_EulerCarac[i].szFile << endl;
						}
						else
						{
							cout << "Error\t " << tests_EulerCarac[i].szFile << endl;
							cout << "\t\tGrowth series" << endl;
						}
					}
					
					mpz_class iTotalDenom( 0 );
					for( unsigned int j( 0 ); j < iDenom.size( ); j++ )
						iTotalDenom += iDenom[j];
					
					if( ci.get_iDimension( ) % 2 ) // n is odd, the denominator should vanish in 1
					{
						if( iTotalDenom == 0 )
							cout << "OK\tGrowth series test\t" << tests_EulerCarac[i].szFile << endl;
						else
						{
							cout << "Error\t " << tests_EulerCarac[i].szFile << endl;
							cout << "\t\tGrowth series test: dénominateur ne s'annule pas en 1" << endl;
						}
					}
					else // n even, we should find the euler characteristics
					{
						mpz_class iTotalNum( 1 ), iSum;
						
						for( auto cyclo : iCyclotomic )
						{
							iSum = 0;
							for( auto coeff : Polynomials::iCyclotomicPolynomials[cyclo] )
								iSum += coeff;
							
							iTotalNum *= iSum;
						}
						
						if( Number_rational( iTotalDenom, iTotalNum ) == ci.get_brEulerCaracteristic( ) )
							cout << "OK\tGrowth series test\t" << tests_EulerCarac[i].szFile << endl;
						else
						{
							cout << "Error\t " << tests_EulerCarac[i].szFile << endl;
							cout << "\t\tGrowth series test: différent de caractéristique d'Euler" << endl;
						}
					}
				}
				
				if( tests_EulerCarac[i].bTestCompacity )
				{
					iCompacity = ci.get_iIsCocompact( );
					if( ( iCompacity == 1 && tests_EulerCarac[i].bIsCompact ) || ( iCompacity == 0 && !tests_EulerCarac[i].bIsCompact ) )
						cout << "OK\tCompacity\t\t" << tests_EulerCarac[i].szFile << endl;
					else
					{
						cout << "Error\t " << tests_EulerCarac[i].szFile << endl;
						
						cout << "\t\tCompacité attendue: " << ( tests_EulerCarac[i].bIsCompact ? "oui" : "non" ) << endl;
						
						cout << "\t\tCompacité calculée: ";
						if( iCompacity == 0 )
							cout << "non" << endl;
						else if( iCompacity == 1 )
							cout << "oui" << endl;
						else
							cout << "?" << endl;
					}
				}
				
				if( tests_EulerCarac[i].bTestArithmeticity )
				{
					iArithmeticity = ci.get_iIsArithmetic( );
					if( ( iArithmeticity == 1 && tests_EulerCarac[i].bIsArithmetic ) || ( iArithmeticity == 0 && !tests_EulerCarac[i].bIsArithmetic ) )
						cout << "OK\tArithmeticity\t\t" << tests_EulerCarac[i].szFile << endl;
					else
					{
						cout << "Error\t " << tests_EulerCarac[i].szFile << endl;
						
						cout << "\t\tExpected arithmeticity: " << ( tests_EulerCarac[i].bIsCompact ? "oui" : "non" ) << endl;
						
						cout << "\t\tComputed arithmeticity: ";
						if( iArithmeticity == 0 )
							cout << "no" << endl;
						else if( iArithmeticity == 1 )
							cout << "yes" << endl;
						else
							cout << "?" << endl;
					}
				}
				
				if( tests_EulerCarac[i].bTestFVector )
				{
					if( ci.get_iFVector( ) == tests_EulerCarac[i].iFVector )
						cout << "OK\tf-vector\t\t" << tests_EulerCarac[i].szFile << endl;
					else
					{
						cout << "Error\t " << tests_EulerCarac[i].szFile << endl;
						
						cout << "\t\tf-vector attendu: ( ";
						vector< unsigned int > iFV( ci.get_iFVector( ) );
						copy( iFV.begin( ), iFV.end( ), ostream_iterator<unsigned int>( cout, " " ) );
						cout << ")\n\t\tf-vector calculé: ( ";
						copy( tests_EulerCarac[i].iFVector.begin( ), tests_EulerCarac[i].iFVector.end( ), ostream_iterator<unsigned int>( cout, " " ) );
						cout << ")" << endl;
					}
				}
				
				if( iFiniteVolume == 1 && !bCanBeFiniteCovolume )
				{
					cout << "Error\t " << tests_EulerCarac[i].szFile << endl;
						
					cout << "\t\tTest CanBeFiniteCovolume" << endl;
				}
				else if( ( iFiniteVolume == 0 && !bCanBeFiniteCovolume ) || ( iFiniteVolume == 1 && bCanBeFiniteCovolume ) )
					cout << "OK\tCanBeFiniteCovolume\t" << tests_EulerCarac[i].szFile << endl;
				
				if( !( ( iFiniteVolume == 1 && tests_EulerCarac[i].bIsFiniteVolume ) || ( iFiniteVolume == 0 && !tests_EulerCarac[i].bIsFiniteVolume ) ) )
				{
					cout << "Error\t " << tests_EulerCarac[i].szFile << endl;
						
					cout << "\t\tTest covolume fini attendu: " << ( tests_EulerCarac[i].bIsFiniteVolume ? "oui" : "non" ) << "\n\t\tTest covolume fini calculé: ";
					if( iFiniteVolume == 0 )
						cout << "non" << endl;
					else if( iFiniteVolume == 1 )
						cout << "oui" << endl;
					else
						cout << "?" << endl;
				}
				else
					cout << "OK\tTest covolume fini\t" << tests_EulerCarac[i].szFile << endl;
			} // #pragma omp critical
		}
	}
}

vector< Test > Tests::get_tests() const
{
	return tests_EulerCarac;
}