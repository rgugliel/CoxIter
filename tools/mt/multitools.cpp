#include "multitools.h"

bool MultiTools::bReadGraphsList( const string& strList )
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
		if( regexp.preg_match_all( "([^[:space:]]+)", szLine, regexpRes ) == 0 ) // si on ne peut pas lire le nom du fichier
		{
			cout << "Ligne non lue: " << szLine << endl;
			continue;
		}
		
		strGraphsFilename.push_back( regexpRes[1][0] );
	}
	
	fileIn.close( );
	
	return true;
}

void MultiTools::createCycles_main( )
{
	// ---------------------------------------------
	// Test et constructions préliminaires
	unsigned int i, j;
	CoxIter ci( "../../../graphs/Tumarkin pyramids n+2/3-tum04_07_01-k=2-l=3-m=3-n=4.coxiter", false, "", false, true, true, false );
	Arithmeticity arithmeticity;
	
	if( !ci.runAllComputations( ) )
	{
		cout << "Erreur: " << ci.get_strError( ) << endl;
		return;
	}
	
	ci.readGraph();
	
	arithmeticity.test( ci, false );
	if( ci.get_iIsArithmetic( ) != 1 )
		cout << "ATTENTION: non-arithmétique" << endl;
	
	iMax = ci.get_iVerticesCount( );
	vector< vector< unsigned int > > iCoxeterMatrix( ci.get_iCoxeterMatrix( ) );
	anGramMatix = vector< vector< ArithmeticityNumber > >( vector< vector< ArithmeticityNumber > >( iMax, vector< ArithmeticityNumber >( iMax, ArithmeticityNumber( 0 ) ) ) );
	
	// ------------------------------------------------------
	// Construction des éléments
	for( i = 0; i < iMax; i++ )
	{
		for( j = 0; j <= i; j++ )
		{
			if( i == j )
				anGramMatix[i][i] = ArithmeticityNumber( 1 );
			else
			{
				// We forget the signs
				if( iCoxeterMatrix[i][j] == 3 )
					anGramMatix[i][j] = ArithmeticityNumber( 1, false, false, false );
				else if( iCoxeterMatrix[i][j] == 0 )
					anGramMatix[i][j] = ArithmeticityNumber( 1, false, false, false );
				else if( iCoxeterMatrix[i][j] == 2 )
					anGramMatix[i][j] = ArithmeticityNumber( 0, false, false, false );
				else if( iCoxeterMatrix[i][j] == 4 )
					anGramMatix[i][j] = ArithmeticityNumber( 1, true, false, false );
				else if( iCoxeterMatrix[i][j] == 5 )
					anGramMatix[i][j] = ArithmeticityNumber( 1, false, false, true );
				else if( iCoxeterMatrix[i][j] == 6 )
					anGramMatix[i][j] = ArithmeticityNumber( 1, false, true, false );
				else
					throw( 0 );
				
				anGramMatix[j][i] = anGramMatix[i][j];
			}
		}
	}
	
	for( unsigned int k( 0 ); k <= iMax; k++ )
	{
		cout << k << endl;
		iWeights = vector< unsigned int >( k, 1 );
		createCycles_enumerate( 0 );
	}
	
	sort( cyclesVectors.begin( ), cyclesVectors.end( ) );
	cyclesVectors = vector< CycleVector >( cyclesVectors.begin( ), unique( cyclesVectors.begin( ), cyclesVectors.end( ) ) );
	for( vector< CycleVector >::const_iterator it( cyclesVectors.begin( ) ); it != cyclesVectors.end( ); ++it )
	{
		cout << (*it).an << " * v" << (*it).iVector << endl;
	}
}

void MultiTools::createCycles_enumerate( const unsigned int& iStartingIndex )
{
	if( iStartingIndex == iWeights.size( ) )
		return;
	
	for( unsigned int i( 1 ); i <= iMax; i++ )
	{
		iWeights[ iStartingIndex ] = i;
		
		createCycles_enumerate( iStartingIndex + 1 );
		
		if( iWeights.size( ) == iStartingIndex + 1 )
			createCycles_tovector( );
	}
}

void MultiTools::createCycles_tovector( )
{
	if( anGramMatix[0][ iWeights[0] - 1 ] == 0 )
		return;
	
	ArithmeticityNumber an( anGramMatix[0][ iWeights[0] - 1 ] );
	unsigned int iPrevious( iWeights[0] );

	for( vector< unsigned int >::const_iterator it( iWeights.begin( ) + 1 ); it != iWeights.end( ); ++it )
	{
		if( anGramMatix[iPrevious - 1][*it - 1] == 0 )
			return;
		
		an *= anGramMatix[iPrevious - 1][*it - 1];
		iPrevious = *it;
	}
	
	an.reduceIntergerPart( );
	cyclesVectors.push_back( CycleVector( an, iWeights[ iWeights.size( ) - 1 ] ) );
	
	// TODO
	/*
	if( iWeights[ iWeights.size( ) - 1 ] == 1 &&  )
	{
		cout << "Cycle: ";
		for( auto p : iWeights )
			cout << ( p + 1 ) << ", ";
		cout << endl;
	}*/
}

vector< string > MultiTools::get_strGraphsFilename( ) const
{
	return strGraphsFilename;
}

CycleVector::CycleVector( ArithmeticityNumber an, unsigned int iVector )
: an( an ), iVector( iVector )
{

}

bool CycleVector::operator==( const CycleVector& cv )
{
	return ( an == cv.an && iVector == cv.iVector );
}

bool operator<( const CycleVector& cv1, const CycleVector& cv2 )
{
	if( cv1.iVector < cv2.iVector )
		return true;
	else if( cv1.iVector > cv2.iVector )
		return false;
	
	if( cv1.an.parti < cv2.an.parti )
		return true;
	else if( cv2.an.parti > cv2.an.parti )
		return false;
	
	if( !cv1.an.part2 && cv2.an.part2 )
		return true;
	else if( cv2.an.part2 && !cv2.an.part2 )
		return false;
	
	if( !cv1.an.part3 && cv2.an.part3 )
		return true;
	
	return false;
}