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

#include "app.h"

App::App( )
: bCoutFile( false ), 
bOutputGraphToDraw( false ), 
bOutputGraph( false ),
bCheckCanBeFiniteCovolume( false ),
bCheckCocompacity( false ), 
bCheckFiniteCovolume( false ), 
bCheckArithmeticity( false ),
bComputeEuler( true ),
bComputeGrowthRate( false ),
bComputeGrowthSeries( false ),
bComputeSignature( false ),
bDoComputations( true ), 
bDebug( false ), 
bGBD( false ), 
bOpenMP( true ),
bPrintCoxeterGraph( false ),
bPrintCoxeterMatrix( false ), 
bPrintGramMatrix( false ),
bPrintHelp( false ),
strOuputMathematicalFormat( "generic" )
{

}

bool App::bReadMainParameters( int argc, char **argv )
{
	string szTemp, szPrevType, strComplete;
	
	if( argc == 1 && isatty(fileno(stdin)) == 0 ) // If no parameter was given but a file was passed via stdin
	{
		bCheckCocompacity = true;
		bCheckFiniteCovolume = true;
		bComputeGrowthSeries = true;
#ifdef _COMPILE_WITH_PARI_
			bComputeGrowthRate = true;
			bComputeSignature = true;
#endif
		return true;
	}
	
	for( int i = 0; i < argc; ++i ) 
	{
		szTemp = std::string( argv[i] );
		
		strComplete += szTemp + " ";
		
		if( szTemp == "-i" ) // input
			szPrevType = "i";
		else if( szTemp == "-arithmeticity" || szTemp == "-a" )
		{
			bCheckCocompacity = true;
			bCheckArithmeticity = true;
			szPrevType = "arithmeticity";
		}
		else if( szTemp == "-c" || szTemp == "-compacity" || szTemp == "-compactness" || szTemp == "-compact" || szTemp == "-cocompact" )
		{
			bCheckCocompacity = true;
			szPrevType = "compacity";
		}
		else if( szTemp == "-cf" || szTemp == "-of" ) // force le cout dans un fichier
		{
			bCoutFile = true;
			szPrevType = "cf";
		}
		else if( szTemp == "-debug" )
		{
			bDebug = true;
			szPrevType = "debug";
		}
		else if( szTemp == "-drawgraph" || szTemp == "-dg" ) // output
		{
			bOutputGraphToDraw = true;
			szPrevType = "drawgraph";
		}
		else if( szTemp == "-drop" || szTemp == "-remove" ) // we remove a vertex
			szPrevType = "drop";
		else if( szTemp == "-ffv" )
		{
			bCheckCanBeFiniteCovolume = true;
			szPrevType = "ffv";
		}
		else if( szTemp == "-fcv" )
		{
			bCheckFiniteCovolume = true;
			szPrevType = "fv";
		}
		else if( szTemp == "-full" )
		{
			bCheckCocompacity = true;
			bCheckFiniteCovolume = true;
			bComputeGrowthSeries = true;
#ifdef _COMPILE_WITH_PARI_
			bComputeGrowthRate = true;
			bComputeSignature = true;
#endif
			szPrevType = "full";
		}
		else if( szTemp == "-fv" )
		{
			bCheckFiniteCovolume = true;
			szPrevType = "fv";
		}
		else if( szTemp == "-g" || szTemp == "-growth" || szTemp == "-poincarre" )
		{
			bComputeGrowthSeries = true;
			szPrevType = "growth";
		}
		else if( szTemp == "-growthrate" || szTemp == "-gr" )
		{
#ifdef _COMPILE_WITH_PARI_
			bComputeGrowthSeries = true;
			bComputeGrowthRate = true;
			szPrevType = "growthrate";
#endif
		}
		else if( szTemp == "-gbd" )
		{
			bGBD = true;
			szPrevType = "gbd";
		}
		else if( szTemp == "-help" )
		{
			bPrintHelp = true;
			szPrevType = "help";
		}
		else if( szTemp == "-nc" ) // don't do the computations using BigNumbers
		{
			bDoComputations = false;
			szPrevType = "nc";
		}
		else if( szTemp == "-ne" || szTemp == "-neuler" ) // don't compute the Euler characteristic
		{
			bComputeEuler = false;
			szPrevType = "ne";
		}
		else if( szTemp == "-nopenmp" || szTemp == "-nparallel" )
		{
			bOpenMP = false;
			szPrevType = "npenmp";
		}
		else if( szTemp == "-o" ) // Files for the output
		{
			szPrevType = "o";
		}
		else if( szTemp == "-pcg" ) // print Coxeter matrix
		{
			bPrintCoxeterGraph = true;
			szPrevType = "pcg";
		}
		else if( szTemp == "-pcm" ) // print Coxeter matrix
		{
			bPrintCoxeterMatrix = true;
			szPrevType = "pcm";
		}
		else if( szTemp == "-pgm" ) // print Gram matrix?
		{
			bPrintGramMatrix = true;
			szPrevType = "pgm";
		}
		else if( szTemp == "-oformat" || szTemp == "-outputformat" )
		{
			szPrevType = "oformat";
		}
		else if( szTemp == "-s" || szTemp == "-signature" )
		{
#ifdef _COMPILE_WITH_PARI_
			bComputeSignature = true;
			szPrevType = "signature";
#endif
		}
		else if( szTemp == "-writegraph" || szTemp == "-wg" ) // write the graph
		{
			bOutputGraph = true;
			szPrevType = "wg";
		}
		else
		{
			if( szPrevType == "i" )
				strInFilename = szTemp;
			else if( szPrevType == "drop" )
				strVerticesRemove.push_back( szTemp );
			else if( szPrevType == "o" )
				strOutFilenameBasis = szTemp;
			else if( szPrevType == "oformat" )
			{
				transform(szTemp.begin(), szTemp.end(), szTemp.begin(), ::tolower);
				if( szTemp == "gap" || szTemp == "latex" || szTemp == "mathematica" || szTemp == "pari" )
					strOuputMathematicalFormat = szTemp;
			}
			
			szPrevType = "";
		}
	}

	// --------------------------------------------------
	// Complete parameters
	string strParams;
	for( int i( 0 ); i < argc; ++i ) 
		strParams += " " + std::string( argv[i] );
	
	// --------------------------------------------------
	// Choosing the vertices
	PCRERegexp regexp;
	PCREResult regexpRes;
	int iRegexpCount;
	
	if( ( iRegexpCount = regexp.preg_match_all( "-S=\\{([[:alnum:][:space:],_-]+)\\}", strParams, regexpRes, PCRE_CASELESS ) ) > 0 )
	{
		str_replace( regexpRes[1][0], " ", "" );
		vector< string> strV( explode( ",", regexpRes[1][0] ) );
		unsigned int strVCount( strV.size() );
		
		for( unsigned int i(0); i < strVCount; i++ )
			strVertices.push_back( strV[i] );
	}
	else if( ( iRegexpCount = regexp.preg_match_all( "-S=([[:alnum:],_-]+)", strParams, regexpRes, PCRE_CASELESS ) ) > 0 )
	{
		for( int i(0); i < iRegexpCount; i++ )
			strVertices.push_back( regexpRes[1][i] );
	}
	
	// --------------------------------------------------
	// GBD
	if( ( iRegexpCount = regexp.preg_match_all( "-is=([[:alnum:],_-]+)", strParams, regexpRes, PCRE_CASELESS ) ) > 0 )
	{
		str_replace( regexpRes[1][0], " ", "" );
		vector< string> strV( explode( ",", regexpRes[1][0] ) );
		strGBDvertex_t0 = strV[0];
		if( strV.size() > 1 )
			strGBDvertex_s0 = strV[1];
			
		str_replace( strParams, regexpRes[0][0], "" );
		bGBD = true;
	}
	else if( ( iRegexpCount = regexp.preg_match_all( "-gbd ([[:alnum:]_-]+)", strParams, regexpRes, PCRE_CASELESS ) ) > 0 )
	{
		strGBDvertex_t0 = regexpRes[1][0];
		str_replace( strParams, regexpRes[0][0], "" );
		bGBD = true;
	}
	
	return true;
}

// TODO: si exposant > 50, afficher erreur pour growth rate (car facteurs pas forcément premiers entre eux)

void App::doGBD( CoxIter& ci )
{
	unsigned int iVerticesCountStart( ci.get_iVerticesCount() );
	
	if( strGBDvertex_t0 == "" )
	{
		cout << "Error: A vertex must be given" << endl;
		return;
	}
	
	GBD gbd( &ci );
	
	if( strGBDvertex_s0 != "" )
	{
		cout << "This is an experimental feature (to be properly tested)" << endl;
		cout << "------------------------------------------------------\n" << endl;
		
		cout << "Infinite sequence:" << endl;
		
		if( !gbd.bIsVertexAdmissible( strGBDvertex_t0 ) )
		{
			cout << "\tError: " << gbd.get_strError( ) << endl;
			return ;		
		}
		
		if( !gbd.bIsVertexAdmissible( strGBDvertex_s0 ) )
		{
			cout << "\tError: " << gbd.get_strError( ) << endl;
			return ;		
		}
		
		if( !ci.get_iDimension() )
		{
			cout << "\tError: The dimension must be specified" << endl;
			return ;		
		}
		
		if( ci.get_iCoxeterMatrixEntry( ci.get_iVertexIndex( strGBDvertex_t0 ), ci.get_iVertexIndex( strGBDvertex_s0 ) ) >= 2 )
		{
			cout << "\tError: The two hyperplanes must be (ultra)parallel " << endl;
			return ;
		}
		
		ci.IS_computations( strGBDvertex_t0, strGBDvertex_s0 );
		
		auto iFVUnits( ci.get_iISFVectorsUnits() ), iFVPowers( ci.get_iISFVectorsPowers() );
		
		cout << "\tf-vector after n doubling:\n\t(" << implode(", ", iFVUnits) << ", 1) + 2^(n-1)*(" << implode(", ", iFVPowers) << ", 0)" << endl;
	}
	else
		cout << "Index two subgroup:" << endl;
	
	// -----------------------------------------
	// Doing the GBD
	if( !gbd.removeVertex( strGBDvertex_t0 ) )
		cout << "\tError: " << gbd.get_strError( ) << endl;	
	else
		cout << "\tNumber of new hyperplanes after first doubling: " << ( ci.get_iVerticesCount() - iVerticesCountStart ) << endl;	
	
	cout << endl;
}

void App::run( )
{
	if( bPrintHelp )
	{
		printHelp();
		return;
	}
	
	chrono::time_point< std::chrono::system_clock > timeStart, timeEnd;
	
	bool bEulerSuccess( true ), bCanBeFiniteCovolume( false ), bSignatureComputed( false );
	
	
	CoxIter ci;
	ci.set_bCheckCocompactness( bCheckCocompacity );
	ci.set_bCheckCofiniteness( bCheckFiniteCovolume );
	ci.set_bDebug( bDebug );
	ci.set_bWriteInfo( true );
	ci.set_strOuputMathematicalFormat( strOuputMathematicalFormat );
	ci.set_strVerticesToConsider( strVertices );
	ci.set_strVerticesToRemove( strVerticesRemove );
	
	if( bCoutFile )
		ci.set_sdtoutToFile( strOutFilenameBasis + ".output" );
	
	Arithmeticity arithmeticity;
	
	#ifdef _COMPILE_WITH_PARI_
	GrowthRate_Result grr;
	grr.bComputed = false;
	grr.iPerron = -1;
	grr.iPisot = -1;
	grr.iSalem = -1;
	
	array< unsigned int, 3 > iSignature;
	#endif
	
	if( isatty( fileno(stdin) ) == 0 )
	{
		if( !ci.parseGraph( std::cin ) )
		{
			cout << "Error while reading graph: " << ci.get_strError( ) << endl;
			return;
		}	
	}
	else if( strInFilename != "" )
	{
		// Reading of the graph
		if( !ci.bReadGraphFromFile( strInFilename ) )
		{
			cout << "Error while reading file: " << ci.get_strError( ) << endl;
			return;
		}
	}
	else
	{
		cout << "No input file given\n" << endl;
		
		printHelp();
		return;
	}

	if( bGBD )
		doGBD( ci );
	
	if( bPrintGramMatrix )
		ci.printGramMatrix( );
	
	if( bPrintCoxeterMatrix )
		ci.printCoxeterMatrix( );
		
	if( bPrintCoxeterGraph )
		ci.printCoxeterGraph( );
	
	if( bOutputGraph && !ci.bWriteGraph( strOutFilenameBasis ) )
		cout << "Error while writing file: " << ci.get_strError( ) << endl;
		
	if( bOutputGraphToDraw )
	{
		if( ci.bWriteGraphToDraw( strOutFilenameBasis ) )
		{
			string strCommand( "dot -Tjpg -o\"" + strOutFilenameBasis + ".jpg\" \"" + strOutFilenameBasis + ".graphviz\"" );
			#ifdef _DOT_PROGRAM_FOUND_
				FILE *fin;
				if( ( fin = popen( strCommand.c_str(), "r" ) ) )
				{
					cout << "Image created: \n\t" <<  strOutFilenameBasis << ".jpg\n" << endl;
					pclose( fin );
				}
				else
					cout << "GraphViz command: \n\t" << strCommand << "\n" <<  endl;
			#else
				cout << "GraphViz command: \n\t" << strCommand << "\n" <<  endl;
			#endif
		}
		else
			cout << "Error while writing file: " << ci.get_strError( ) << endl;
	}
	
	if( bGBD )
		return;
	
	// -----------------------------------------------------------------
	// composantes connexes sphériques et euclidiennes
	cout << "Finding connected subgraphs......" << endl;
	timeStart = chrono::system_clock::now();
	ci.exploreGraph( );
	
	try{
		if( bCheckCanBeFiniteCovolume )
			bCanBeFiniteCovolume = ci.bCanBeFiniteCovolume( );
	}
	catch( const string& strE )
	{
		bCheckCanBeFiniteCovolume = false;
		cout << "\nError:\n\t" << strE << "\n" << endl;
	}
	
	// -----------------------------------------------------------------
	// calcul des produits (la majorité du temps de calcul concerne ce bloc)
	if( bComputeEuler || bComputeGrowthSeries || bCheckCocompacity || bCheckFiniteCovolume )
	{
		cout << "Finding graphs products......" << endl;
		ci.computeGraphsProducts( );
	}
	
	unsigned int iDimension( ci.get_iDimension( ) );
	
	if( bCheckFiniteCovolume )
		ci.isFiniteCovolume( );
	
	// -----------------------------------------------------------------
	// calcul de la caractéristique d'Euler, f-vecteur et compacité
	cout << "Computations......" << endl; 
	if( bComputeEuler && !ci.bEulerCharacteristicFVector( ) )
	{
		bEulerSuccess = false;
		cout << "\n\n######################################################################" << endl;
		cout << "\tAn error occurred." << endl;
		cout << "\tCheck the graph encoding." << endl;
		cout << "\tYou can run CoxIter with option '-debug' to see if the graph contains a spherical subgraph which has too big rank." << endl;
		cout << "######################################################################\n" << endl;
	}
	
	if( bCheckCocompacity )
		ci.iIsGraphCocompact( );
	
	if( bCheckArithmeticity )
		arithmeticity.test( ci, true );
	
	if( bComputeGrowthSeries )
	{
		ci.growthSeries();
	}
	
	if( ci.get_iHasDottedLineWithoutWeight() != 0 )
		bComputeSignature = false;
	
	if( bComputeSignature )
	{
		#ifdef _COMPILE_WITH_PARI_
		try{
			Signature s;
			iSignature = s.iComputeSignature( ci.get_strGramMatrix_PARI() );
			bSignatureComputed = true;
		}
		catch( const string& strE )
		{
			cout << "\n---------------------------------------------------------" << endl;
			cout << "Error while computing the signature:\n\t" << strE << endl;
			cout << "---------------------------------------------------------\n" << endl;
		}
		#endif
	}
	
	if( bComputeGrowthRate )
	{
		#ifdef _COMPILE_WITH_PARI_
		try{
			GrowthRate gr;
			grr = gr.grrComputations( ci.get_iGrowthSeries_denominator() );
		}
		catch( const string& strE )
		{
			cout << "\n---------------------------------------------------------" << endl;
			cout << "Error while computing the growth rate:\n\t" << strE << endl;
			cout << "---------------------------------------------------------\n" << endl;
		}
		#endif
	}
	
	timeEnd = chrono::system_clock::now();
	cout << "\tComputation time: " << chrono::duration <double, milli>(timeEnd-timeStart).count( ) / 1000 << "s\n" << endl;

	// -----------------------------------------------------------------
	// Affichage des informations
	cout << "Information" << endl;
	
	if( ci.get_bDimensionGuessed() )
		cout << "\tGuessed dimension: " << ci.get_iDimension() << endl;
	
	cout << "\tCocompact: " << ( ci.get_iIsCocompact( ) >= 0 ? ( ci.get_iIsCocompact( ) == 0 ? "no" : "yes" ) : "?" )  << endl;
	if( bCheckCanBeFiniteCovolume ) 
		cout << "\tCan be of finite covolume: " << ( bCanBeFiniteCovolume ? "yes" : "no" ) << endl;
	
	cout << "\tFinite covolume: " << ( ci.get_iIsFiniteCovolume( ) >= 0 ? ( ci.get_iIsFiniteCovolume( ) == 0 ? "no" : "yes" ) : "?" )  << endl;
	if( bCheckArithmeticity )
	{
		cout << "\tArithmetic: ";
		if( ci.get_iIsArithmetic( ) == 1 )
			cout << "yes" << endl;
		else if( ci.get_iIsArithmetic( ) == 0 )
			cout << "no" << endl;
		else
			cout << "?" << ( arithmeticity.get_strError( ) != "" ? "(" + arithmeticity.get_strError( ) + ")" : ( ci.get_bHasDottedLine( ) ? " (GRAPH HAS DOTTED EDGE)" : "" ) ) << endl;
	}
	
	// f-vector, alternating sum of the components of the f-vector
	if( bComputeEuler && bEulerSuccess )
	{
		if( iDimension )
		{
			vector< unsigned int > iFVector( ci.get_iFVector( ) );
			
			cout << "\tf-vector: (";
			for( unsigned int i( 0 ); i <= iDimension; i++ )
				cout << ( i ? ", " : "" ) << iFVector[i];
			cout << ")" << endl;
			
			cout << "\tNumber of vertices at infinity: " << ci.get_iVerticesAtInfinityCount( ) << endl;
			
			cout << "\tAlternating sum of the components of the f-vector: " << ci.get_iFVectorAlternateSum( ) << endl;
		}
		
		if( bDoComputations )
			cout << "\tEuler characteristic: " << ci.get_brEulerCaracteristic( ) << endl;
		else
			cout << "--------- BEGIN: Euler characteristic ---------\n" << ci.get_strEulerCharacteristic_computations( ) << "\n" << "--------- END: Euler characteristic ---------" << endl;
	}
	
	// volume
	if( bComputeEuler && iDimension && bEulerSuccess && !( iDimension % 2 ) && ci.get_iIsFiniteCovolume( ) == 1 )
	{
		cout << "\tCovolume: ";
		
		MPZ_rational cov( ( iDimension / 2 ) % 2 ? -1 : 1 );
		for( unsigned int i( 1 ); i <= iDimension; i++ )
		{
			cov *= 2;
			cov /= i;
			if( i <= ( iDimension / 2 ) )
				cov *= i;
		}
		
		cout << "pi^" << ( iDimension / 2 ) << " * " << cov * ci.get_brEulerCaracteristic( ) << endl;
	}
	
	if( bSignatureComputed )
	{
		#ifdef _COMPILE_WITH_PARI_
		cout << "\tSignature (numerically): " << iSignature[0] << "," << iSignature[1] << "," << iSignature[2] << endl;
		#endif
	}
	
	if( bComputeGrowthSeries )
	{
		cout << "\nGrowth series: " << endl;
		ci.printGrowthSeries();
		cout << endl;
		
		#ifdef _COMPILE_WITH_PARI_
		if( bComputeGrowthRate && grr.bComputed && ci.get_bGrowthSeriesReduced() )
		{
			cout << "\nGrowth rate: " << grr.strGrowthRate << endl;
			cout << "\tPerron number: " << ( grr.iPerron < 0 ? "?" : ( grr.iPerron > 0 ? "yes" : "no" ) ) << endl;
			cout << "\tPisot number: " << ( grr.iPisot < 0 ? "?" : ( grr.iPisot > 0 ? "yes" : "no" ) ) << endl;
			cout << "\tSalem number: " << ( grr.iSalem < 0 ? "?" : ( grr.iSalem > 0 ? "yes" : "no" ) ) << endl;
		}
		#endif
	}
	
	if( ci.get_iIsArithmetic() == -1 )
	{
		vector< string > strCycles( arithmeticity.get_strListCycles() );
		
		if( strCycles.size() )
		{
			cout << "\nThe group is arithmetic if and only if all the following values lie in Z: \n";
			auto strWeights( ci.get_strWeights() );
			
			cout << implode( "\n", strCycles ) << endl;
			cout << "with" << endl;
			for( auto it : strWeights )
				cout << "l" << iLinearizationMatrix_row( it.first, ci.get_iVerticesCount() ) << "m" << iLinearizationMatrix_col( it.first, ci.get_iVerticesCount() ) << " = " << it.second << endl;
		}
	}
	
	cout << endl;
}

void App::printHelp() const
{
	cout << "  _____          _____ _            \n"
" / ____|        |_   _| |\n"           
"| |     _____  __ | | | |_ ___ _ __ \n"
"| |    / _ \\ \\/ / | | | __/ _ \\ '__|\n"
"| |___| (_) >  < _| |_| ||  __/ |\n"   
" \\_____\\___/_/\\_\\_____|\\__\\___|_|   \n" << endl;

	cout << "CoxIter is a program to compute invariants of hyperbolic Coxeter groups\n\n"
	"The basis usage is as follows:\n\t./coxiter < file-describing-the-graph\n\n"
	"For example, you could do\n\t./coxiter < ../graphs/8-Bugaenko.coxiter\n\n"
	"You can also choose to compute/test only some of the invariants:\n"
	"\t-a  : test whether the group is arithmetic\n"
	"\t-c  : test whether the group is cocompact\n"
	"\t-fv : test whether the group has finite covolume\n"
	"\t-g  : growth series\n"
	"\t-gr : growth rate\n"
	"\t-s  : signature\n\n"
	"There is many more options regarding the format of the output and the possible computations.\nThe full documentation is available here:\n\thttp://coxiter.rgug.ch/" << endl;
}
