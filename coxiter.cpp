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

#include "coxiter.h"

CoxIter::CoxIter()
:	bCheckCocompactness( false ),
	bCheckCofiniteness( false ),
	bCoutFile( false ),
	bDebug( false ),
	bGramMatrixField( false ),
	bGrowthSeriesComputed( false ),
	bHasBoldLine( false ),
	bHasDottedLine( false ),
	iHasDottedLineWithoutWeight( 0 ),
	bWriteInfo( false ),
	bWriteProgress( false ),
	bGraphExplored( false ),
	bGraphsProductsComputed( false ),
	bUseOpenMP( true ),
	brEulerCaracteristic( 0 ),
	graphsList_spherical( 0 ),
	graphsList_euclidean( 0 ),
	iDimension( 0 ),
	iEuclideanMaxRankFound( 0 ),
	iSphericalMaxRankFound( 0 ),
	bDimension_guessed( false ),
	iFVectorAlternateSum( 0 ),
	iIsArithmetic( -1 ),
	iIsCocompact( -2 ),
	iIsFiniteCovolume( -2 ),
	iVerticesAtInfinityCount( 0 ),
	iVerticesCount( 0 ),
	outCout( 0 ),
	sBufOld( 0 ),
	strError( "" ),
	strOuputMathematicalFormat( "generic" )
{ 
	#ifndef _COMPILE_WITH_OPENMP_
	this->bUseOpenMP = false;
	#endif
}

CoxIter::CoxIter( const vector< vector< unsigned int > >& iMatrix, const unsigned int& iDimension )
:	bCheckCocompactness( false ),
	bCheckCofiniteness( false ),
	bCoutFile( false ),
	bGramMatrixField( false ),
	bGraphExplored( false ),
	bGraphsProductsComputed( false ),
	bGrowthSeriesComputed( false ),
	bHasDottedLine( false ),
	iHasDottedLineWithoutWeight( 0 ),
	bWriteInfo( false ), 
	bWriteProgress( false ),
	bDebug( false ),
	bUseOpenMP( true ),
	brEulerCaracteristic( 0 ),
	graphsList_spherical( 0 ),
	graphsList_euclidean( 0 ),
	iDimension( iDimension ),
	iEuclideanMaxRankFound( 0 ),
	iSphericalMaxRankFound( 0 ),
	bDimension_guessed( false ),
	iFVectorAlternateSum( 0 ),
	iIsCocompact( -1 ),
	iIsFiniteCovolume( -1 ),
	iVerticesAtInfinityCount( 0 ),
	iVerticesCount( 0 ),
	outCout( 0 ),
	sBufOld( 0 ),
	strError( "" ),
	strOuputMathematicalFormat( "" )
{
	iVerticesCount = iMatrix.size();
	
	initializations();
	
	iCoxeterMatrix = iMatrix;
	
	#ifndef _COMPILE_WITH_OPENMP_
	this->bUseOpenMP = false;
	#endif
}

CoxIter::~CoxIter()
{
	if( graphsList_spherical )
		delete graphsList_spherical;
	
	if( graphsList_euclidean )
		delete graphsList_euclidean;
	
	// if cout is redirected to a file
	if( bCoutFile )
	{
		outCout->close();
		cout.rdbuf( sBufOld ); // we restore the cout
	}
}

bool CoxIter::bRunAllComputations()
{
	if( !iCoxeterMatrix.size() )
		return false;
	
	if( !bGraphExplored )
		exploreGraph();
	
	if( !bGraphsProductsComputed )
		computeGraphsProducts();
	
	if( !bEulerCharacteristicFVector() )
		return false;
	
	if( bCheckCofiniteness )
		isFiniteCovolume();
	
	if( bCheckCocompactness )
		iIsGraphCocompact();
	
	return true;
}

#ifndef _COMPILE_WITHOUT_REGEXP_
bool CoxIter::parseGraph( istream& streamIn )
{
	string strLine;
	PCRERegexp regexp;
	PCREResult regexpRes;
	
	// loops variable, first vertice, second vertice, weight, number of vertices, index of the current row
	unsigned int i, i1, i2, i3, iVerticesFileCount, iRowIndex(1);
	
	vector< unsigned int >::const_iterator it;
	
	vector< unsigned int > iOrders; // orders found
	
	// ---------------------------------------------------------------------------
	// Reading the number of vertices and, eventually, dimension
	if( getline( streamIn, strLine ) )
	{
		if( regexp.preg_match_all( "([[:digit:]]+)[[:space:]]?([[:digit:]]*)", strLine, regexpRes ) == 1 )
		{
			iVerticesFileCount = iVerticesCount = stoi( regexpRes[1][0] );
			iDimension = regexpRes[2][0] != "" ? stoi( regexpRes[2][0] ) : 0;
		}
		else
		{
			strError = "First line with number of vertices missing";
			return false;
		}
	}
	else
	{
		strError = "EMPTY_FILE";
		return false;
	}
	
	// ---------------------------------------------------------------------------
	// first line
	if( !getline( streamIn, strLine ) )
	{
		strError = "EMPTY_FILE";
		return false;
	}
	
	// names of the vertices
	if( regexp.preg_match_all( "^vertices labels:[[:space:]]?([[:alnum:]-_ ]+)$", strLine, regexpRes ) )
	{
		vector< string > strVL( explode( " ", regexpRes[1][0] ) );
		if( strVL.size() != iVerticesFileCount )
		{
			strError = "VERTICES_LABEL_COUNT";
			return false;
		}
		
		for( i = 0; i < iVerticesFileCount; i++ )
		{
			map_vertices_labelToIndex[ strVL[i] ] = i;
			map_vertices_indexToLabel.push_back( strVL[i] );
		}
		
		if( map_vertices_labelToIndex.size() != iVerticesFileCount )
		{
			strError = "VERTICES_LABEL_COUNT";
			return false;
		}
		
		if( !getline( streamIn, strLine ) )
		{
			strError = "EMPTY_FILE";
			return false;
		}
	}
	else
	{
		for( i = 0; i < iVerticesFileCount; i++ )
		{
			map_vertices_labelToIndex[ to_string( i + 1 ) ] = i;
			map_vertices_indexToLabel.push_back( to_string( i + 1 ) );
		}
	}
	
	bool bRemoveDottedEdges( false ); // If we want to remove dotted edges
	
	// ---------------------------------------------------------------------------
	// removed vertices
	vector< unsigned int > iVerticesShift( iVerticesFileCount, 0 ) ; // Shifts for the removed vertices
	unsigned int iTruncCount( 0 );
	
	if( strVertices.size() ) // If we want to specify a subset of the vertices
	{
		auto strAllVertices( map_vertices_indexToLabel );
		sort( strAllVertices.begin(), strAllVertices.end() );
		
		set_difference( strAllVertices.begin(), strAllVertices.end(), strVertices.begin(), strVertices.end(), std::back_inserter(strVerticesRemove) );
		sort( strVerticesRemove.begin(), strVerticesRemove.end() );
	}
	
	for( vector< string >::const_iterator itStr( strVerticesRemove.begin() ); itStr != strVerticesRemove.end(); ++itStr )
	{
		if( *itStr == "dotted" && map_vertices_labelToIndex.find( "dotted" ) == map_vertices_labelToIndex.end() ) // remove dotted edges?
		{
			bRemoveDottedEdges = true;
			continue;
		}
		
		if( map_vertices_labelToIndex.find( *itStr ) == map_vertices_labelToIndex.end() )
		{
			strError = "This vertex does not exist: " + *itStr;
			return false;
		}
		
		iTruncCount++;
		
		for( i = map_vertices_labelToIndex[ *itStr ]; i < iVerticesFileCount; i++ )
			iVerticesShift[ i ]++;
	}
	iVerticesCount -= iTruncCount;
	
	// ---------------------------------------------------------------------------
	// initializations
	initializations(); // now that we know the real number of vertices

	// ---------------------------------------------------------------------------
	// reading the graph
	do
	{	
		regexpRes.clear();
		
		// Usual row: "first vertice" "second vertice" "weight"
		if( regexp.preg_match_all( "([[:alnum:]_-]+)[[:space:]]([[:alnum:]_-]+)[[:space:]]([[:digit:]]+)([[:space:]]+#[[:space:]]*([^\n]+))?", strLine, regexpRes ) )
		{
			if( map_vertices_labelToIndex.find( regexpRes[1][0] ) == map_vertices_labelToIndex.end() )
			{
				strError = "The following vertex is unknown: " + regexpRes[1][0];
				return false;
			}
			
			if( map_vertices_labelToIndex.find( regexpRes[2][0] ) == map_vertices_labelToIndex.end() )
			{
				strError = "The following vertex is unknown: " + regexpRes[2][0];
				return false;
			}
			
			i1 = map_vertices_labelToIndex[ regexpRes[1][0] ];
			i2 = map_vertices_labelToIndex[ regexpRes[2][0] ];
			i3 = stoi( regexpRes[3][0] );
			
			if( i3 == 1 && bRemoveDottedEdges )
				i3 = 2;
			
			iRowIndex ++;

			// Removed vertex?
			if( binary_search( strVerticesRemove.begin(), strVerticesRemove.end(), regexpRes[1][0] ) || binary_search( strVerticesRemove.begin(), strVerticesRemove.end(), regexpRes[2][0] ) )
				continue;
			
			// on tient compte du décalage lié à la troncation
			i1 -= iVerticesShift[ i1 ];
			i2 -= iVerticesShift[ i2 ];
			
			// on garde le poids (pour corps engendré par les coefficients de la matrice de Gram)
			if( find( iOrders.begin(), iOrders.end(), i3 ) == iOrders.end() )
				iOrders.push_back( i3 );
			
			if( i3 == 1 ) // Weight of the dotted line given?
			{
				if( regexpRes.size() > 5 )
				{
					unsigned int iIndex( iLinearizationMatrix_index( min(i1,i2), max(i1,i2), iVerticesCount ) );
					strWeights[ iIndex ] = regexpRes[5][0];
				}
				else
					iHasDottedLineWithoutWeight = 1;
			}
			
			// si on avait déjà cette arête avec un ordre différent
			if( iCoxeterMatrix[i1][i2] != 2 && iCoxeterMatrix[i1][i2] != i3 )
			{
				strError = "Edge has multiple orders (" + regexpRes[1][0] + "," + regexpRes[2][0] + ")";
				return false;
			}
			
			iCoxeterMatrix[i1][i2] = i3;
			iCoxeterMatrix[i2][i1] = i3;
			
			if( i3 == 1 ) // dotted
				bHasDottedLine = true;
			else if( i3 == 0 )
				bHasBoldLine = true;
		}
		else if( strLine != "" )
		{
			if( bWriteInfo )
				cout << "Unread line (incorrect format): " << "#" << strLine << "#" << iRowIndex << endl;
			
			iRowIndex ++;
			continue;
		}
	}while( getline( streamIn, strLine ) );
	
	// ---------------------------------------------------------------------------
	// Labels and co
	auto v_ItL( map_vertices_indexToLabel );
	map_vertices_indexToLabel.clear();
	map_vertices_labelToIndex.clear();
	
	unsigned int j(0);
	for( unsigned int i(0); i < iVerticesFileCount; i++ )
	{
		if( ( !i && !iVerticesShift[i] ) || ( i && iVerticesShift[i] == iVerticesShift[i-1] ) )
		{
			map_vertices_indexToLabel.push_back( v_ItL[i]  );
			map_vertices_labelToIndex[ v_ItL[i] ] = j++;
		}
	}
	
	// ---------------------------------------------------------------------------
	// some information
	if( bWriteInfo )
	{
		cout << "Reading graph: " << endl;
		cout << "\tNumber of vertices: " << iVerticesCount << endl;
		cout << "\tDimension: " << ( iDimension ? to_string( iDimension ) : "?" ) << endl;
		
		cout << "\tVertices: ";
		for( vector< string >::const_iterator itStr( map_vertices_indexToLabel.begin() ); itStr != map_vertices_indexToLabel.end(); ++itStr )
			cout << ( itStr != map_vertices_indexToLabel.begin() ? ", " : "" ) << *itStr;
		cout << endl;
	}
	
	// ---------------------------------------------------------------------------
	// Field generated by the entries of the Gram matrix
	for( it = iOrders.begin(); it != iOrders.end(); ++it )
	{
		if( *it == 1 ) // dotted line
			break;
		else if( *it == 4 )
			strGramMatrixField += ( strGramMatrixField == "" ? "sqrt(2)" : ", sqrt(2)" );
		else if( *it == 5 )
			strGramMatrixField += ( strGramMatrixField == "" ? "sqrt(5)" : ", sqrt(5)" );
		else if( *it == 6 )
			strGramMatrixField += ( strGramMatrixField == "" ? "sqrt(3)" : ", sqrt(3)" );
		else if( *it >= 7 ) // le static_cast est là pour VC++
			strGramMatrixField += ( strGramMatrixField == "" ? "cos(pi/" + to_string( static_cast<long long>(*it) ) + ")" : ", cos(pi/" + to_string( static_cast<long long>(*it) ) + ")" );
	}
	
	if( it == iOrders.end() )
	{
		strGramMatrixField = strGramMatrixField == "" ? "Q" : ( "Q[" + strGramMatrixField + "]" );
		bGramMatrixField = true;
		
		if( bWriteInfo )
			cout << "\tField generated by the entries of the Gram matrix: " << strGramMatrixField << endl;
	}
	else
	{
		if( bWriteInfo )
			cout << "\tField generated by the entries of the Gram matrix: ?" << endl;
	}
	
	if( bWriteInfo )
		cout << "File read\n" << endl;
	
	return true;
}


bool CoxIter::bReadGraphFromFile( const string& strInputFilename )
{
	// ---------------------------------------------------------------------------
	// try to open the file
	ifstream fileIn( strInputFilename.c_str() );
	if( fileIn.fail() )
	{
		strError = "Cannot open file";
		return false;
	}
	
	if( !parseGraph( fileIn ) )
		return false;
	
	fileIn.close();
	
	return true;
}
#endif

void CoxIter::initializations()
{
	size_t i;
	
	// ------------------------------------------------------
	// de initializations
	if( graphsList_spherical )
		delete graphsList_spherical;
	
	if( graphsList_euclidean )
		delete graphsList_euclidean;
	
	graphsProductsCount_spherical.clear();
	graphsProductsCount_euclidean.clear();
	
	iFactorials.clear();
	iPowersOf2.clear();
	
	bGraphExplored = false;
	bGraphsProductsComputed = false;
	
	// ------------------------------------------------------
	// initializations
	iCoxeterMatrix = vector< vector<unsigned int> >(iVerticesCount, vector<unsigned int>( iVerticesCount, 2 ) );
	bVerticesVisited = vector< bool >( iVerticesCount, false );
	bEdgesVisited = vector< vector<bool> >(iVerticesCount, vector<bool>( iVerticesCount, false ) );
	
	graphsList_spherical = new GraphsList( iVerticesCount, &map_vertices_indexToLabel );
	graphsList_euclidean = new GraphsList( iVerticesCount, &map_vertices_indexToLabel );
	
	graphsProductsCount_euclidean = vector< map<vector< vector< short unsigned int > >, unsigned int> >( iVerticesCount + 1, map<vector< vector< short unsigned int > >, unsigned int>() );
	graphsProductsCount_spherical = vector< map<vector< vector< short unsigned int > >, unsigned int> >( iVerticesCount + 1, map<vector< vector< short unsigned int > >, unsigned int>() );
	
	// ------------------------------------------------------------
	// sauvegarde de quelques calculs
	iFactorials = vector< mpz_class > ( iVerticesCount + 2, 1 );
	iPowersOf2 = vector< mpz_class > ( iVerticesCount + 2, 1 );
	for( i = 1; i <= iVerticesCount + 1; i++ )
	{
		iFactorials[i] = iFactorials[i-1] * i;
		iPowersOf2[i] = mpz_class(2) * iPowersOf2[i-1];
	}
}

bool CoxIter::bWriteGraph( const string& strOutFilenameBasis )
{
	if( strOutFilenameBasis == "" )
	{
		strError = "No file specified for writing the graph";
		return false;
	}
	
	map_vertices_labels_create();
	
	string strFilename( strOutFilenameBasis + ".coxiter" );
	ofstream out( strFilename.c_str() );
	if( !out.is_open() )
	{
		strError = "Cannot open the file for writing the graph";
		return false;
	}
	
	out << iVerticesCount << ( iDimension ? " " + to_string( iDimension ) : "" ) << endl;
	out << "vertices labels: ";
	for( vector< string >::const_iterator it( map_vertices_indexToLabel.begin() ); it != map_vertices_indexToLabel.end(); ++it )
		out << ( it == map_vertices_indexToLabel.begin() ? "" : " " ) << *it;
	out << endl;
	
	for( unsigned int i(0); i < iVerticesCount; i++ )
	{
		for( unsigned int j( 0 ); j < i; j++ )
		{
			if( iCoxeterMatrix[i][j] != 2 )
				out << map_vertices_indexToLabel[j] << " " << map_vertices_indexToLabel[i] << " " << iCoxeterMatrix[i][j] << endl;
		}
	}
	
	out.close();
	
	return true;
}

void CoxIter::map_vertices_labels_create()
{
	if( map_vertices_indexToLabel.size() )
		return; // nothing to do
		
	for( unsigned int i( 0 ); i < iVerticesCount; i++ )
	{
		map_vertices_labelToIndex[ to_string( i + 1 ) ] = i;
		map_vertices_indexToLabel.push_back( to_string( i + 1 ) );
	}
}

void CoxIter::map_vertices_labels_reinitialize()
{
	map_vertices_labelToIndex.clear();
	map_vertices_indexToLabel.clear();
	
	for( unsigned int i( 0 ); i < iVerticesCount; i++ )
	{
		map_vertices_labelToIndex[ to_string( i + 1 ) ] = i;
		map_vertices_indexToLabel.push_back( to_string( i + 1 ) );
	}
}


bool CoxIter::bWriteGraphToDraw( const string& strOutFilenameBasis )
{
	unsigned int i, j;
	
	map_vertices_labels_create();
	
	// ----------------------------------------------------------------------
	// ouverture du fichier
	if( strOutFilenameBasis == "" )
	{
		strError = "No file specified for writing the graph";
		return false;
	}
	
	string strFilename( strOutFilenameBasis + ".graphviz" );
	ofstream out( strFilename.c_str() );
	if( !out.is_open() )
	{
		strError = "Cannot open the file for writing the graph";
		return false;
	}
	
	// ----------------------------------------------------------------------
	// écriture à proprement parler
	out << "graph G { " << endl;
	for( i = 0; i < iVerticesCount; i++ )
		out << "\t\"" << map_vertices_indexToLabel[i] << "\";" << endl;
	
	for( i = 0; i < iVerticesCount; i++ )
	{
		for( j = i + 1; j < iVerticesCount; j++ ) 
		{
			if( iCoxeterMatrix[i][j] > 3 )
				out << "\t\"" << map_vertices_indexToLabel[i] << "\" -- \"" << map_vertices_indexToLabel[j] << "\" [label=\"" << iCoxeterMatrix[i][j] << "\"];" << endl;
			else if( iCoxeterMatrix[i][j] > 2 )
				out << "\t\"" << map_vertices_indexToLabel[i] << "\" -- \"" << map_vertices_indexToLabel[j] << "\";" << endl;
			else if( iCoxeterMatrix[i][j] == 1 )
				out << "\t\"" << map_vertices_indexToLabel[i] << "\" -- \"" << map_vertices_indexToLabel[j] << "\" [style=dotted];" << endl;
			else if( iCoxeterMatrix[i][j] == 0 )
				out << "\t\"" << map_vertices_indexToLabel[i] << "\" -- \"" << map_vertices_indexToLabel[j] << "\" [label=\"inf\"];" << endl;
		}
	}
	
	out << "}";
	out.close();
	
	// indique la commande GraphViz (l'image générée aide à vérifier que le graphe a été correctement encodé)
	if( bWriteInfo )
		cout << "GraphViz command: \n\t" <<  "dot -Tjpg -o\"" << strOutFilenameBasis << ".jpg\" \"" << strFilename << "\"\n" <<  endl;
	
	return true;
}

void CoxIter::exploreGraph()
{
	vector<short unsigned int> iVertices;
	short unsigned int i, j, k, l;
	
	if( !iVerticesCount )
		throw( string( "CoxIter::exploreGraph: No graph given" ) );
	
	if( bGraphExplored )
		return;
	
	// -------------------------------------------------------------------
	// pour chaque sommet, on cherche toutes les chaînes qui partent, ce qui donne les An, Bn, Dn, En, Hn
	iPath.clear();
	
	if( bWriteProgress )
		cout << "\tFirst part of spherical graphs" << endl;
	
	for( i = 0; i < iVerticesCount; i++ )
	{
		iPath.clear();
		bVerticesVisited = vector< bool >( iVerticesCount, false );
		bEdgesVisited = vector< vector<bool> >(iVerticesCount, vector<bool>( iVerticesCount, false ) );
		
		DFS( i, i );
		
		if( bWriteProgress && i && !( i % 10 ) )
			cout << "\t\t" << ( 100 * i / iVerticesCount ) << "%" << endl;
	}
	if( bWriteProgress )
		cout << "\n\tFinished" << endl;
	
	// -------------------------------------------------------------------
	// recherche des A_1, G_2^k avec k >= 4, F_4
	if( bWriteProgress )
		cout << "\tSecond part of spherical graphs" << endl;
	
	vector<bool> bVerticesLinkable, bVerticesLinkableTemp;
	for( i = 0; i < iVerticesCount; i++ )
	{
		bVerticesLinkable = vector<bool>( iVerticesCount, true );
		for( j = 0; j < iVerticesCount; j++ )
		{
			if( iCoxeterMatrix[i][j] != 2 )
				bVerticesLinkable[j] = false;
		}
		
		// ajout du sommet (A_1)
		bVerticesLinkable[i] = false;
		graphsList_spherical->addGraph( vector<short unsigned int>(1, i), bVerticesLinkable, 0, true );
		
		// on regarde si on trouve avec ce sommet: Gn, TA1 ,TC2
		for( j = 0; j < iVerticesCount; j++ )
		{
			if( iCoxeterMatrix[i][j] >= 4 || !iCoxeterMatrix[i][j] )
			{
				// ------------------------------------------------------------------
				// G2 et TA1
				if( i < j )
				{
					iVertices.clear();
					iVertices.push_back( i );
					iVertices.push_back( j );
					
					bVerticesLinkableTemp = bVerticesLinkable;
					for( k = 0; k < iVerticesCount; k++ )
					{
						if( iCoxeterMatrix[j][k] != 2 )
							bVerticesLinkableTemp[k] = false;
					}
			
					if( iCoxeterMatrix[i][j] ) // ici, c'est un graphe sphérique
						graphsList_spherical->addGraph( iVertices, bVerticesLinkableTemp, 6, true, 0, 0, iCoxeterMatrix[i][j] );
					else // ici, graphe euclidien (TA1)
						graphsList_euclidean->addGraph( iVertices, bVerticesLinkableTemp, 0, false, 0, 0, 0 );
				}
				
				// ------------------------------------------------------------------
				// TC2 = [ 4, 4 ]
				if( iCoxeterMatrix[i][j] == 4 )
				{
					for( k = 0; k < iVerticesCount; k++ )
					{	
						if( iCoxeterMatrix[k][j] == 4 && i != k && iCoxeterMatrix[i][k] == 2 )
						{
							bVerticesLinkableTemp = bVerticesLinkable;
							for( l = 0; l < iVerticesCount; l++ )
							{
								if( iCoxeterMatrix[k][l] != 2 )
									bVerticesLinkableTemp[l] = false;
								if( iCoxeterMatrix[j][l] != 2 )
									bVerticesLinkableTemp[l] = false;
							}
							
							graphsList_euclidean->addGraph( vector<short unsigned int>( 1, j ), bVerticesLinkableTemp, 2, false, i, k );
						}
					}
				}
			}
		}
		
		if( bWriteProgress && i && !( i % 10 ) )
			cout << "\t\t" << ( 100 * i / iVerticesCount ) << "%" << endl;
	}
	if( bWriteProgress )
		cout << "\n\tFinished" << endl;
	
	bGraphExplored = true;
}

void CoxIter::DFS( unsigned int iRoot, unsigned int iFrom )
{
	// -------------------------------------------------------------------
	// initializations
	bool bSubcall( false ); // to know if we call DFS
	
	unsigned int i;
	
	/*
	 * 	We don't want cycles
	 * 		We mark neighbours of iFrom as visited (to avoir cycles)
	 * 		We stock this in verticesVisited to restore it at the end
	 */
	vector< unsigned int > visitedVertices;
	
	if( iRoot != iFrom )
	{
		for( i = 0; i < iVerticesCount; i++ )
		{
			if( iCoxeterMatrix[iFrom][i] != 2 )
			{
				if( !bVerticesVisited[i] )
					visitedVertices.push_back( i );
				
				bVerticesVisited[i] = true;
			}
		}
	}
	else
		bVerticesVisited[iRoot] = true; // obviously...
	
	// -------------------------------------------------------------------
	// DFS
	iPath.push_back( iRoot ); // we add iRoot to the path
	
	for( i = 0; i < iVerticesCount; i++ )
	{
		// if we have an edge AND if i was not traversed AND if the edge (iRoot,i) was not traversed
		if( iCoxeterMatrix[iRoot][i] == 3 && !bVerticesVisited[i] && !bEdgesVisited[iRoot][i] )
		{
			bEdgesVisited[iRoot][i] = bEdgesVisited[i][iRoot] = true;
			bSubcall = true;
			DFS( i, iRoot );
		}
	}
	
	// -------------------------------------------------------------------
	// un-initializations
	for( vector<unsigned int>::iterator it( visitedVertices.begin() ); it != visitedVertices.end(); ++it )
		bVerticesVisited[ *it ] = false;
	
	bVerticesVisited[iRoot] = false;
	
	if( iFrom != iRoot ) // If it is a recursive call
		bEdgesVisited[iRoot][iFrom] = bEdgesVisited[iFrom][iRoot] = false;
	
	// If DFS was not called, then the path is maximal
	if( !bSubcall )
		addGraphsFromPath();
	
	iPath.pop_back();
}

void CoxIter::addGraphsFromPath()
{
	// sommets que l'on ne peut pas lier au graphe (n sommets);
	vector<bool> bVerticesLinkable( iVerticesCount, true );
	
	// sommets que l'on ne peut pas lier au graphe (n-1 sommets), sommets que l'on ne peut pas lier au graphe (n-2 sommets)
	vector<bool> bVerticesLinkable_0_nMin1( iVerticesCount, true ), bVerticesLinkable_0_nMin2( iVerticesCount, true );
	
	// sommets que l'on ne peut pas lier au graphe (1 --> n), sommets que l'on ne peut pas lier au graphe (1 --> n-1), , sommets que l'on ne peut pas lier au graphe (2 --> n)
	vector<bool> bVerticesLinkable_1_n( iVerticesCount, true ), bVerticesLinkable_1_nMin1( iVerticesCount, true ), bVerticesLinkable_2_n( iVerticesCount, true );
	
	// vecteur temporaire
	vector<bool> bVerticesLinkableTemp, bVerticesLinkableTempTemp;
	
	// chemin en cours de construction
	vector<short unsigned int> iPathTemp;
	
	// i, j, k, l: variables de boucles
	short unsigned int i, j, k, l, iMax( iPath.size() ), iOrder;
	
	for( i = 0; i < iMax; i++ )
	{
		iPathTemp.push_back( iPath[i] );
		
		// --------------------------------------------------------------------
		// mise à jour des voisinages occupés
		for( j = 0; j < iVerticesCount; j++ )
		{
			if( iCoxeterMatrix[ iPath[i] ][ j ] != 2 )
				bVerticesLinkable[ j ] = false;
			
			if( i >= 1 && iCoxeterMatrix[ iPath[i] ][ j ] != 2 )
				bVerticesLinkable_1_n[ j ] = false;
			
			if( i >= 2 && iCoxeterMatrix[ iPath[i] ][ j ] != 2 )
				bVerticesLinkable_2_n[ j ] = false;
			
			if( i >= 1 && iCoxeterMatrix[ iPath[i - 1] ][ j ] != 2 )
				bVerticesLinkable_0_nMin1[ j ] = false;
			
			if( i >= 2 && iCoxeterMatrix[ iPath[i - 1] ][ j ] != 2 )
				bVerticesLinkable_1_nMin1[ j ] = false;
			
			if( i >= 2 && iCoxeterMatrix[ iPath[i - 2] ][ j ] != 2 )
				bVerticesLinkable_0_nMin2[ j ] = false;
		}
		
		// --------------------------------------------------------------------
		// An
		if( i != 0 ) // on ajoute pas les sommets
			graphsList_spherical->addGraph( iPathTemp, bVerticesLinkable, 0, true );
		
		// --------------------------------------------------------------------
		// TAn, n >= 2
		if( i >= 1 )
		{
			bVerticesLinkable_1_nMin1[ iPathTemp[ i - 1 ] ] = false;
			
			for( j = 0; j < iVerticesCount; j++ )
			{
				// si la partie centrale ne pose pas de problème ET qu'on est lié à chaque extrémité
				if( bVerticesLinkable_1_nMin1[j] && !bVerticesLinkable[j] && iCoxeterMatrix[j][ iPathTemp[0] ] == 3 &&  iCoxeterMatrix[j][ iPathTemp[i] ] == 3 )
				{
					// mise à jour des linkables avec le somme trouvé
					bVerticesLinkableTemp = bVerticesLinkable;
					for( k = 0; k < iVerticesCount; k++ )
					{
						if( iCoxeterMatrix[k][j] != 2 )
							bVerticesLinkableTemp[k] = false;
					}
					
					graphsList_euclidean->addGraph( iPathTemp, bVerticesLinkableTemp, 0, false, j, 0, 1 ); // TODO OPTIMIZATION modifier ce 1 (relatif à une meilleure valeur que "0" par défaut pour les dernières variables)
				}
			}
		}
		
		// --------------------------------------------------------------------
		// Dn, TBn, TDn
		if( i >= 2 )
		{
			bVerticesLinkable_0_nMin2[ iPathTemp[ i - 2 ] ] = false;

			// --------------------------------------------------------------------
			// Dn, TBn (n >= 4), TDn (n >= 4)
			// on regarde les voisins de l'avant dernier sommet
			for( j = 0; j < iVerticesCount; j++ )
			{
				// si y'a une arête ET si on est pas déjà dans le chemin ET si pas interdit ET pas lien entre deux extrémités
				if( iCoxeterMatrix[j][ iPathTemp[ i-1 ] ] == 3 && ( iPathTemp[i-1] != j && iPathTemp[i] != j ) && bVerticesLinkable_0_nMin2[j] && iCoxeterMatrix[j][ iPathTemp[i] ] == 2 )
				{
					bVerticesLinkableTemp = bVerticesLinkable;
					for( k = 0; k < iVerticesCount; k++ )
					{
						if( iCoxeterMatrix[k][j] != 2 )
							bVerticesLinkableTemp[k] = false;
					}
					
					graphsList_spherical->addGraph( iPathTemp, bVerticesLinkableTemp, 3, true, j ); // Dn
					
					// --------------------------------------------------------------------
					// ici, on va tenter de trouver un TD_n (n >= 4) (i.e. prolonger par une arrête au 2ème sommet)
					for( k = 0; k < iVerticesCount; k++ )
					{
						if( k != iPathTemp[0] && k != j && k != iPathTemp[i] && iCoxeterMatrix[ iPathTemp[1] ][k] == 3 && bVerticesLinkable_2_n[k] && iCoxeterMatrix[ iPathTemp[0] ][k] == 2 && iCoxeterMatrix[k][j] == 2 )
						{
							bVerticesLinkableTempTemp = bVerticesLinkableTemp;
							for( l = 0; l < iVerticesCount; l++ )
							{
								if( iCoxeterMatrix[k][l] != 2 )
									bVerticesLinkableTempTemp[l] = false;
							} 
							
							graphsList_euclidean->addGraph( iPathTemp, bVerticesLinkableTempTemp, 3, false, k, j ); // TDn
						}
					}
					
					// --------------------------------------------------------------------
					// ici, on va tenter de trouver un TB_n (n >= 4) (i.e. prolonger par une arête de poids 4 à gauche)
					for( k = 0; k < iVerticesCount; k++ )
					{
						if( iCoxeterMatrix[ iPathTemp[0] ][k] == 4 && bVerticesLinkable_1_n[k] && iCoxeterMatrix[j][k] == 2 )
						{
							bVerticesLinkableTempTemp = bVerticesLinkableTemp;
							for( l = 0; l < iVerticesCount; l++ )
							{
								if( iCoxeterMatrix[k][l] != 2 )
									bVerticesLinkableTempTemp[l] = false;
							}
							
							graphsList_euclidean->addGraph( iPathTemp, bVerticesLinkableTempTemp, 1, false, j, k, 1 ); // TBn
						}
					}
				}

				// --------------------------------------------------------------------
				// TB3
				if( i == 2 )
				{
					if( iCoxeterMatrix[iPathTemp[1]][j] == 4 && iCoxeterMatrix[iPathTemp[0]][j] == 2 && iCoxeterMatrix[iPathTemp[2]][j] == 2 )
					{
						bVerticesLinkableTemp = bVerticesLinkable;
						for( k = 0; k < iVerticesCount; k++ )
						{
							if( iCoxeterMatrix[k][j] != 2 )
								bVerticesLinkableTemp[k] = false;
						}
						
						graphsList_euclidean->addGraph( iPathTemp, bVerticesLinkableTemp, 1, false, j ); // TB3
					}
				}
			}
		}
		
		// --------------------------------------------------------------------
		// E6, E7, E8, TE6, TE7, TE8
		if( i >= 4 && i <= 7 )
		{
			AnToEn_AnToTEn( iPathTemp, bVerticesLinkable );
		}
		
		// --------------------------------------------------------------------
		// Bn, F4 et Hn, TG2, TCn et TF4
		if( i >= 1 )
		{
			for( j = 0; j < iVerticesCount; j++ ) // on regarde si on peut prolonger la chaîne de 1 avec une arête de poids 4
			{
				/*
				 * 	Le premier paquet de conditions donne:
				 * 		Prolonger par une arrête de poids 4 --> Cn ; sphérique
				 * 		Prolonger par une arrête de poids 5 si (2 ou 3 sommets) --> (H3 ou H4) ; sphérique
				 * 		Prolonger par une arrête de poids 6 --> [ 3, 6 ] ; euclidien
				*/
				if( ( ( iCoxeterMatrix[ iPath[i] ][j] == 4 ) || ( iCoxeterMatrix[ iPath[i] ][j] == 5 && ( i == 1 || i == 2 ) ) || ( iCoxeterMatrix[ iPath[i] ][j] == 6 && i == 1 ) ) && bVerticesLinkable_0_nMin1[j] )
				{
					iOrder = iCoxeterMatrix[ iPath[i] ][j];
					bVerticesLinkableTemp = bVerticesLinkable;
					for( k = 0; k < iVerticesCount; k++ )
					{
						if( iCoxeterMatrix[k][j] != 2 )
							bVerticesLinkableTemp[k] = false;
					}
					
					if( iCoxeterMatrix[ iPath[i] ][j] < 6 ) // sphérique
						graphsList_spherical->addGraph( iPathTemp, bVerticesLinkableTemp, ( iOrder == 4 ? 1 : 7 ), true, j );
					else
						graphsList_euclidean->addGraph( iPathTemp, bVerticesLinkableTemp, 6, false, j, 0, 1 );
					
					// ------------------------------------------
					// on va tenter de prolonger cela en un TCn, n \geq 3
					if( iCoxeterMatrix[ iPath[i] ][j] == 4 )
					{
						for( k = 0; k < iVerticesCount; k++ )
						{
							if( iCoxeterMatrix[k][ iPathTemp[0] ] == 4 && k != j && bVerticesLinkable_1_n[k] && iCoxeterMatrix[k][j] == 2 )
							{
								for( l = 0; l < iVerticesCount; l++ )
								{
									if( iCoxeterMatrix[k][l] != 2 )
										bVerticesLinkableTemp[l] = false;
								}
								
								graphsList_euclidean->addGraph( iPathTemp, bVerticesLinkableTemp, 2, false, k, j );
							}
						}
					}
					
					// ------------------------------------------
					// ici, on a un B3, que l'on va tenter de prolonger en F4 ou un B4 que l'on va tenter de prolonger en un TF4 
					if( ( i == 1 || i == 2 ) &&  iCoxeterMatrix[ iPath[i] ][j] == 4 )
						B3ToF4_B4ToTF4( bVerticesLinkable_0_nMin1, iPathTemp, j );
					
				} // if( ( ( iCoxeterMatrix[ iPath[i] ][j] == 4 ) || ( iCoxeterMatrix[ iPath[i] ][j] == 5 && ( i == 1 || i == 2 ) ) || ( iCoxeterMatrix[ iPath[i] ][j] == 6 && i == 1 ) ) && bVerticesLinkable_0_nMin1[j] )
			}
		}
		
	}
}

void CoxIter::AnToEn_AnToTEn( const vector< short unsigned int >& iPathTemp, const vector< bool >& bVerticesLinkable )
{
	unsigned int iPathLength( iPathTemp.size() );
	
	/*
	 * 	2 pour le sommets 3 (i.e. cas sphérique: E6, E7, E8 ou car euclidien \tilde E8)
	 * 	3 pour le sommet 4 (i.e. cas euclidien: \tilde E7)
	 */
	
	bool bSpherical( iPathLength <= 7 ? true : false );
	unsigned int iStart( bSpherical || iPathLength == 8 ? 2 : 3 );
	
	// E6, E7, E8, \tilde E8
	AnToEn_AnToTEn( iPathTemp, bVerticesLinkable, bSpherical, iStart );
	
	if( iPathLength == 7 ) // \tile E7
		AnToEn_AnToTEn( iPathTemp, bVerticesLinkable, false, 3 );
}

void CoxIter::AnToEn_AnToTEn( const vector< short unsigned int >& iPathTemp, const vector< bool >& bVerticesLinkable, const bool& bSpherical, const short unsigned int& iStart )
{
	unsigned int iPathLength( iPathTemp.size() ), j, k, l;
	vector<bool> bVerticesLinkableTemp, bVerticesLinkableTempTemp;
	
	/*
	 * 	Ici, on a donc un An (n=5, 6 ou 7) avec 1 -- 2 -- 3 -- 4 -- 5 ...
	 * 	On va cherche si iStart a un voisin admissible
	 */
	for( unsigned int i(0); i < iVerticesCount; i++ )
	{
		// si le sommet est pas utilisbale (s'il l'est c'est qu'il n'est pas voisin de la base) ET si y'a un lien
		if( false == bVerticesLinkable[i] && iCoxeterMatrix[i][ iPathTemp[iStart] ] == 3 )
		{
			// on va chercher si c'est uniquement à cause d'un des sommets différents de iStart que le sommet n'est pas admissible
			for( j = 0; j < iPathLength; j++ )
			{
				if( j != iStart && iCoxeterMatrix[ iPathTemp[j] ][ i ] != 2 )
					break;
			}
			
			if( j == iPathLength ) // admissible
			{	
				bVerticesLinkableTemp = bVerticesLinkable;
				for( k = 0; k < iVerticesCount; k++ )
				{
					if( iCoxeterMatrix[i][k] != 2 )
						bVerticesLinkableTemp[k] = false;
				}
				
				if( bSpherical )
				{
					graphsList_spherical->addGraph( iPathTemp, bVerticesLinkableTemp, 4, true, i ); // En
					
					// on a un E6 qu'on va tenter de prolonger en un TE6
					if( iPathLength == 5 )
					{
						for( j = 0; j < iVerticesCount; j++ )
						{
							if( iCoxeterMatrix[i][j] == 3 && bVerticesLinkable[j] )
							{
								bVerticesLinkableTempTemp = bVerticesLinkableTemp;
								for( l = 0; l < iVerticesCount; l++ )
								{
									if( iCoxeterMatrix[j][l] != 2 )
										bVerticesLinkableTempTemp[l] = false;
								}
								
								graphsList_euclidean->addGraph( iPathTemp, bVerticesLinkableTempTemp, 4, false, i, j ); // En	
							}
						}
					}
				}
				else
					graphsList_euclidean->addGraph( iPathTemp, bVerticesLinkableTemp, 4, false, i ); // TEn
			}
		}
	}
}

/*
 * 	B3ToF4_B4ToTF4
 * 		On tente de prolonger un B3 en un F4 ou un B4 en \tilde F4 (euclidien)             
 * 										4
 * 			Le B3 donné est: iPathTemp[0] --------- iPathTemp[1] ------------ iVEnd
 * 
 * 			Le B4 donné est: iPathTemp[0] --------- iPathTemp[1] ------------ iPathTemp[2] ------------ iVEnd
 * 													4
 * 
 * 		Paramètres:
 * 			bVerticesBeginLinkable: ce qui est linkable à cause du premier (deux premiers) sommet(s)
 * 			iPathTemp: deux(trois) premiers sommets
 * 			iVEnd: 3ème(4ème) sommet
 * */
void CoxIter::B3ToF4_B4ToTF4( const vector<bool> &bVerticesBeginLinkable, vector<short unsigned int> iPathTemp, const short unsigned int &iVEnd )
{
	bool bSpherical( iPathTemp.size() == 2 ? true : false ); // true si sphérique (on cherche F4), false si euclidien (on cherche TF4)
	unsigned int i, j, iV2( iPathTemp[1] );
	vector<bool> bVerticesLinkable;

	iPathTemp.push_back( iVEnd );
	
	// on va parcourir les sommets et regarder les voisins de poids 3 de iVEnd
	for( i = 0; i < iVerticesCount; i++ )
	{
		if( iCoxeterMatrix[iVEnd][i] == 3 && bVerticesBeginLinkable[i] && iCoxeterMatrix[i][iV2] == 2 )
		{
			bVerticesLinkable = bVerticesBeginLinkable;
				
			if( !bSpherical && iCoxeterMatrix[iPathTemp[2]][i] != 2 ) // If i is connected to iPathTemp[2], this won't work
				continue;
			
			for( j = 0; j < iVerticesCount; j++ )
			{
				if( iCoxeterMatrix[j][iV2] != 2 || iCoxeterMatrix[j][iVEnd] != 2 || iCoxeterMatrix[j][i] != 2 )
					bVerticesLinkable[j] = false;
			}
			
			if( bSpherical )
				graphsList_spherical->addGraph( iPathTemp, bVerticesLinkable, 5, true, i ); // Fn
			else
				graphsList_euclidean->addGraph( iPathTemp, bVerticesLinkable, 5, false, i ); // TFn
		}
	}
}

void CoxIter::printPath()
{
	if( iPath.size() == 1 )
		return ;
	
	unsigned int iMax( iPath.size() );
	for( unsigned int i(0); i < iMax; i++ )
		cout << iPath[i] << " ; ";
	cout << endl;
}

/*! 	\fn vector2str
* 	\brief Vector --> string
* 	\param iPolynomial( const vector< mpz_class >& iPolynomial ) Polynomial
* 
* 	\remark: When mpz_class will have a to_string( mpz ) we'll move that to tools/polynomials.h TODO
*/
string vector2str( const vector< mpz_class >& iPolynomial )
{
	bool bFirst( true );
	unsigned int iSize( iPolynomial.size() );
	string strRes;
	mpz_class mpzTemp;
	
	for( unsigned int i( 0 ); i < iSize; i++ )
	{
		if( iPolynomial[i] != 0 )
		{
			if( bFirst )
			{
				strRes += iPolynomial[i].get_str() + ( i ? " * x" + string( i > 1 ? "^" + to_string( i ) : "" ) : "" );
				bFirst = false;
			}
			else
			{
				if( ( iPolynomial[i] != 1 && iPolynomial[i] != -1 ) || !i )
				{
					mpzTemp = abs( iPolynomial[i] );
					strRes += ( iPolynomial[i] > 0 ? " + " : " - " ) + mpzTemp.get_str() + ( i ? " * x" + string( i > 1 ? "^" + to_string( i ) : "" ) : "" );
				}
				else
					strRes += ( iPolynomial[i] > 0 ? " + " : " - " ) + string( i > 1 ? "x^" + to_string( i ) : "x" );
			}
		}
	}
	
	return strRes;
}

const vector< vector< GraphsProductSet > >* CoxIter::get_ptr_graphsProducts() const
{
	return &graphsProducts;
}

string CoxIter::get_strGrowthSeries_raw()
{
	string strGrowth;
	
	if( strOuputMathematicalFormat == "generic" )
	{
		strGrowth = "g(x) = (" + growthSeries_raw + ")^-1;";
	}
	else if( strOuputMathematicalFormat == "mathematica" )
	{
		strGrowth = "Symb[s_, x_] := Product[Sum[x^i, {i, 0, s[[i]] - 1}], {i, 1, Length[s]}];\n";
		strGrowth += "g[x_] := (" + growthSeries_raw + ")^-1;";
	}
	else if( strOuputMathematicalFormat == "pari" )
	{
		strGrowth = "Symb = (S, y) -> prod(i=1, length(S), sum(i=0,S[i]-1,y^i));\n";
		strGrowth += "g(x) = (" + growthSeries_raw + ")^-1;";
	}
	
	return strGrowth;
}

string CoxIter::get_strGrowthSeries()
{
	string strGrowth;
	
	if( strOuputMathematicalFormat == "generic" )
	{
		strGrowth = "f(x) = ";
		if( growthSeries_iCyclotomicNumerator.size() )
		{
			strGrowth += "C(";
			unsigned int iMax( growthSeries_iCyclotomicNumerator.size() );
			for( unsigned int i(0); i < iMax; i++ )
				strGrowth += ( i ? "," : "" ) + to_string( growthSeries_iCyclotomicNumerator[i] );
			strGrowth += ")";
		}
		
		strGrowth += "/(" + vector2str( growthSeries_iPolynomialDenominator ) + ")";
	}
	else if( strOuputMathematicalFormat == "gap" )
	{
		unsigned int iCycloSize( growthSeries_iCyclotomicNumerator.size() );
		unsigned int iDenominatorSize( growthSeries_iPolynomialDenominator.size() );
		
		strGrowth += "f := Product( [";
		for( unsigned int i(0); i < iCycloSize; i++ )
			strGrowth += ( i ? "," : "" ) + to_string( growthSeries_iCyclotomicNumerator[i] );
		strGrowth += "], i -> CyclotomicPolynomial(Rationals,i))/ValuePol( [";
		for( unsigned int i(0); i < iDenominatorSize; i++ )
			cout << ( i ? "," : "" ) << growthSeries_iPolynomialDenominator[i];
		cout << "], X(Rationals));";
	}
	else if( strOuputMathematicalFormat == "mathematica" )
	{
		unsigned int iCycloSize( growthSeries_iCyclotomicNumerator.size() );
		
		strGrowth = "Cyclo[s_, x_] := Product[Cyclotomic[s[[i]], x], {i, 1, Length[s]}];";
		
		strGrowth += "f[x_] := Cyclo[{";
		for( unsigned int i(0); i < iCycloSize; i++ )
			strGrowth += ( i ? "," : "" ) + to_string( growthSeries_iCyclotomicNumerator[i] );
		strGrowth += "},x]/(" + vector2str( growthSeries_iPolynomialDenominator ) + ");";
	}
	else if( strOuputMathematicalFormat == "pari" )
	{
		unsigned int iCycloSize( growthSeries_iCyclotomicNumerator.size() );
		strGrowth = "Cyclo = (S, y) -> prod(i=1, length(S), polcyclo(S[i],y));\n";
		strGrowth += "f(x) = Cyclo( [";
		for( unsigned int i(0); i < iCycloSize; i++ )
			strGrowth += ( i ? "," : "" ) + to_string( growthSeries_iCyclotomicNumerator[i] );
		strGrowth += "],x)/(" + vector2str( growthSeries_iPolynomialDenominator ) + ");";
	}
	
	return strGrowth;
}

void CoxIter::printGrowthSeries()
{
	if( !bGrowthSeriesComputed )
		growthSeries();
	
	if( strOuputMathematicalFormat == "generic" )
	{
		cout << "f(x) = ";
		if( growthSeries_iCyclotomicNumerator.size() )
		{
			cout << "C(";
			unsigned int iMax( growthSeries_iCyclotomicNumerator.size() );
			for( unsigned int i(0); i < iMax; i++ )
				cout << ( i ? "," : "" ) << growthSeries_iCyclotomicNumerator[i];
			cout << ")";
		}
		
		cout << "/(";
		Polynomials::polynomialDisplay( growthSeries_iPolynomialDenominator );
		cout << ")";
		
		if( bDebug )
			cout << "\ng(x) = (" << growthSeries_raw << ")^-1;";
	}
	else if( strOuputMathematicalFormat == "gap" )
	{
		unsigned int iCycloSize( growthSeries_iCyclotomicNumerator.size() );
		unsigned int iDenominatorSize( growthSeries_iPolynomialDenominator.size() );
		
		cout << "f := Product( [";
		for( unsigned int i(0); i < iCycloSize; i++ )
			cout << ( i ? "," : "" ) << growthSeries_iCyclotomicNumerator[i];
		cout << "], i -> CyclotomicPolynomial(Rationals,i))/ValuePol( [";
		for( unsigned int i(0); i < iDenominatorSize; i++ )
			cout << ( i ? "," : "" ) << growthSeries_iPolynomialDenominator[i];
		cout << "], X(Rationals));";
		
		
		if( bDebug )
			cout << "\ng(x) = (" << growthSeries_raw << ")^-1;";
	}
	else if( strOuputMathematicalFormat == "mathematica" )
	{
		unsigned int iCycloSize( growthSeries_iCyclotomicNumerator.size() );
		
		cout << "Cyclo[s_, x_] := Product[Cyclotomic[s[[i]], x], {i, 1, Length[s]}];" << endl;
		if( bDebug )
			cout << "Symb[s_, x_] := Product[Sum[x^i, {i, 0, s[[i]] - 1}], {i, 1, Length[s]}];" << endl;

		cout << "f[x_] := Cyclo[{";
		for( unsigned int i(0); i < iCycloSize; i++ )
			cout << ( i ? "," : "" ) << growthSeries_iCyclotomicNumerator[i];
		cout << "},x]";
		
		cout << "/(";
		Polynomials::polynomialDisplay( growthSeries_iPolynomialDenominator );
		cout << ");";
		
		if( bDebug )
			cout << "\ng[x_] := (" << growthSeries_raw << ")^-1;";
	}
	else if( strOuputMathematicalFormat == "pari" )
	{
		unsigned int iCycloSize( growthSeries_iCyclotomicNumerator.size() );
		
		cout << "Cyclo = (S, y) -> prod(i=1, length(S), polcyclo(S[i],y));" << endl;
		if( bDebug )
			cout << "Symb = (S, y) -> prod(i=1, length(S), sum(i=0,S[i]-1,y^i));" << endl;
		
		
		cout << "f(x) = Cyclo( [";
		for( unsigned int i(0); i < iCycloSize; i++ )
			cout << ( i ? "," : "" ) << growthSeries_iCyclotomicNumerator[i];
		cout << "],x)/(";
		Polynomials::polynomialDisplay( growthSeries_iPolynomialDenominator );
		cout << ");";
		
		if( bDebug )
			cout << "\ng(x) = (" << growthSeries_raw << ")^-1;";
	}
}

ostream& operator<<( ostream &o, const CoxIter &g )
{
	o << "Graphes sphériques: " << endl;
	o << *g.graphsList_spherical;
	
	o << "Graphes euclidien: " << endl;
	o << *g.graphsList_euclidean;
	
	return o;
}

int CoxIter::iIsGraphCocompact()
{
	if( iIsCocompact >= 0 )
		return iIsCocompact;
	
	if( !bGraphsProductsComputed )
		computeGraphsProducts();
	
	if( !bCheckCocompactness || !graphsProducts.size() || !graphsProducts[0].size() )
	{
		iIsCocompact = -1;
		return -1;
	}
	
	if( !graphsProducts[1].size() )
	{
		iIsCocompact = 0;
		return 0;
	}
	
	if( bHasBoldLine )
	{
		iIsCocompact = 0;
		return 0;
	}
	
	// ----------------------------------------------------
	// the test
	if( bUseOpenMP && iVerticesCount >= 15 )
		iIsCocompact = b_isGraph_cocompact_finiteVolume_parallel( 1 ) ? 1 : 0;
	else
		iIsCocompact = b_isGraph_cocompact_finiteVolume_sequential( 1 ) ? 1 : 0;
	
	return iIsCocompact;
}

int CoxIter::isFiniteCovolume()
{
	if( iIsFiniteCovolume >= 0 )
		return iIsFiniteCovolume;
	
	if( !bGraphsProductsComputed )
		computeGraphsProducts();
	
	// ----------------------------------------------------
	// some stupid tests
	if( !graphsProducts.size() || !graphsProducts[0].size() )
	{
		iIsFiniteCovolume = -1;
		return -1;
	}
	
	if( !graphsProducts[0].size() )
	{
		iIsFiniteCovolume = 0;
		return 0;
	}
	
	if( !graphsProducts[2].size() && !graphsProducts[1].size() ) // No vertices
	{
		iIsFiniteCovolume = 0;
		return 0;
	}
	
	/*
	* The duplication of the data is not beautiful but it allows
	* us to simplify the code of the function b_isGraph_cocompact_finiteVolume.
	* Moreover, we want the program to be quick but memory is not really an issue.
	* 
	* [2] Contains a combined list of finite and infinite vertices
	*/
	graphsProducts[2].insert( graphsProducts[2].end(), graphsProducts[1].begin(), graphsProducts[1].end() );

	// ----------------------------------------------------
	// the test
	if( bUseOpenMP && iVerticesCount >= 15 )
		iIsFiniteCovolume = b_isGraph_cocompact_finiteVolume_parallel( 2 ) ? 1 : 0;
	else
		iIsFiniteCovolume = b_isGraph_cocompact_finiteVolume_sequential( 2 ) ? 1 : 0;
	
	return iIsFiniteCovolume;
}

bool CoxIter::b_isGraph_cocompact_finiteVolume_sequential( unsigned int iIndex )
{
	unsigned int iExtendedCount, iMax( graphsProducts[0].size() ), i;
	
	vector< Graph* > vDiffSubNotBig, vDiffBigNotSub;
	vector< Graph* >::const_iterator itGSub, itGBig;
	
	bool bExtendable;
	
	for( i = 0; i < iMax; i++ )
	{
		iExtendedCount = 0;
		
		for( vector< GraphsProductSet >::const_iterator gpBig( graphsProducts[iIndex].begin() ); gpBig != graphsProducts[iIndex].end(); ++gpBig )
		{
			vDiffSubNotBig.clear();
			vDiffBigNotSub.clear();
			
			set_difference( graphsProducts[0][i].graphs.begin(), graphsProducts[0][i].graphs.end(), gpBig->graphs.begin(), gpBig->graphs.end(), back_inserter( vDiffSubNotBig ) );
			set_difference( gpBig->graphs.begin(), gpBig->graphs.end(), graphsProducts[0][i].graphs.begin(), graphsProducts[0][i].graphs.end(), back_inserter( vDiffBigNotSub ) );
			
			bExtendable = true;
			for( itGSub = vDiffSubNotBig.begin(); itGSub != vDiffSubNotBig.end(); ++itGSub )
			{
				for( itGBig = vDiffBigNotSub.begin(); itGBig != vDiffBigNotSub.end(); ++itGBig )
				{
					if( (*itGSub)->bIsSubgraphOf( *itGBig ) )
						break;
				}
				
				if( itGBig == vDiffBigNotSub.end() )
				{
					bExtendable = false;
					break;
				}
			}
			
			if( bExtendable )
				iExtendedCount++;
		}
		
		if( iExtendedCount != 2 )
		{
			if( bDebug )
			{	
				cout << "----------------------------------------------------------" << endl;
				cout << ( iIndex == 1 ? "Compactness" : "Finite covolume" ) << " test" << endl;
				cout << "Trying to extend the product: " << endl;
				cout << graphsProducts[0][i] << endl;
				cout << "Succeeded in " << iExtendedCount << " ways instead of 2" << endl;
				
				for( vector< GraphsProductSet >::const_iterator gpBig( graphsProducts[iIndex].begin() ); gpBig != graphsProducts[iIndex].end(); ++gpBig )
				{
					if( graphsProducts[0][i].b_areVerticesSubsetOf( *gpBig ) )
						cout << "Candidate: \n" << *gpBig << endl;
				}
				cout << "----------------------------------------------------------" << endl;
			}
			
			return 0;
		}
	}
	
	return 1;
}

bool CoxIter::b_isGraph_cocompact_finiteVolume_parallel( unsigned int iIndex )
{
	unsigned int iExtendedCount, iMax( graphsProducts[0].size() ), i;
	
	vector< Graph* > vDiffSubNotBig, vDiffBigNotSub;
	vector< Graph* >::const_iterator itGSub, itGBig;
	
	bool bExtendable, bExit(false);
	
	#pragma omp parallel if( bUseOpenMP && iVerticesCount >= 15 )
	{
		#pragma omp single nowait
		{
			for( i = 0; i < iMax && !bExit; i++ )
			{
				#pragma omp task private(itGBig, itGSub, bExtendable, vDiffSubNotBig, vDiffBigNotSub, iExtendedCount) shared(iMax, iIndex, bExit) firstprivate(i)
				{
					iExtendedCount = 0;
				
					for( vector< GraphsProductSet >::const_iterator gpBig( graphsProducts[iIndex].begin() ); gpBig != graphsProducts[iIndex].end(); ++gpBig )
					{
						vDiffSubNotBig.clear();
						vDiffBigNotSub.clear();
						set_difference( graphsProducts[0][i].graphs.begin(), graphsProducts[0][i].graphs.end(), gpBig->graphs.begin(), gpBig->graphs.end(), back_inserter( vDiffSubNotBig ) );
						set_difference( gpBig->graphs.begin(), gpBig->graphs.end(), graphsProducts[0][i].graphs.begin(), graphsProducts[0][i].graphs.end(), back_inserter( vDiffBigNotSub ) );

						bExtendable = true;
						for( itGSub = vDiffSubNotBig.begin(); itGSub != vDiffSubNotBig.end(); ++itGSub )
						{
							for( itGBig = vDiffBigNotSub.begin(); itGBig != vDiffBigNotSub.end(); ++itGBig )
							{
								if( (*itGSub)->bIsSubgraphOf( *itGBig ) )
									break;
							}
							
							if( itGBig == vDiffBigNotSub.end() )
							{
								bExtendable = false;
								break;
							}
						}
						
						if( bExtendable )
							iExtendedCount++;
					}
					
					if( iExtendedCount != 2 )
					{
						if( bDebug )
						{	
							#pragma omp critical
							{
								cout << "----------------------------------------------------------" << endl;
								cout << ( iIndex == 1 ? "Compactness" : "Finite covolume" ) << " test" << endl;
								cout << "Trying to extend the product: " << endl;
								cout << graphsProducts[0][i] << endl;
								cout << "Succeeded in " << iExtendedCount << " ways instead of 2" << endl;
								
								for( vector< GraphsProductSet >::const_iterator gpBig( graphsProducts[iIndex].begin() ); gpBig != graphsProducts[iIndex].end(); ++gpBig )
								{
									if( graphsProducts[0][i].b_areVerticesSubsetOf( *gpBig ) )
										cout << "Candidate: \n" << *gpBig << endl;
								}
								cout << "----------------------------------------------------------" << endl;
							}
						}
						
						#pragma omp atomic write
							bExit = true;
					}
				}
				
				if( bExit )
					break;
			}
		}
	}
	
	return ( bExit ? 0 : 1 );
}

void CoxIter::computeGraphsProducts()
{
	if( bGraphsProductsComputed )
		return;
	
	if( !bGraphExplored )
		exploreGraph();
	
	if( bDebug )
	{
		cout << "Connected spherical graphs" << endl;
		cout << *this->graphsList_spherical << endl;
	}
	
	graphsProducts = vector< vector< GraphsProductSet > >( 3 );
	vector< bool > bGPVerticesNonLinkable( vector< bool >( iVerticesCount, false ) );
	GraphsProduct gp; ///< Current graphs product
	
	// --------------------------------------------------------------
	// produits de graphes sphériques
	GraphsListIterator grIt_spherical( this->graphsList_spherical );
	
	#pragma omp parallel if( bUseOpenMP && iVerticesCount >= 15 )
	{
		#pragma omp single nowait
		while( grIt_spherical.ptr )
		{
			#pragma omp task firstprivate(grIt_spherical, bGPVerticesNonLinkable, gp)
			{
				computeGraphsProducts( grIt_spherical, &graphsProductsCount_spherical, true, gp, bGPVerticesNonLinkable );
			}
			
			++grIt_spherical;
		}
	}
	
	if( bWriteProgress )
		cout << "\tFinished" << endl;
	
	// --------------------------------------------------------------
	// produits de graphes euclidiens
	if( bDebug )
	{
		cout << "Connected euclidean graphs" << endl;
		cout << *this->graphsList_euclidean;
	}
	
	GraphsListIterator grIt_euclidean( this->graphsList_euclidean );
	#pragma omp parallel if( bUseOpenMP && iVerticesCount >= 15 )
	{
		#pragma omp single nowait
		while( grIt_euclidean.ptr )
		{
			#pragma omp task firstprivate(grIt_euclidean, bGPVerticesNonLinkable, gp)
			{
				computeGraphsProducts( grIt_euclidean, &graphsProductsCount_euclidean, false, gp, bGPVerticesNonLinkable );
			}
			
			++grIt_euclidean;
		}
	}
	
	if( bDebug )
	{
		cout << "\nProduct of euclidean graphs" << endl;
		printEuclideanGraphsProducts( &graphsProductsCount_euclidean );
	}
	
	bGraphsProductsComputed = true;
	
	// ---------------------------------------------------------
	// We guess the dimension
	if( !iDimension )
	{
		iDimension = max( iEuclideanMaxRankFound + 1, iSphericalMaxRankFound );
		bDimension_guessed = true;
		
		if( iEuclideanMaxRankFound == iSphericalMaxRankFound )
		{
			graphsProducts[0] = graphsProducts[1];
			graphsProducts[1].clear();
		}
		else if( iEuclideanMaxRankFound > iSphericalMaxRankFound )
		{
			graphsProducts[0].clear();
			graphsProducts[1].clear();
		}
		else if( iSphericalMaxRankFound > iEuclideanMaxRankFound + 1 )
		{
			graphsProducts[2].clear();
		}
	}
}

void CoxIter::computeGraphsProducts( GraphsListIterator grIt, vector< map<vector< vector< short unsigned int > >, unsigned int> > *graphsProductsCount, const bool& bSpherical, GraphsProduct& gp, vector< bool >& bGPVerticesNonLinkable )
{
	vector< short unsigned int >::iterator iIt;
	vector< short unsigned int > iVerticesFlagged;
	unsigned int iGraphRank(0);
	static unsigned int iMaxRank( iDimension ? iDimension : iVerticesCount );
	
	vector< vector< short unsigned int > > vFootPrintTest;
	
	while( grIt.ptr && ( gp.iRank + iGraphRank <= iMaxRank ) )
	{
		// ---------------------------------------------------
		// est ce que le graphe est admissible?
		for( iIt = grIt.ptr->iVertices.begin(); iIt != grIt.ptr->iVertices.end(); ++iIt )
		{
			if( bGPVerticesNonLinkable[ *iIt ] ) // si pas linkable
				break;
		}
		
		// si le graphe est admissible
		if( iIt == grIt.ptr->iVertices.end() )
		{
			// le graphe est ajouté au produit
			gp.graphs.push_back( grIt.ptr );

			// taille du graphe courant
			iGraphRank = bSpherical ? grIt.ptr->iVertices.size() : ( grIt.ptr->iVertices.size() - 1 );
			gp.iRank += iGraphRank;
			
			// Create the footprint of the product. The goal is to decide if we already have this product
			vFootPrintTest = gp.createFootPrint();
			
			#pragma omp critical
			{
				if( bCheckCocompactness || bCheckCofiniteness )
				{
					if( iDimension ) // If we know the dimension, everything is easier
					{
						if( bSpherical )
						{
							// Keeping track of spherical subgraphs
							if( ( gp.iRank == ( iDimension - 1 ) ||  gp.iRank == iDimension ) )
								graphsProducts[ gp.iRank + 1 - iDimension ].push_back( GraphsProductSet( gp ) );
						}
						
						// Euclidean subgraphs
						if( !bSpherical && gp.iRank == ( iDimension - 1 ) && bCheckCofiniteness )
							graphsProducts[2].push_back( GraphsProductSet( gp ) );
					}
					else
					{
						if( bSpherical )
						{
							if( gp.iRank == iSphericalMaxRankFound + 1 )
							{
								graphsProducts[0] = graphsProducts[1];
								graphsProducts[1].clear();
								graphsProducts[1].push_back( GraphsProductSet( gp ) );
							}
							else if( gp.iRank > iSphericalMaxRankFound + 1 )
							{
								graphsProducts[0].clear();
								graphsProducts[1].clear();
								graphsProducts[1].push_back( GraphsProductSet( gp ) );
							}
							else if( gp.iRank + 1 >= iSphericalMaxRankFound )
								graphsProducts[ gp.iRank + 1 - iSphericalMaxRankFound ].push_back( GraphsProductSet( gp ) );
						}
						else
						{	
							if( bCheckCofiniteness )
							{
								if( gp.iRank > iEuclideanMaxRankFound )
									graphsProducts[2].clear();
								
								if( gp.iRank >= iEuclideanMaxRankFound )
									graphsProducts[2].push_back( GraphsProductSet( gp ) );
							}
						}
					}
				}
				
				if( bSpherical && gp.iRank >= iSphericalMaxRankFound )
					iSphericalMaxRankFound = gp.iRank;
				
				if( !bSpherical && gp.iRank >= iEuclideanMaxRankFound )
					iEuclideanMaxRankFound = gp.iRank;
				
				if( (*graphsProductsCount)[ gp.iRank ].find( vFootPrintTest ) == (*graphsProductsCount)[ gp.iRank ].end() )
					(*graphsProductsCount)[ gp.iRank ][ vFootPrintTest ] = 1;
				else
					(*graphsProductsCount)[ gp.iRank ][ vFootPrintTest ]++;
			}
			
			// mise à jour des sommets que l'on ne peut plus prendre
			for( unsigned int i = 0; i < iVerticesCount; i++ )
			{
				if( !grIt.ptr->bVerticesLinkable[i] && !bGPVerticesNonLinkable[i] )
				{
					iVerticesFlagged.push_back( i );
					bGPVerticesNonLinkable[i] = true;		
				}
			}
			
			// récursion
			computeGraphsProducts( ++grIt, graphsProductsCount, bSpherical, gp, bGPVerticesNonLinkable );
			
			// -----------------------------------------------
			// dé-initialisations
			
			// on remet la liste à son état d'avant la récursion	
			for( iIt = iVerticesFlagged.begin(); iIt != iVerticesFlagged.end(); ++ iIt )
				bGPVerticesNonLinkable[ *iIt ] = false;
			
			gp.iRank -= iGraphRank;
			
			// le graphe est enlevé
			gp.graphs.pop_back();
			
			if( !gp.graphs.size() )
				break;
		}
		else
			++grIt;
	}
}

bool CoxIter::bCanBeFiniteCovolume()
{
	// -----------------------------------------------------------
	// Some verifications
	if( !iDimension )
		throw( string( "CoxIter::bCanBeFiniteCovolume: Dimension not specified" ) );
	
	if( !bGraphExplored )
		exploreGraph();
	
	// -----------------------------------------------------------
	// Initializations
	graphsProducts_bCanBeFiniteCovolume = vector< vector< GraphsProductSet > >(1);
	
	GraphsListIterator grIt_euclidean( this->graphsList_euclidean );
	vector< bool > bGPVerticesNonLinkable( vector< bool >( iVerticesCount, false ) );
	GraphsProduct gp; ///< Current graphs product
	
	// -----------------------------------------------------------
	// We find the products of euclidean graphs
	#pragma omp parallel if( bUseOpenMP && iVerticesCount >= 15 )
	{
		#pragma omp single nowait
		while( grIt_euclidean.ptr )
		{
			#pragma omp task firstprivate(grIt_euclidean, bGPVerticesNonLinkable, gp)
			{
				bCanBeFiniteCovolume_computeGraphsProducts( grIt_euclidean, gp, bGPVerticesNonLinkable );
			}
			
			++grIt_euclidean;
		}
	}
	
	// -----------------------------------------------------------
	// Check the condition: every connected affine graph of rank at least 2 is a subgraph of an affine graph of rank n-1
	grIt_euclidean = GraphsListIterator( GraphsListIterator( this->graphsList_euclidean, 3 ) ); // TODO: partir à 2?
	bool bCanBeExtended;
	
	while( grIt_euclidean.ptr )
	{
		bCanBeExtended = false;
		for( auto graphProd : graphsProducts_bCanBeFiniteCovolume[0] ) // TODO OPTIMIZATION paralleliser?
		{
			for( auto gr : graphProd.graphs )
			{
				// For two affine graphs G1 and G2, G1 is a subgraph of G2 iff G1=G2
				if( *grIt_euclidean.ptr == *gr )
				{
					bCanBeExtended = true;
					break;
				}
			}
			
			if( bCanBeExtended )
				break;
		}
		
		if( !bCanBeExtended )
		{
			if( bDebug )
			{
				cout << "Can be of finite covolume: no" << endl;
				cout << "\tCannot extend the affine graph: " << endl;
				cout << "\t" << *grIt_euclidean.ptr << "\n" << endl;
			}
			
			return false;
		}
		
		++grIt_euclidean;
	}
	
	return true;
}

void CoxIter::bCanBeFiniteCovolume_computeGraphsProducts(GraphsListIterator grIt, GraphsProduct& gp, vector< bool >& bGPVerticesNonLinkable)
{
	vector< short unsigned int >::iterator iIt;
	vector< short unsigned int > iVerticesFlagged;
	unsigned int iGraphRank(0);
	
	while( grIt.ptr && ( gp.iRank + iGraphRank <= iVerticesCount ) )
	{
		// ---------------------------------------------------
		// est ce que le graphe est admissible?
		for( iIt = grIt.ptr->iVertices.begin(); iIt != grIt.ptr->iVertices.end(); ++iIt )
		{
			if( bGPVerticesNonLinkable[ *iIt ] ) // si pas linkable
				break;
		}
		
		// si le graphe est admissible
		if( iIt == grIt.ptr->iVertices.end() )
		{
			// le graphe est ajouté au produit
			gp.graphs.push_back( grIt.ptr );

			// taille du graphe courant
			iGraphRank = grIt.ptr->iVertices.size() - 1;
			gp.iRank += iGraphRank;
			
			#pragma omp critical
			{
				if( gp.iRank == ( iDimension - 1 ) )
					graphsProducts_bCanBeFiniteCovolume[0].push_back( GraphsProductSet( gp ) );
			}
			
			// mise à jour des sommets que l'on ne peut plus prendre
			for( unsigned int i = 0; i < iVerticesCount; i++ )
			{
				if( !grIt.ptr->bVerticesLinkable[i] && !bGPVerticesNonLinkable[i] )
				{
					iVerticesFlagged.push_back( i );
					bGPVerticesNonLinkable[i] = true;		
				}
			}
			
			// récursion
			bCanBeFiniteCovolume_computeGraphsProducts( ++grIt, gp, bGPVerticesNonLinkable );
			
			// -----------------------------------------------
			// dé-initialisations
			
			// on remet la liste à son état d'avant la récursion	
			for( iIt = iVerticesFlagged.begin(); iIt != iVerticesFlagged.end(); ++ iIt )
				bGPVerticesNonLinkable[ *iIt ] = false;
			
			gp.iRank -= iGraphRank;
			
			// le graphe est enlevé
			gp.graphs.pop_back();
			
			if( !gp.graphs.size() )
				break;
		}
		else
			++grIt;
	}
}

vector< vector< short unsigned int > > CoxIter::bCanBeFiniteCovolume_complete()
{
	// -----------------------------------------------------------
	// Some verifications
	if( !iDimension )
		throw( string( "CoxIter::bCanBeFiniteCovolume_complete: Dimension not specified" ) );
	
	if( !bGraphExplored )
		exploreGraph();
	
	// -----------------------------------------------------------
	// Initializations
	graphsProducts_bCanBeFiniteCovolume = vector< vector< GraphsProductSet > >( iDimension + 1 );
	
	GraphsListIterator grIt_euclidean( this->graphsList_euclidean );
	vector< bool > bGPVerticesNonLinkable( vector< bool >( iVerticesCount, false ) );
	GraphsProduct gp; ///< Current graphs product
	
	// -----------------------------------------------------------
	// We find the products of euclidean graphs
	#pragma omp parallel if( bUseOpenMP && iVerticesCount >= 15 )
	{
		#pragma omp single nowait
		while( grIt_euclidean.ptr )
		{
			
			#pragma omp task firstprivate(grIt_euclidean, bGPVerticesNonLinkable, gp)
			{
				bCanBeFiniteCovolume_complete_computeGraphsProducts( grIt_euclidean, gp, bGPVerticesNonLinkable );
			}
			
			++grIt_euclidean;
		}
	}
	
	vector< vector< short unsigned int > > iGraphsNotExtendable( 0 );
	
	// -----------------------------------------------------------
	// Check the condition: every connected affine graph of rank at least 2 is a subgraph of an affine graph of rank n-1
	bool bCanBeExtended, bSubproduct, bFound;
	
	for( unsigned int i(2); i < iDimension - 1; i++ )
	{
		for( auto gpSmall : graphsProducts_bCanBeFiniteCovolume[i] )
		{
			// If gpSmall is not a subgraph of gpBig for every gpBig (affine of rank iDimension-1), then the graph is not of finite covolume
			bCanBeExtended = false;
			
			for( auto gpBig : graphsProducts_bCanBeFiniteCovolume[iDimension - 1] )
			{
				// First test
				if( !gpSmall.b_areVerticesSubsetOf( gpBig ) )
					continue;
				
				bSubproduct = true;
				for( auto gSmall : gpSmall.graphs )
				{
					bFound = false;
					for( auto gBig : gpBig.graphs )
					{
						if( *gSmall == *gBig )
						{
							bFound = true;
							break;
						}
					}
					
					if( !bFound )
					{
						bSubproduct = false;
						break;
					}
				}
				
				if( bSubproduct )
				{
					bCanBeExtended = true;
					break;
				}
			}
			
			if( !bCanBeExtended )
			{
				iGraphsNotExtendable.push_back( gpSmall.get_iVertices() );
			}
		}
	}
	
	return iGraphsNotExtendable;
}

void CoxIter::bCanBeFiniteCovolume_complete_computeGraphsProducts(GraphsListIterator grIt, GraphsProduct& gp, vector< bool >& bGPVerticesNonLinkable)
{
	vector< short unsigned int >::iterator iIt;
	vector< short unsigned int > iVerticesFlagged;
	unsigned int iGraphRank(0);
	
	while( grIt.ptr && ( gp.iRank + iGraphRank <= iVerticesCount ) )
	{
		// ---------------------------------------------------
		// est ce que le graphe est admissible?
		for( iIt = grIt.ptr->iVertices.begin(); iIt != grIt.ptr->iVertices.end(); ++iIt )
		{
			if( bGPVerticesNonLinkable[ *iIt ] ) // si pas linkable
				break;
		}
		
		// si le graphe est admissible
		if( iIt == grIt.ptr->iVertices.end() )
		{
			// le graphe est ajouté au produit
			gp.graphs.push_back( grIt.ptr );

			// taille du graphe courant
			iGraphRank = grIt.ptr->iVertices.size() - 1;
			gp.iRank += iGraphRank;
			
			#pragma omp critical
			{
				if( 2 <= gp.iRank && gp.iRank <= ( iDimension - 1 ) )
					graphsProducts_bCanBeFiniteCovolume[gp.iRank].push_back( GraphsProductSet( gp ) );
			}
			
			// mise à jour des sommets que l'on ne peut plus prendre
			for( unsigned int i = 0; i < iVerticesCount; i++ )
			{
				if( !grIt.ptr->bVerticesLinkable[i] && !bGPVerticesNonLinkable[i] )
				{
					iVerticesFlagged.push_back( i );
					bGPVerticesNonLinkable[i] = true;		
				}
			}
			
			// récursion
			bCanBeFiniteCovolume_complete_computeGraphsProducts( ++grIt, gp, bGPVerticesNonLinkable );
			
			// -----------------------------------------------
			// dé-initialisations
			
			// on remet la liste à son état d'avant la récursion	
			for( iIt = iVerticesFlagged.begin(); iIt != iVerticesFlagged.end(); ++ iIt )
				bGPVerticesNonLinkable[ *iIt ] = false;
			
			gp.iRank -= iGraphRank;
			
			// le graphe est enlevé
			gp.graphs.pop_back();
			
			if( !gp.graphs.size() )
				break;
		}
		else
			++grIt;
	}
}

mpz_class CoxIter::i_orderFiniteSubgraph(const unsigned int& iType, const unsigned int& iDataSupp )
{
	if( iType == 0 ) // A_n
		return iFactorials[ iDataSupp + 1 ];
	else if( iType == 1 ) // Bn
		return ( iFactorials[ iDataSupp ] * iPowersOf2[ iDataSupp ] );
	else if( iType == 3 ) // Dn
		return ( iFactorials[ iDataSupp ] * iPowersOf2[ iDataSupp - 1 ] );
	else  if( iType == 4 )
	{
		if( iDataSupp == 6 )
			return 51840;
		else if( iDataSupp == 7 )
			return 2903040;
		else if( iDataSupp == 8 )
			return 696729600;
	}
	else if( iType == 5 ) // F4
		return 1152;
	else if( iType == 6 ) // G_2^n
		return ( 2 * iDataSupp );
	else if( iType == 7 )
	{
		if( iDataSupp == 3 )
			return 120;
		else if( iDataSupp )
			return 14400;
	}
	else
		throw( 0 );
	
	return 0;
}

void CoxIter::growthSeries_mergeTerms( vector< mpz_class >& iPolynomial, vector< unsigned int >& iSymbol, vector< mpz_class > iTemp_polynomial, const vector< unsigned int >& iTemp_symbol, mpz_class biTemp )
{
	unsigned int iSymbol_max( iSymbol.size() ? iSymbol.size() - 1 : 0 );
	unsigned int iTemp_symbol_max( iTemp_symbol.size() - 1 );
	
	vector< unsigned int > iTemp_symbolDenominatorTemp( iTemp_symbol );
	
	if( iTemp_symbol_max < iSymbol_max )
		iTemp_symbolDenominatorTemp.insert( iTemp_symbolDenominatorTemp.end(), iSymbol_max - iTemp_symbol_max, 0 );
	
	// First step to compute the lcm of the two symbols
	for( unsigned int i(1); i <= iTemp_symbol_max; i++ )
	{
		for( unsigned int j( i > iSymbol_max ? 1 : iSymbol[i] + 1 ); j <= iTemp_symbol[i]; j++ )
			Polynomials::polynomialDotSymbol( iPolynomial, i );
	}
	
	// Second step to compute the lcm of the two symbols
	for( unsigned int i(1); i <= iSymbol_max; i++ )
	{
		iTemp_symbolDenominatorTemp[i] = i > iTemp_symbol_max ? iSymbol[i] : max( iTemp_symbol[i], iSymbol[i] );
		
		for( unsigned int j( i > iTemp_symbol_max ? 1 : iTemp_symbol[i] + 1 ); j <= iSymbol[i]; j++ )
			Polynomials::polynomialDotSymbol( iTemp_polynomial, i );
	}
	
	// we eventually add some zeroes
	if( iPolynomial.size() < iTemp_polynomial.size() )
		iPolynomial.insert( iPolynomial.end(), iTemp_polynomial.size() - iPolynomial.size(), 0 );

	unsigned int iTempPolynomialDegree( iTemp_polynomial.size() - 1 );
	
	// Addition of the two numerators
	for( unsigned int i( 0 ); i <= iTempPolynomialDegree; i++ )
		iPolynomial[i] += iTemp_polynomial[i] * biTemp;
	
	// ----------------------------------------------------
	// Final stuff
	iSymbol = iTemp_symbolDenominatorTemp;

	// We remove final 0
	while( iPolynomial[ iTempPolynomialDegree ] == 0 )
		iTempPolynomialDegree--;
	iPolynomial.erase( iPolynomial.begin() + iTempPolynomialDegree + 1, iPolynomial.end() );
	
}

void CoxIter::growthSeries()
{
	if( !bUseOpenMP || iVerticesCount < 10 )
		growthSeries_sequential();
	else
		growthSeries_parallel();
	
	if( bDebug )
		growthSeries_details();
}

void CoxIter::growthSeries_details()
{
	if( !bGraphExplored )
		exploreGraph();
	
	if( !bGraphsProductsComputed )
		computeGraphsProducts();
	
	unsigned int iSizeMax( graphsProductsCount_spherical.size() );
	
	unsigned int iExponent; // Temporary exponent
	
	growthSeries_raw = "1";
	
	string strGrowth_temp, strSymbol;
	
	for( unsigned int iSize( 1 ); iSize < iSizeMax; iSize++ ) // For each size
	{
		strGrowth_temp = "";
		for( auto iProduct : graphsProductsCount_spherical[iSize] ) // For each product of that size
		{
			// ----------------------------------------------------
			// Preliminary stuff
			
			// We compute the symbol and the exponent of this product
			growthSeries_symbolExponentFromProduct( iProduct.first, strSymbol, iExponent );
			
			mpz_class biTemp( (int) iProduct.second * ( (iSize % 2) ? -1 : 1 ) );
			
			strGrowth_temp += ( strGrowth_temp == "" ? "" : " + " ) + ( iProduct.second == 1 ? "" : to_string( iProduct.second ) + " * " ) + "x^" + to_string( iExponent ) + "/";
			
			if( strOuputMathematicalFormat == "mathematica" )
				strGrowth_temp += "Symb[{" + strSymbol + "},x]";
			else if( strOuputMathematicalFormat == "pari" )
				strGrowth_temp += "Symb([" + strSymbol + "],x)";
			else
				strGrowth_temp += "[" + strSymbol + "]";
		}
		
		if( strGrowth_temp != "" )
			growthSeries_raw += ( (iSize % 2) ? " - (" : " + (" ) + strGrowth_temp + ")";
	}
}

void CoxIter::growthSeries_sequential()
{
	if( !bGraphExplored )
		exploreGraph();
	
	if( !bGraphsProductsComputed )
		computeGraphsProducts();
	
	unsigned int iSizeMax( graphsProductsCount_spherical.size() );
	
	vector< unsigned int > iSymbol; // Temporary symbol
	unsigned int iSymbolMax; // Size of the temporary symbol
	unsigned int iExponent; // Temporary exponent
	
	vector< unsigned int > growthSeries_iSymbolNumerator;
	growthSeries_iPolynomialDenominator = vector< mpz_class >( { 1 } );
	growthSeries_iCyclotomicNumerator.clear();
	growthSeries_bFractionReduced = true;
	
	vector< unsigned int > iSymbolDenominatorTemp;
	unsigned int iSymbolDenominatorMax( 0 );
	
	for( unsigned int iSize( 1 ); iSize < iSizeMax; iSize++ ) // For each size
	{
		for( auto iProduct : graphsProductsCount_spherical[iSize] ) // For each product of that size
		{
			// ----------------------------------------------------
			// Preliminary stuff
			
			// We compute the symbol and the exponent of this product
			growthSeries_symbolExponentFromProduct( iProduct.first, iSymbol, iExponent );
			
			iSymbolMax = iSymbol.size() - 1;
			iSymbolDenominatorTemp = iSymbol;
			
			mpz_class biTemp( (int) iProduct.second * ( (iSize % 2) ? -1 : 1 ) );
			
			vector< mpz_class > iTempPolynomial( vector< mpz_class >( iExponent, 0 ) );
			iTempPolynomial.push_back( 1 ); // x^iExponent
			
			for( unsigned int i( iSymbolMax + 1 ); i <= iSymbolDenominatorMax; i++ ) // we add some zeroes
				iSymbolDenominatorTemp.push_back( 0 );
			
			// ----------------------------------------------------
			// Update
			/*
			 * We have here the current rational function: growthSeries_iPolynomialDenominator / iSymbolDenominator
			 * We want to add the rational function: (-1)^iSize * iProduct.second * x^iExponent / iSymbol
			 */
			
			// First step to compute the lcm of the two symbols iSymbolDenominator and iSymbol
			for( unsigned int i(1); i <= iSymbolMax; i++ )
			{
				for( unsigned int j( i > iSymbolDenominatorMax ? 1 : growthSeries_iSymbolNumerator[i] + 1 ); j <= iSymbol[i]; j++ )
					Polynomials::polynomialDotSymbol( growthSeries_iPolynomialDenominator, i );
			}
			
			// Second step to compute the lcm of the two symbols iSymbolDenominator and iSymbol
			for( unsigned int i(1); i <= iSymbolDenominatorMax; i++ )
			{
				//cout << "\tiTemp_symbol (" << i << "): " << implode( ",", iSymbol ) << endl;
			
				iSymbolDenominatorTemp[i] = max( iSymbol[i], growthSeries_iSymbolNumerator[i] );
				
				for( unsigned int j( i > iSymbolMax ? 1 : iSymbol[i] + 1 ); j <= growthSeries_iSymbolNumerator[i]; j++ )
					Polynomials::polynomialDotSymbol( iTempPolynomial, i );
			}
			
			// we eventually add some zeroes
			if( growthSeries_iPolynomialDenominator.size() < iTempPolynomial.size() )
				growthSeries_iPolynomialDenominator.insert( growthSeries_iPolynomialDenominator.end(), iTempPolynomial.size() - growthSeries_iPolynomialDenominator.size(), 0 );
			
			unsigned int iTempPolynomialDegree( iTempPolynomial.size() - 1 );
			
			// Addition of the two numerators
			for( unsigned int i( 0 ); i <= iTempPolynomialDegree; i++ )
				growthSeries_iPolynomialDenominator[i] += iTempPolynomial[i] * biTemp;
			
			// ----------------------------------------------------
			// Final stuff
			growthSeries_iSymbolNumerator = iSymbolDenominatorTemp;
			iSymbolDenominatorMax = growthSeries_iSymbolNumerator.size() - 1;

			// We remove final 0
			while( growthSeries_iPolynomialDenominator[ iTempPolynomialDegree ] == 0 )
				iTempPolynomialDegree--;
			growthSeries_iPolynomialDenominator.erase( growthSeries_iPolynomialDenominator.begin() + iTempPolynomialDegree + 1, growthSeries_iPolynomialDenominator.end() );
		}
	}
	
	// --------------------------------------------------------------
	// Symbols --> Cyclotomic polynomials
	vector< unsigned int > iCyclotomicTemp;
	
	for( unsigned int i( iSymbolDenominatorMax ); i >= 2; i-- )
	{
		if( growthSeries_iSymbolNumerator[i] )
		{
			auto iDivisors( iListDivisors( i, true ) );
			iDivisors.push_back(i);
		
			for( unsigned int j(1); j <= growthSeries_iSymbolNumerator[i]; j++ )
				iCyclotomicTemp.insert( iCyclotomicTemp.end(), iDivisors.begin(), iDivisors.end() );
		}
	}
	
	// --------------------------------------------------------------
	// Simplifications
	unsigned int iCyclotomicTempSize( iCyclotomicTemp.size() ), iCyclotomicMax( Polynomials::iCyclotomicPolynomials.size() - 1 );
	for( unsigned int i(0); i < iCyclotomicTempSize; i++ )
	{
		if( iCyclotomicMax < iCyclotomicTemp[i] || !Polynomials::dividePolynomialByPolynomial( growthSeries_iPolynomialDenominator, Polynomials::iCyclotomicPolynomials[ iCyclotomicTemp[i] ] ) )
			growthSeries_iCyclotomicNumerator.push_back( iCyclotomicTemp[i] );
		
		if( iCyclotomicMax < iCyclotomicTemp[i] )
			growthSeries_bFractionReduced = false;
	}
	
	// --------------------------------------------------------------
	// Final stuff
	
	// We remove final 0
	while( !growthSeries_iSymbolNumerator[iSymbolDenominatorMax] )
		iSymbolDenominatorMax--;
	growthSeries_iSymbolNumerator.erase( growthSeries_iSymbolNumerator.begin() + iSymbolDenominatorMax + 1, growthSeries_iSymbolNumerator.end() );
	
	sort( growthSeries_iCyclotomicNumerator.begin(), growthSeries_iCyclotomicNumerator.end() );
	
	bGrowthSeriesComputed = true;
}


void CoxIter::growthSeries_parallel()
{
	if( !bGraphExplored )
		exploreGraph();
	
	if( !bGraphsProductsComputed )
		computeGraphsProducts();
	
	unsigned int iSizeMax( graphsProductsCount_spherical.size() );

	growthSeries_iPolynomialDenominator.clear();
	growthSeries_iCyclotomicNumerator.clear();
	growthSeries_bFractionReduced = true;
	
	// -----------------------------------------------------------------
	// Private and local variables
	vector< unsigned int > iTemp_symbolDenominatorTemp;
	vector< unsigned int > iSymbol; // Temporary symbol
	unsigned int iExponent; // Temporary exponent
	unsigned int iThreadId;
	
	// -----------------------------------------------------------------
	// Shared and local variables
	int iOMPMaxThreads( omp_get_max_threads() );
	
	vector< vector< mpz_class > > gs_iPolynomialDenominator( iOMPMaxThreads, vector< mpz_class >( {0} ) );
	vector< vector< unsigned int > > gs_iSymbolNumerator( iOMPMaxThreads, vector< unsigned int >( 0 ) );
	
	gs_iPolynomialDenominator[0][0] = 1; // Master thread, empty set --> trivial subgroup
	
	#pragma omp parallel for default(none) shared(iSizeMax, gs_iSymbolNumerator, gs_iPolynomialDenominator) private(iExponent, iSymbol, iThreadId, iTemp_symbolDenominatorTemp) schedule(static,1)
	for( unsigned int iSize = iSizeMax - 1; iSize >= 1; iSize-- ) // For each size
	{
		iThreadId = omp_get_thread_num();
		
		for( auto iProduct : graphsProductsCount_spherical[iSize] ) // For each product of that size
		{
			// ----------------------------------------------------
			// Preliminary stuff
			
			// We compute the symbol and the exponent of this product
			growthSeries_symbolExponentFromProduct( iProduct.first, iSymbol, iExponent );
			mpz_class biTemp( (int) iProduct.second * ( (iSize % 2) ? -1 : 1 ) );
			
			vector< mpz_class > iTemp_polynomial( vector< mpz_class >( iExponent, 0 ) );
			iTemp_polynomial.push_back( 1 ); // x^iExponent
			
			growthSeries_mergeTerms( gs_iPolynomialDenominator[iThreadId], gs_iSymbolNumerator[iThreadId], iTemp_polynomial, iSymbol, biTemp );
		}
	}
	
	// --------------------------------------------------------------
	// Reduction
	growthSeries_iPolynomialDenominator = gs_iPolynomialDenominator[0];
	for( int iThread(1); iThread < iOMPMaxThreads; iThread++ )
	{
		if( gs_iSymbolNumerator[iThread].size() )
			growthSeries_mergeTerms( growthSeries_iPolynomialDenominator, gs_iSymbolNumerator[0], gs_iPolynomialDenominator[iThread], gs_iSymbolNumerator[iThread] );
	}
	
	// --------------------------------------------------------------
	// Symbols --> Cyclotomic polynomials
	vector< unsigned int > iCyclotomicTemp;
	unsigned int iSymbolDenominatorMax( gs_iSymbolNumerator[0].size() - 1 );
	
	for( unsigned int i( iSymbolDenominatorMax ); i >= 2; i-- )
	{
		if( gs_iSymbolNumerator[0][i] )
		{
			auto iDivisors( iListDivisors( i, true ) );
			iDivisors.push_back(i);
		
			for( unsigned int j(1); j <= gs_iSymbolNumerator[0][i]; j++ )
				iCyclotomicTemp.insert( iCyclotomicTemp.end(), iDivisors.begin(), iDivisors.end() );
		}
	}
	
	// --------------------------------------------------------------
	// Simplifications
	unsigned int iCyclotomicTempSize( iCyclotomicTemp.size() ), iCyclotomicMax( Polynomials::iCyclotomicPolynomials.size() - 1 );
	for( unsigned int i(0); i < iCyclotomicTempSize; i++ )
	{
		if( iCyclotomicMax < iCyclotomicTemp[i] || !Polynomials::dividePolynomialByPolynomial( growthSeries_iPolynomialDenominator, Polynomials::iCyclotomicPolynomials[ iCyclotomicTemp[i] ] ) )
			growthSeries_iCyclotomicNumerator.push_back( iCyclotomicTemp[i] );
		
		if( iCyclotomicMax < iCyclotomicTemp[i] )
			growthSeries_bFractionReduced = false;
	}
	
	// --------------------------------------------------------------
	// Final stuff
	sort( growthSeries_iCyclotomicNumerator.begin(), growthSeries_iCyclotomicNumerator.end() );
	bGrowthSeriesComputed = true;
}

void CoxIter::get_iGrowthSeries( vector< unsigned int >& iCyclotomicNumerator, vector< mpz_class >& iPolynomialDenominator, bool& bReduced )
{
	if( !bGrowthSeriesComputed )
		growthSeries();

	iCyclotomicNumerator = growthSeries_iCyclotomicNumerator;
	iPolynomialDenominator = growthSeries_iPolynomialDenominator;
	bReduced = growthSeries_bFractionReduced;
}

bool CoxIter::get_bGrowthSeriesReduced()
{
	if( !bGrowthSeriesComputed )
		growthSeries();

	return growthSeries_bFractionReduced;
}

vector< mpz_class > CoxIter::get_iGrowthSeries_denominator()
{
	if( !bGrowthSeriesComputed )
		growthSeries();
	
	return growthSeries_iPolynomialDenominator;
}

void CoxIter::growthSeries_symbolExponentFromProduct(const vector< vector< short unsigned int > >& iProduct, string& strSymbol, unsigned int& iExponent) const
{
	short unsigned int j, jMax, k;
	
	vector< unsigned int > iSymbol;
	iExponent = 0;
	
	for( unsigned int i( 0 ); i < 8; i++ ) // For each type of irreducible spherical graph
	{
		jMax = iProduct[i].size();
		
		// pour chaque taille
		for( j = 0; j < jMax; j++ )
		{
			if( !iProduct[i][j] )
				continue;
			
			vector< unsigned int > iSymbolTemp;
			
			switch( i )
			{
				case 0: // An
					iExponent += iProduct[i][j] * ( j + 1 ) * ( j + 2 ) / 2;
					for( k = 2; k <= j + 2; k++ )
						iSymbolTemp.push_back( k );
					break;
				case 1: // Bn
					iExponent += iProduct[i][j] * ( j + 1 ) * ( j + 1 );
					for( k = 1; k <= j + 1; k++ )
						iSymbolTemp.push_back( 2 * k );
					break;
				case 3: // Dn
					iExponent += iProduct[i][j] * j * (j + 1);
					for( k = 1; k <= j; k++ )
						iSymbolTemp.push_back( 2 * k );
					iSymbolTemp.push_back( j + 1 );
					break;
				case 4: // En:
					if( j == 5 ) // E6
					{
						iSymbolTemp = vector< unsigned int >{ 2,5,6,8,9,12 };
						iExponent += iProduct[i][j] * 36;
					}
					else if( j == 6 )
					{
						iSymbolTemp = vector< unsigned int >{ 2,6,8,10,12,14,18 };
						iExponent += iProduct[i][j] * 63;
					}
					else if( j == 7 )
					{
						iSymbolTemp = vector< unsigned int >{ 2,8,12,14,18,20,24,30 };
						iExponent += iProduct[i][j] * 120;
					}
					break;
				case 5: // F4
					iSymbolTemp = vector< unsigned int >{ 2,6,8,12 };
					iExponent += iProduct[i][j] * 24;
					break;
				case 6: // G_2^m
					iExponent += iProduct[i][j] * (j + 1);
					iSymbolTemp.push_back( 2 );
					iSymbolTemp.push_back( j + 1 );
					break;
				case 7:
					if( j == 2 ) // H3
					{
						iSymbolTemp = vector< unsigned int >{ 2,6,10 };
						iExponent += iProduct[i][j] * 15;
					}
					else if( j == 3 )
					{
						iSymbolTemp = vector< unsigned int >{ 2,12,20,30 };
						iExponent += iProduct[i][j] * 60;
					}
					break;
			}
			
			for( k = 0; k < iProduct[i][j]; k++ )
				iSymbol.insert( iSymbol.end(), iSymbolTemp.begin(), iSymbolTemp.end() );
		}
	}
	
	sort( iSymbol.begin(), iSymbol.end() );
	
	strSymbol = implode( ",", iSymbol );
}

void CoxIter::growthSeries_symbolExponentFromProduct(const vector< vector< short unsigned int > >& iProduct, vector< unsigned int >& iSymbol, unsigned int& iExponent) const
{
	unsigned int j, jMax, k;

	iSymbol.clear();
	iExponent = 0;
	
	for( unsigned int i( 0 ); i < 8; i++ ) // For each type of irreducible spherical graph
	{
		jMax = iProduct[i].size();
		
		// pour chaque taille
		for( j = 0; j < jMax; j++ )
		{
			if( !iProduct[i][j] )
				continue;
			
			vector< unsigned int > iSymbolTemp;
			
			switch( i )
			{
				case 0: // An
					iExponent += iProduct[i][j] * ( j + 1 ) * ( j + 2 ) / 2;
					for( k = 2; k <= j + 2; k++ )
						iSymbolTemp.push_back( k );
					break;
				case 1: // Bn
					iExponent += iProduct[i][j] * ( j + 1 ) * ( j + 1 );
					for( k = 1; k <= j + 1; k++ )
						iSymbolTemp.push_back( 2 * k );
					break;
				case 3: // Dn
					iExponent += iProduct[i][j] * j * (j + 1);
					for( k = 1; k <= j; k++ )
						iSymbolTemp.push_back( 2 * k );
					iSymbolTemp.push_back( j + 1 );
					break;
				case 4: // En:
					if( j == 5 ) // E6
					{
						iSymbolTemp = vector< unsigned int >{ 2,5,6,8,9,12 };
						iExponent += iProduct[i][j] * 36;
					}
					else if( j == 6 )
					{
						iSymbolTemp = vector< unsigned int >{ 2,6,8,10,12,14,18 };
						iExponent += iProduct[i][j] * 63;
					}
					else if( j == 7 )
					{
						iSymbolTemp = vector< unsigned int >{ 2,8,12,14,18,20,24,30 };
						iExponent += iProduct[i][j] * 120;
					}
					break;
				case 5: // F4
					iSymbolTemp = vector< unsigned int >{ 2,6,8,12 };
					iExponent += iProduct[i][j] * 24;
					break;
				case 6: // G_2^m
					iExponent += iProduct[i][j] * (j + 1);
					iSymbolTemp.push_back( 2 );
					iSymbolTemp.push_back( j + 1 );
					break;
				case 7:
					if( j == 2 ) // H3
					{
						iSymbolTemp = vector< unsigned int >{ 2,6,10 };
						iExponent += iProduct[i][j] * 15;
					}
					else if( j == 3 )
					{
						iSymbolTemp = vector< unsigned int >{ 2,12,20,30 };
						iExponent += iProduct[i][j] * 60;
					}
					break;
			}
			
			for( auto symb : iSymbolTemp )
			{
				for( k = iSymbol.size(); k <= symb; k++ )
					iSymbol.push_back( 0 );
				
				iSymbol[ symb ] += iProduct[i][j];
			}
		}
	}
}

bool CoxIter::bEulerCharacteristicFVector()
{
	// variables de boucles
	size_t i, j, k, iMax;
	map<vector< vector< short unsigned int > >, unsigned int>::iterator itMap;
	
	mpz_class biTemp, biOrderTemp;
	MPZ_rational brAlternateTemp;
	
	bool bPositive(true);
	
	iFVector = vector< unsigned int >( iDimension + 1, 0 );
	int iFVectorIndex( iDimension );
	
	unsigned int iCurrentVerticesCount( 0 );
	
	iFVectorAlternateSum = 0;
	brEulerCaracteristic = 1;
	strEulerCharacteristic_computations = "1";
	
	if( bDebug )
		cout << "\nProducts of spherical graphs" << endl;
	
	iFVector[ iDimension ] = 1;
	
	// par taille de nombre de sommets
	for( vector< map<vector< vector< short unsigned int > >, unsigned int> >::iterator itMaps( graphsProductsCount_spherical.begin() ); itMaps != graphsProductsCount_spherical.end(); ++itMaps )
	{
		brAlternateTemp = 0;
		
		// on parcourt les produits pour la taille donnée
		for( itMap = itMaps->begin(); itMap != itMaps->end(); ++itMap )
		{
			biTemp = 1;

			if( bDebug )
				cout << "\t" << iCurrentVerticesCount << ": ";
			
			// pour chaque type de graphe
			for( i = 0; i < 8; i++ )
			{
				iMax = (itMap->first[i]).size();
				
				// pour chaque taille
				for( j = 0; j < iMax; j++ )
				{
					if( itMap->first[i][j] )
					{
						if( bDebug )
							cout << (char)(i + 65) << "_" << ( j + 1 ) << "^" << itMap->first[i][j] << " | ";
						
						biOrderTemp = i_orderFiniteSubgraph( i, j + 1 );
						for( k = 1; k <= itMap->first[i][j]; k++ )
							biTemp *= biOrderTemp;
					}
				}
			}

			brAlternateTemp += MPZ_rational( itMap->second, biTemp );
			
			if( iDimension )
			{
				if( iFVectorIndex < 0 )
					return false;
				
				iFVector[ iFVectorIndex ] += itMap->second;
			}
			if( bDebug )
				cout << "N: " << itMap->second << " / Order: " << biTemp.get_str() << endl;
		}
		
		iCurrentVerticesCount++;
		
		if( bPositive )
			brEulerCaracteristic += brAlternateTemp;
		else
			brEulerCaracteristic -= brAlternateTemp;
		

		bPositive = !bPositive;
		iFVectorIndex--;
	}
	
	if( iDimension && graphsProductsCount_euclidean.size() + 1 > iDimension )
	{
		// si la dimension est spécifiée, on va mettre à jour le f-vecteur et la somme alternée avec le nombre de sommets à l'infini
		iVerticesAtInfinityCount = 0;
		for( itMap = graphsProductsCount_euclidean[ iDimension - 1 ].begin(); itMap != graphsProductsCount_euclidean[ iDimension - 1 ].end(); ++itMap )
			iVerticesAtInfinityCount += itMap->second;
		
		iFVector[0] += iVerticesAtInfinityCount;
		
		for( i = 0; i < iDimension; i++ )
			iFVectorAlternateSum += ( ( i % 2 ) ? -1 : 1 ) * iFVector[i];
	}
	
	return true;
}

// ##################################################################################################################################3
// Affichages

void CoxIter::printEuclideanGraphsProducts( vector< map<vector< vector< short unsigned int > >, unsigned int> >* graphsProductsCount )
{
	// variables de boucles
	size_t i, j, iMax;
	map<vector< vector< short unsigned int > >, unsigned int>::iterator itMap;
	
	// par taille de nombre de sommets
	for( vector< map<vector< vector< short unsigned int > >, unsigned int> >::iterator itMaps( graphsProductsCount->begin() ); itMaps != graphsProductsCount->end(); ++itMaps )
	{
		// on parcourt les produits pour la taille donnée
		for( itMap = itMaps->begin(); itMap != itMaps->end(); ++itMap )
		{
			cout << "\t";
			// pour chaque type de graphe
			for( i = 0; i < 8; i++ )
			{
				iMax = (itMap->first[i]).size();
				
				// pour chaque taille
				for( j = 0; j < iMax; j++ )
				{
					if( itMap->first[i][j] )
						cout << "T" << (char)(i + 65) << "_" << j << "^" << itMap->first[i][j] << " | ";
				}
			}
			
			cout << "N: " << itMap->second << endl;
		}
	}
}

void CoxIter::printCoxeterMatrix()
{
	cout << "Coxeter matrix" << endl;
	
	cout << "\tVertices: ";
	for( vector< string >::const_iterator it( map_vertices_indexToLabel.begin() ); it != map_vertices_indexToLabel.end(); ++it )
		cout << ( it == map_vertices_indexToLabel.begin() ? "" : ", " ) << *it;
	cout << endl;
	
	unsigned int i, j;
	for( i = 0; i < iVerticesCount; i++ )
	{
		cout << "\t";
		for( j = 0; j < iVerticesCount; j++ )
		{
			cout << ( j ? "," : "" ) << ( i == j ? 1 : ( iCoxeterMatrix[i][j] < 2 ? 0 : iCoxeterMatrix[i][j] ) );
		}
		cout << endl;
	}
}

void CoxIter::printGramMatrix()
{
	if( strOuputMathematicalFormat == "gap" )
		printGramMatrix_GAP();
	else if( strOuputMathematicalFormat == "latex" )
		printGramMatrix_LaTeX();
	else if( strOuputMathematicalFormat == "mathematica" )
		printGramMatrix_Mathematica();
	else if( strOuputMathematicalFormat == "pari" )
		printGramMatrix_PARI();
	else
		cout << "Gram matrix  \n\t" << get_strGramMatrix() << "\n" << endl;
	
	unsigned int i, j;
	for( i = 0; i < iVerticesCount; i++ )
	{
		for( j = 0; j < i; j++ )
		{
			if( iCoxeterMatrix[i][j] == 1 && strWeights.find( iLinearizationMatrix_index(j,i, iVerticesCount) ) == strWeights.end() )
			{
				cout << "l" << j<< "m" << i << ": weight of the dotted line between hyperplanes " << map_vertices_indexToLabel[j] << " and " << map_vertices_indexToLabel[i] << endl;
			}
		}
	}
	
	cout << endl;
}

void CoxIter::printGramMatrix_GAP()
{
	cout << "Gram matrix (GAP): \n\t" << get_strGramMatrix_GAP() << "\n" << endl;
}

void CoxIter::printGramMatrix_Mathematica()
{
	cout << "Gram matrix (Mathematica): \n\t" << get_strGramMatrix_Mathematica() << "\n" << endl;
}

void CoxIter::printGramMatrix_PARI()
{
	cout << "Gram matrix (PARI): \n\t" << get_strGramMatrix_PARI() << "\n" << endl;
}

void CoxIter::printGramMatrix_LaTeX()
{
	cout << "Gram matrix (LaTeX): \n\t" << get_strGramMatrix_LaTeX() << "\n" << endl;
}


void CoxIter::printEdgesVisitedMatrix()
{
	unsigned int i, j;
	cout << "Matrix of visited edges" << endl;
	
	for( i = 0; i < iVerticesCount; i++ )
	{
		for( j = 0; j < iVerticesCount; j++ )
			cout << ( bEdgesVisited[i][j] ? 1 : 0 )  << " ";
		cout << endl;
	}
}

bool CoxIter::bIsVertexValid( const string& strVertexLabel ) const
{
	return ( find( map_vertices_indexToLabel.begin(), map_vertices_indexToLabel.end(), strVertexLabel ) != map_vertices_indexToLabel.end() );
}

void CoxIter::map_vertices_labels_removeReference( const unsigned int& iIndex )
{
	if( iIndex > map_vertices_indexToLabel.size() )
		return;
	
	for( unsigned int i( iIndex + 1 ); i < iVerticesCount; i++ )
		map_vertices_labelToIndex[ map_vertices_indexToLabel[i] ]--;
	
	map_vertices_labelToIndex.erase( map_vertices_indexToLabel[ iIndex ] );
	map_vertices_indexToLabel.erase( map_vertices_indexToLabel.begin() + iIndex );
}

void CoxIter::map_vertices_labels_addReference( const string& strLabel )
{
	map_vertices_labelToIndex[ strLabel ] = map_vertices_indexToLabel.size();
	map_vertices_indexToLabel.push_back( strLabel );
}

// ##################################################################################################################################3
// Accesseurs
unsigned int CoxIter::get_iVertexIndex( const string& strVertexLabel ) const
{
	map< string, unsigned int >::const_iterator it( map_vertices_labelToIndex.find( strVertexLabel ) );
	
	if( it == map_vertices_labelToIndex.end() )
		throw( string( "INVALID_VERTEX_NAME" ) );
	
	return it->second;
}

string CoxIter::get_strVertexLabel( const unsigned int& iVertex ) const
{
	if( iVertex >= iVerticesCount )
		throw( 0 );
	
	return map_vertices_indexToLabel[ iVertex ];
}

vector< vector< unsigned int > > CoxIter::get_iCoxeterMatrix() const
{
	return iCoxeterMatrix;
}

std::map< unsigned int, string > CoxIter::get_strWeights() const
{
	return strWeights;
}

string CoxIter::get_strCoxeterMatrix() const
{
	string strCox;
	
	for( vector< vector< unsigned int > >::const_iterator itRow( iCoxeterMatrix.begin() ); itRow != iCoxeterMatrix.end(); ++itRow )
	{
		if( itRow != iCoxeterMatrix.begin() )
			strCox += "\n";
		
		for( vector< unsigned int >::const_iterator itCol( itRow->begin() ); itCol != itRow->end(); ++itCol )
			strCox += ( itCol == itRow->begin() ? "" : "," ) + to_string( *itCol );
	}
	
	return strCox;
}

vector< vector< string > > CoxIter::get_array_str_2_GramMatrix() const
{
	size_t i, j;
	vector< vector< string > > strGramMatrix( vector< vector< string > >( iVerticesCount, vector< string >( iVerticesCount, "" ) ) );
	
	for( i = 0; i < iVerticesCount; i++ )
	{
		for( j = 0; j <= i; j++ )
		{
			if( i == j )
				strGramMatrix[i][i] = "2" ;
			else if( iCoxeterMatrix[i][j] == 0 )
				strGramMatrix[i][j] = "-2" ;
			else if( iCoxeterMatrix[i][j] == 1 )
				strGramMatrix[i][j] = "2*l" + to_string( static_cast<long long>( min( i, j ) ) ) + "m" + to_string( static_cast<long long>( max( i, j ) ) );
			else
			{
				if( iCoxeterMatrix[i][j] == 2 )
					strGramMatrix[i][j] = "0";
				else if( iCoxeterMatrix[i][j] == 3 )
					strGramMatrix[i][j] = "-1";
				else if( iCoxeterMatrix[i][j] == 4 )
					strGramMatrix[i][j] = "-sqrt(2)";
				else if( iCoxeterMatrix[i][j] == 5 )
					strGramMatrix[i][j] = "-(1+sqrt(5))/2";
				else if( iCoxeterMatrix[i][j] == 6 )
					strGramMatrix[i][j] = "-sqrt(3)";
				else
					strGramMatrix[i][j] = "-2*cos(%pi/" + to_string( static_cast<long long>( iCoxeterMatrix[i][j] ) ) + ")";
			}
			
			strGramMatrix[j][i] = strGramMatrix[i][j];
		}
	}
	
	return strGramMatrix;
}

string CoxIter::get_strGramMatrix() const
{
	size_t i, j;
	string strGramMatrix( "" );
	
	for( i = 0; i < iVerticesCount; i++ )
	{
		strGramMatrix += ( i ? ", [" : "[ " );
		for( j = 0; j < iVerticesCount; j++ )
		{
			if( j > 0 )
				strGramMatrix += ", ";
			
			if( i == j )
				strGramMatrix += "1" ;
			else if( iCoxeterMatrix[i][j] == 0 )
				strGramMatrix += "-1" ;
			else if( iCoxeterMatrix[i][j] == 1 )
			{
				map< unsigned int, string >::const_iterator itF( strWeights.find( iLinearizationMatrix_index( min(i,j), max(i,j), iVerticesCount ) ) );
				
				if( itF != strWeights.end() )
					strGramMatrix += itF->second;
				else
					strGramMatrix += "l" + to_string( static_cast<long long>( min( i, j ) ) ) + "m" + to_string( static_cast<long long>( max( i, j ) ) );
			}
			else
			{
				if( iCoxeterMatrix[i][j] == 2 )
					strGramMatrix += "0";
				else if( iCoxeterMatrix[i][j] == 3 )
					strGramMatrix += "-1/2";
				else if( iCoxeterMatrix[i][j] == 4 )
					strGramMatrix += "-sqrt(2)/2";
				else if( iCoxeterMatrix[i][j] == 5 )
					strGramMatrix += "-(1+sqrt(5))/4";
				else if( iCoxeterMatrix[i][j] == 6 )
					strGramMatrix += "-sqrt(3)/2";
				else
					strGramMatrix += "-cos(%pi/" + to_string( static_cast<long long>( iCoxeterMatrix[i][j] ) ) + ")";
			}
		}
		strGramMatrix += "]";
	}
	
	return strGramMatrix;
}

string CoxIter::get_strGramMatrix_LaTeX() const
{
	size_t i, j;
	
	string strGramMatrix( "G = \\left( \\begin{array}{*{" + to_string( iVerticesCount ) + "}{c}}" );
	for( i = 0; i < iVerticesCount; i++ )
	{
		strGramMatrix += ( i ? "\\\\" : "" );
		for( j = 0; j < iVerticesCount; j++ )
		{
			if( j > 0 )
				strGramMatrix += " & ";
			
			if( i == j )
				strGramMatrix += "1" ;
			else if( iCoxeterMatrix[i][j] == 0 )
				strGramMatrix += "-1" ;
			else if( iCoxeterMatrix[i][j] == 1 )
			{
				map< unsigned int, string >::const_iterator itF( strWeights.find( iLinearizationMatrix_index( min(i,j), max(i,j), iVerticesCount ) ) );
				
				if( itF != strWeights.end() )
					strGramMatrix += itF->second;
				else
					strGramMatrix += "l" + to_string( static_cast<long long>( min( i, j ) ) ) + "m" + to_string( static_cast<long long>( max( i, j ) ) );
			}
			else
			{
				if( iCoxeterMatrix[i][j] == 2 )
					strGramMatrix += "0";
				else if( iCoxeterMatrix[i][j] == 3 )
					strGramMatrix += "-\\frac{1}{2}";
				else if( iCoxeterMatrix[i][j] == 4 )
					strGramMatrix += "-\\frac{\\sqrt 2}{2}";
				else if( iCoxeterMatrix[i][j] == 6 )
					strGramMatrix += "-\\frac{\\sqrt 3}{2}";
				else
					strGramMatrix += "-\\cos\\big(\\frac{\\pi}{" + to_string( static_cast<long long>( iCoxeterMatrix[i][j] ) ) + "}\\big)";
			}
		}
	}
	
	strGramMatrix += "\\end{array} \\right)";

	return strGramMatrix;
}

string CoxIter::get_strGramMatrix_Mathematica() const
{
	size_t i, j;
	
	string strGramMatrix( "G := {" );
	for( i = 0; i < iVerticesCount; i++ )
	{
		strGramMatrix += ( i ? ", {" : "{ " );
		for( j = 0; j < iVerticesCount; j++ )
		{
			if( j > 0 )
				strGramMatrix += ", ";
			
			if( i == j )
				strGramMatrix += "1" ;
			else if( iCoxeterMatrix[i][j] == 0 )
				strGramMatrix += "-1" ;
			else if( iCoxeterMatrix[i][j] == 1 )
			{
				map< unsigned int, string >::const_iterator itF( strWeights.find( iLinearizationMatrix_index( min(i,j), max(i,j), iVerticesCount ) ) );
				
				if( itF != strWeights.end() )
					strGramMatrix += itF->second;
				else
					strGramMatrix += "l" + to_string( static_cast<long long>( min( i, j ) ) ) + "m" + to_string( static_cast<long long>( max( i, j ) ) );
			}
			else
			{
				if( iCoxeterMatrix[i][j] == 2 )
					strGramMatrix += "0";
				else if( iCoxeterMatrix[i][j] == 3 )
					strGramMatrix += "-1/2";
				else if( iCoxeterMatrix[i][j] == 4 )
					strGramMatrix += "-Sqrt[2]/2";
				else if( iCoxeterMatrix[i][j] == 6 )
					strGramMatrix += "-Sqrt[3]/2";
				else
					strGramMatrix += "-Cos[Pi/" + to_string( static_cast<long long>( iCoxeterMatrix[i][j] ) ) + "]";
			}
		}
		strGramMatrix += "}";
	}
	
	strGramMatrix += "};";

	return strGramMatrix;
}

string CoxIter::get_strGramMatrix_PARI() const
{
	size_t i, j;
	string strGramMatrix( "G = [" );
	
	for( i = 0; i < iVerticesCount; i++ )
	{
		strGramMatrix += ( i ? "; " : " " );
		for( j = 0; j < iVerticesCount; j++ )
		{
			if( j > 0 )
				strGramMatrix += ", ";
			
			if( i == j )
				strGramMatrix += "1" ;
			else if( iCoxeterMatrix[i][j] == 0 )
				strGramMatrix += "-1" ;
			else if( iCoxeterMatrix[i][j] == 1 )
			{
				map< unsigned int, string >::const_iterator itF( strWeights.find( iLinearizationMatrix_index( min(i,j), max(i,j), iVerticesCount ) ) );
				
				if( itF != strWeights.end() )
					strGramMatrix += itF->second;
				else
					strGramMatrix += "l" + to_string( static_cast<long long>( min( i, j ) ) ) + "m" + to_string( static_cast<long long>( max( i, j ) ) );
			}
			else
			{
				if( iCoxeterMatrix[i][j] == 2 )
					strGramMatrix += "0";
				else if( iCoxeterMatrix[i][j] == 3 )
					strGramMatrix += "-1/2";
				else if( iCoxeterMatrix[i][j] == 4 )
					strGramMatrix += "-sqrt(2)/2";
				else if( iCoxeterMatrix[i][j] == 5 )
					strGramMatrix += "-(1+sqrt(5))/4";
				else if( iCoxeterMatrix[i][j] == 6 )
					strGramMatrix += "-sqrt(3)/2";
				else
					strGramMatrix += "-cos(Pi/" + to_string( static_cast<long long>( iCoxeterMatrix[i][j] ) ) + ")";
			}
		}
	}
	
	return ( strGramMatrix + "];" );
}

string CoxIter::get_strGramMatrix_GAP() const
{
	size_t i, j;
	string strGramMatrix( "G := [ [" );
	
	for( i = 0; i < iVerticesCount; i++ )
	{
		strGramMatrix += ( i ? "], [" : " " );
		for( j = 0; j < iVerticesCount; j++ )
		{
			if( j > 0 )
				strGramMatrix += ", ";
			
			if( i == j )
				strGramMatrix += "1" ;
			else if( iCoxeterMatrix[i][j] == 0 )
				strGramMatrix += "-1" ;
			else if( iCoxeterMatrix[i][j] == 1 )
			{
				map< unsigned int, string >::const_iterator itF( strWeights.find( iLinearizationMatrix_index( min(i,j), max(i,j), iVerticesCount ) ) );
				
				if( itF != strWeights.end() )
					strGramMatrix += itF->second;
				else
					strGramMatrix += "l" + to_string( static_cast<long long>( min( i, j ) ) ) + "m" + to_string( static_cast<long long>( max( i, j ) ) );
			}
			else
			{
				if( iCoxeterMatrix[i][j] == 2 )
					strGramMatrix += "0";
				else if( iCoxeterMatrix[i][j] == 3 )
					strGramMatrix += "-1/2";
				else if( iCoxeterMatrix[i][j] == 4 )
					strGramMatrix += "-Sqrt(2)/2";
				else if( iCoxeterMatrix[i][j] == 5 )
					strGramMatrix += "-(1+Sqrt(5))/4";
				else if( iCoxeterMatrix[i][j] == 6 )
					strGramMatrix += "-Sqrt(3)/2";
				else
					strGramMatrix += "-Cos(FLOAT.PI/" + to_string( static_cast<long long>( iCoxeterMatrix[i][j] ) ) + ")";
			}
		}
	}
	
	return ( strGramMatrix + "] ];" );
}

string CoxIter::get_strGramMatrixField() const
{
	return ( bGramMatrixField ? strGramMatrixField : "" );
}


MPZ_rational CoxIter::get_brEulerCaracteristic() const
{
	return brEulerCaracteristic;
}

string CoxIter::get_strEulerCaracteristic() const
{
	return brEulerCaracteristic.to_string();
}

string CoxIter::get_strEulerCharacteristic_computations() const
{
	return strEulerCharacteristic_computations;
}

int CoxIter::get_iFVectorAlternateSum() const
{
	return iFVectorAlternateSum;
}

bool CoxIter::get_bWriteInfo() const
{
	return bWriteInfo;
}

void CoxIter::set_bWriteInfo( const bool& bNewValue )
{
	bWriteInfo = bNewValue;
}

bool CoxIter::get_bDebug() const
{
	return bDebug;
}

vector< unsigned int > CoxIter::get_iFVector() const
{
	return iFVector;
}

unsigned int CoxIter::get_iVerticesAtInfinityCount() const
{
	return iVerticesAtInfinityCount;
}

int CoxIter::get_iIsCocompact()
{
	if( iIsCocompact == -2 && bCheckCocompactness )
		return iIsGraphCocompact();
	
	return iIsCocompact;
}

int CoxIter::get_iIsFiniteCovolume()
{
	if( iIsFiniteCovolume == -2 && bCheckCofiniteness )
		return isFiniteCovolume();
	
	return iIsFiniteCovolume;
}

int CoxIter::get_iIsArithmetic() const
{
	return iIsArithmetic;
}

unsigned int CoxIter::get_iDimension() const
{
	return iDimension;
}

bool CoxIter::get_bDimensionGuessed() const
{
	return bDimension_guessed;
}

string CoxIter::get_strError() const
{
	return strError;
}

unsigned int CoxIter::get_iVerticesCount() const
{
	return iVerticesCount;
}

bool CoxIter::get_bHasDottedLine() const
{
	return bHasDottedLine;
}

int CoxIter::get_iHasDottedLineWithoutWeight() const
{
	return iHasDottedLineWithoutWeight;
}

string CoxIter::get_str_map_vertices_indexToLabel( const size_t& i ) const
{
	return map_vertices_indexToLabel[i];
}

GraphsList* CoxIter::get_gl_graphsList_spherical() const
{
	return graphsList_spherical;
}

GraphsList* CoxIter::get_gl_graphsList_euclidean() const
{
	return graphsList_euclidean;
}

bool CoxIter::get_b_hasSphericalGraphsOfRank( const unsigned int& iRank ) const
{
	if( iRank > iVerticesCount )
		return false;
	
	return ( graphsProductsCount_spherical[iRank].size() != 0 );
}

bool CoxIter::get_b_hasEuclideanGraphsOfRank( const unsigned int& iRank ) const
{
	if( iRank > iVerticesCount )
		return false;

	return ( graphsProductsCount_euclidean[iRank - 1].size() != 0 );
}

void CoxIter::set_iIsArithmetic( const unsigned int& iArithmetic )
{
	iIsArithmetic = iArithmetic; // Hope that the value given in paramater is correct
}

void CoxIter::set_bCheckCocompactness(const bool& bValue)
{
	bCheckCocompactness = bValue;
}

void CoxIter::set_bCheckCofiniteness(const bool& bValue)
{
	bCheckCofiniteness = bValue;
}

void CoxIter::set_bDebug(const bool& bValue)
{
	bDebug = bValue;
}

void CoxIter::set_bUseOpenMP(const bool& bValue)
{
	#ifdef _COMPILE_WITH_OPENMP_
	bUseOpenMP = bValue;
	#endif
}

void CoxIter::set_sdtoutToFile( const string& strFilename )
{
	string strOutputCoutFilename( strFilename );
	outCout = new ofstream( strOutputCoutFilename.c_str() );
	
	if( outCout->is_open() )
	{
		sBufOld = cout.rdbuf( outCout->rdbuf() );
		bCoutFile = true;
	}
}


void CoxIter::set_iDimension(const unsigned int& iDimension_)
{
	iDimension = iDimension_;
}

void CoxIter::set_iCoxeterMatrix( const vector< vector< unsigned int > >& iMat )
{
	strWeights.clear();
	iVerticesCount = iMat.size();
	
	initializations();
	iCoxeterMatrix = iMat;
	
	bGraphExplored = false;
	bGraphsProductsComputed = false;
	bGrowthSeriesComputed = false;
	
	iEuclideanMaxRankFound = 0;
	iSphericalMaxRankFound = 0;
	bDimension_guessed = false;
	iHasDottedLineWithoutWeight = -1;
}

void CoxIter::set_strOuputMathematicalFormat(const string& strO)
{
	strOuputMathematicalFormat = strO;
}

void CoxIter::set_strVerticesToConsider( const vector< string >& strVerticesToConsider )
{
	strVertices = strVerticesToConsider;
	sort( strVertices.begin(), strVertices.end() );
	strVertices = vector< string >( strVertices.begin(), unique( strVertices.begin(), strVertices.end() ) );
	
}

void CoxIter::set_strVerticesToRemove(const vector< string >& strVerticesRemove_)
{
	strVerticesRemove = strVerticesRemove_;
	sort( strVerticesRemove.begin(), strVerticesRemove.end() );
	strVerticesRemove = vector< string >( strVerticesRemove.begin(), unique( strVerticesRemove.begin(), strVerticesRemove.end() ) );
}

