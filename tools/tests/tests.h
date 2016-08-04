#ifndef __TESTS_H__
#define __TESTS_H__

#include "../../coxiter.h"
#include "../../growthrate.h"
#include "../../signature.h"
#include "../../../tools/regexp.h"
#include "../../arithmeticity.h"
#include "../../../tools/numbers/number_rational.h"

#include <iostream>
#include <string>
#include <fstream>
#include <map>
#include <vector>
#include <omp.h>

using namespace std;

/*!
 * \file tests.h
 * \author Rafael Guglielmetti
*/

/*! \struct Test
 * \brief Contient les informations pour réaliser un test
*/
struct Test
{
	string szFile; ///< Chemin vers le fichier
	
	bool bTestEuler; ///< Si on teste la caractéristique d'Euler
	Number_rational brResult; ///< Caractéristique d'Euler théorique
	
	bool bTestFVector; ///< Si on teste le f-vecteur
	vector< unsigned int > iFVector; ///< f-vecteur théorique
	
	bool bTestCompacity; ///< True if we test the cocompacity
	bool bIsCompact; ///< True if the group is cocompact, false otherwise
	bool bIsFiniteVolume; ///< True if finite covolume
	bool bTestArithmeticity; ///< True if we test the arithmeticity
	bool bIsArithmetic; ///< True if arithmetic
	
	bool bTestGrowthSeries;
	vector< mpz_class > growthSeries_iPolynomialDenominator; ///< (i-1)th term contains the coefficient of x^i
	vector< unsigned int > growthSeries_iCyclotomicNumerator; ///< Contains a list oif cyclotomic polynomials
	string strGrowthRate; ///< Empty or the growth rate
};

/*! \class Tests
 * \brief Classe de tests; devrait être appelée sur une grande liste de graphes après chaque modification sur le programme
*/
class Tests
{
	private:
		bool szError; ///< Si une erreur s'est produite
		
		vector< Test > tests_EulerCarac; ///< Tableau des graphes à analyser avec ( nom du fichier, caractéristique théorique )
		
		bool bCoutFile; ///< true si sortie vers fichier
		ofstream *outCout; ///< flux de sortie vers fichier
		streambuf *sBufOld; ///< pour restaurer le cout
		
	public:
		~Tests( );
		
		/*!
		 * 	\fn readGraphsFile
		 * 	\brief Lit la liste des graphes à analyser
		 * 
		 * 	Format de chaque ligne: chemin_vers_le_graph caractéristique_théorique
		 * 	
		 * 	\param szInputFilename( string ) Fichier à lire
		 * 	\return True si succès
		 */
		bool readGraphsFile( string szInputFilename );
		
		/*!
		 * 	\fn runTests
		 * 	\brief Fait les calculs pour les graphes dans tests_EulerCarac
		 */
		void runTests( );
		
		vector< Test > get_tests( ) const;
};

#endif