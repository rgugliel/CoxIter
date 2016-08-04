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

#include "growthrate.h"

GrowthRate::GrowthRate()
{
	pari_init( 50000000, 2 );
	
	/* Note: gEpsilon must be BIG compared to the precision up to which we compute the roots (so that
	 * we can detect when a root is too small).
	 * */
	gEpsilon = dbltor( 1e-50 );
	iPariPrecision = 8;
	
	// The following two lines are useless except with the test program.
	
	long prec;
	setrealprecision( 38, &prec );
	
}

GrowthRate::~GrowthRate()
{
	pari_close();
}

// TODO: résultat Salem
GrowthRate_Result GrowthRate::grrComputations( vector< mpz_class > iPolynomial, const bool& bOnlyGrowthRate )
{
	irreducibleFactors( iPolynomial );
	gGrowthRate = dbltor( 1.0 ); // Have to allocate the memory for this outside of some functions
	minimalRoot();

	// ----------------------------------------------------------
	pari_sp av = avma; // current state of the PARI stack
	
	// ----------------------------------------------------
	// Results
	GrowthRate_Result grr;
	grr.iPerron = 1;
	grr.iPisot = -1; // Not decided yed
	grr.iSalem = -1; // Not decided yed
	grr.strGrowthRate = GENtostr( gGrowthRate );
	
	// ----------------------------------------------------
	// Here we know the minimal root and the corresponding polynomial
	GEN gTemp, gRoot; // Temp variable
	
	if( bOnlyGrowthRate )
	{
		avma = av; // cleaning the stack
		grr.bComputed = true;
		return grr;
	}
	
	// ----------------------------------------------------
	// Control variables
	long int iRealRootsCount( sturm( t_POLfactors[ iIndexMaximalRoot ] ) ), iRealRootsFound( 0 ); // Number of real roots and the one detected as such
	long int iRealNegativeRootsCount( sturmpart( t_POLfactors[ iIndexMaximalRoot ], NULL, gen_0 ) ), iRealNegativeRootsFound( 0 ); // Number of negative real roots and the one detected as such
	long int iRootsOnUnitCircle( iNumberRootsUnitCircle( t_POLfactors[ iIndexMaximalRoot ] ) ), iRootsOnUnitCircleFound( 0 );
	long int iNumberRootsInsideUnitCircle( 0 ), iNumberRootsMaybeInsideUnitCircle( 0 );
	
	if( iRootsOnUnitCircle == -1 || degree( t_POLfactors[ iIndexMaximalRoot ] ) % 2 ) // The polynomial is not palindromic or the degree is odd
		grr.iSalem = 0; // Not Salem
	else
		grr.iSalem = degree( t_POLfactors[ iIndexMaximalRoot ] ) - 2 == iRootsOnUnitCircle;
	
	unsigned long iRootsCount( lg( gMaximalRoots ) );
	
	for( unsigned long int i(1); i < iRootsCount; i++ )
	{
		gRoot = gel(gMaximalRoots, i);

		// ----------------------------------------
		// Comparison of |gRoot| with 1
		gTemp = cxnorm( gRoot ); // gTemp := |gRoot|^2
		
		// On the unit circle
		if( mpcmp( gEpsilon, absr( gsub( gTemp, gen_1 ) ) ) > 0 ) // Probably = 1
			iRootsOnUnitCircleFound++;
		
		// Inside the unit circle
		if( mpcmp( gadd( gTemp, gEpsilon ), gen_1 ) < 0 )
			iNumberRootsInsideUnitCircle++;
		if( mpcmp( gsub( gTemp, gEpsilon ), gen_1 ) < 0 )
			iNumberRootsMaybeInsideUnitCircle++;
		
		if( mpcmp( absr( gimag( gRoot ) ), gEpsilon ) > 0 ) // The element is non-real
			continue;
		else
		{
			// ----------------------------------------
			// Perron?
			if( mpcmp( greal( gRoot ), gen_0 ) < 0 ) // If negative
			{
				iRealNegativeRootsFound++;
				
				// We check that |gRoot| + epsilon < gGrowRate
				if( mpcmp( gadd( greal( gRoot ), gEpsilon ), gGrowthRate ) >= 0 && grr.iPerron != -1 )
				{
					if( mpcmp( gsub( greal( gRoot ), gEpsilon ), gGrowthRate ) < 0 )
						grr.iPerron = -1; // We cannot decide
					else
					{
						cout << "tau1="; output( gGrowthRate ); cout << endl;
						cout << "Peron0;\nroot=";
						output( gRoot );
						cout << "\ntau2="; output( gGrowthRate ); 
						cout << endl;
						grr.iPerron = 0; // Not Perron
					}
				}
			}
			
			iRealRootsFound++;
		}
	}
		
	if( iRealRootsCount != iRealRootsFound )
		throw( string( "GrowthRate::grrComputations: Some root have a too small imaginary part" ) );
	
	if( iRealNegativeRootsCount != iRealNegativeRootsFound )
		throw( string( "GrowthRate::grrComputations: Some negative root is too close to zero" ) );
	
	if( iNumberRootsInsideUnitCircle != iNumberRootsMaybeInsideUnitCircle )
	{
		/*
		 * If we were able to determine the number of roots on the unit circle and if these are spare roots,
		 * then we can remove them.
		*/
		
		if( iRootsOnUnitCircle == -1 || iRootsOnUnitCircleFound != iRootsOnUnitCircle || ( iNumberRootsMaybeInsideUnitCircle - iNumberRootsInsideUnitCircle ) != iRootsOnUnitCircle )
			grr.iPisot = -2; // Cannot decide
		else
			iNumberRootsInsideUnitCircle -= iRootsOnUnitCircleFound;
	}
	
	if( grr.iPisot == -1 )
		grr.iPisot = iNumberRootsInsideUnitCircle == ( degree( t_POLfactors[ iIndexMaximalRoot ] ) - 1 ) ? 1 : 0;
	
	// ----------------------------------------------------------
	avma = av; // cleaning the stack
	grr.bComputed = true;
	return grr;
}

void GrowthRate::irreducibleFactors( const vector< mpz_class >& iPolynomial )
{
	// ---------------------------------------------------
	// Factors
	GEN gDenominator, gFactors;
	long int iRCount;
	
	gDenominator = vector2t_POL( iPolynomial ); // Conversion vector< mpz_class > to t_POL
	gFactors = factor( gDenominator ); // Factorization of the polynomial
	gFactors = gel(gFactors,1); // Irreducible factors
	long k( lg(gFactors) - 1 ); // Number of factors
	
	for( unsigned int i(1); i <= k; i++ )
	{
		vector< long int > vec( t_POL2vector( gel(gFactors,i) ) );
		
		if( ( iRCount = sturmpart( gel(gFactors,i), gen_0, gen_1 ) ) ) // If there exists at leat one real root between 0 and 1
		{
			// We want to excluce the case where the only roots in ]0,1] is 1
			if( poleval( gel(gFactors,i), gen_1 ) != gen_0 || iRCount > 1 )
			{
				t_POLfactors.push_back( RgX_recip( gel(gFactors,i) ) );
			}
		}
	}
}

void GrowthRate::minimalRoot()
{
	pari_sp ltop = avma; // current state of the PARI stack
	
	GEN gTemp;
	
	long int iRootsBetween1AndInfinity; // Number of roots between 1 and infinity
	long int iFoundRoots; // Number of roots we wound (to be compared with iRootsBetween0And1)
	
	unsigned int iFactorsCount( t_POLfactors.size() );
	
	for( unsigned int j(0); j < iFactorsCount; j++ )
	{
		auto f( t_POLfactors[j] );
		
		// ------------------------------------------------
		// Number of roots we have to find
		iFoundRoots = 0;
		iRootsBetween1AndInfinity = sturmpart( f, gen_1, NULL );

		GEN gRoots( roots(f, iPariPrecision) );
		
		unsigned long iRootsCount( lg( gRoots ) );
		for( unsigned long int i(1); i < iRootsCount; i++ )
		{
			gTemp = gel(gRoots, i);
			
			if( mpcmp( absr( gimag( gTemp ) ), gEpsilon ) > 0 ) // The element is non-real
				continue;
			
			if( mpcmp( absr( greal( gTemp ) ), gen_1 ) <= 0 ) // If the root is <= 1 in absolute value
				continue;
			
			if( mpcmp( greal( gTemp ), gen_0 ) <= 0 ) // If the root is <= 0
				continue;
			iFoundRoots++;
			
			// If some roots are to close
			if( mpcmp( gEpsilon, absr( gsub( gGrowthRate, gTemp ) ) ) < 0 )
			{
				avma = ltop; // cleaning the stack
				throw( string( "GrowthRate::minimalRoot() : Some real roots are to close" ) );
			}
			
			// Greather than the greatest?
			if( mpcmp( gGrowthRate, greal( gTemp ) ) < 0 )
			{
				iIndexMaximalRoot = j;
				gGrowthRate = greal( gTemp );
				gMaximalRoots = gRoots;
			}
		}
		
		if( iFoundRoots != iRootsBetween1AndInfinity )
		{
			avma = ltop; // cleaning the stack
			throw( string( "GrowthRate::minimalRoot() : Number of roots" ) );
		}
	}
	
	gerepileall( ltop, 2, &gMaximalRoots, &gGrowthRate );
}

long int GrowthRate::iNumberRootsUnitCircle( GEN gPol )
{
	pari_sp ltop = avma;
	if( cmp_RgX( gPol, RgX_recip( gPol ) ) != 0 ) // If the polynomial is not palindromic
	{
		avma = ltop;
		return -1;
	}
	
	long iDegree( degree( gPol ) );
	long iMax( iDegree % 2 ? ( iDegree - 1)/2 : iDegree/2-1 );
	
	// -------------------------------------------
	// some useful polynomials
	GEN gPolx( mkpoln( 2, gen_1, gen_0 ) ); // x
	GEN gPol101( mkpoln( 3, gen_1, gen_0, gen_1 ) ); // x^2 + 1
	GEN gPolResult( mkpoln( 1, gen_0 ) );
	
	if( !( iDegree % 2 ) )
		gPolResult = gmul( gel(gPol,iDegree/2+2), powgi( gPol101, stoi( iDegree/2 ) ) );
	
	for( long int j(0); j <= iMax; j++ )
	{
		GEN gPolTemp( mkpoln( 1, gen_0 ) );
		long int ikMax( iDegree-2*j );
		GEN gkMax( stoi( iDegree-2*j ) );
		
		for( long int k( iDegree % 2 ); k <= ikMax; k += 2 )
		{
			if( (ikMax - k) % 4 == 2 )
				gPolTemp = gadd( gPolTemp, gmul( gneg( gen_1 ), gmul( binomial(gkMax, k), powgi( gPolx, stoi(k) ) ) ) );
			else
				gPolTemp = gadd( gPolTemp, gmul( binomial(gkMax, k ), powgi( gPolx, stoi( k ) ) ) );
		}
		
		gPolResult = gadd( gPolResult, gmul( gel(gPol,j+2), gmul( powgi( gPol101, stoi(j) ), gmul( gPolTemp, gen_2 ) ) ) );
	}
	
	pari_sp lbot = avma;
	gPolResult = gerepile( ltop, lbot, gcopy(gPolResult) ); // we free the memory (except for gPolResult)
	
	long iResult( sturm( gPolResult )  );
	
	avma = ltop; // We free the memory
	
	return iResult;
}