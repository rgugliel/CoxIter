/*
Copyright (C) 2013, 2014, 2015, 2016
Rafael Guglielmetti, rafael.guglielmetti@unifr.ch
*/

/*
This file is part of CoxIter or VINO or the meta package CIVA.

These are free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

CIVA is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with CIVA. If not, see <http://www.gnu.org/licenses/>.
*/

#include "paripolynomials.h"

namespace PariPolynomials
{
	GEN vector2t_POL( const vector< int >& iCoefficients ) // [i] = x^i
	{
		unsigned int iSize( iCoefficients.size() );
		const long iVar = 0; // Polynomial variable
		
		GEN pol( cgetg(iSize+2, t_POL) ); // iSize coefficients, sign and variable number
		pol[1] = evalsigne(1) | evalvarn(iVar) | evallgef(iSize+2); // not equal to zero, variable, coefficients
		
		for( unsigned int i(0); i < iSize; i++ )
			pol[i+2] = (long)stoi( iCoefficients[i] );
		
		return pol;
	}
	
	GEN vector2t_POL( const vector< mpz_class >& iCoefficients ) // [i] = x^i
	{
		unsigned int iSize( iCoefficients.size() );
		const long iVar = 0; // Polynomial variable
		
		GEN pol( cgetg(iSize+2, t_POL) ); // iSize coefficients, sign and variable number
		pol[1] = evalsigne(1) | evalvarn(iVar) | evallgef(iSize+2); // not equal to zero, variable, coefficients
		
		for( unsigned int i(0); i < iSize; i++ )
		{
			if( !iCoefficients[i].fits_slong_p() )
				throw( string( "PariPolynomials::vector2t_POL : Component does not fit into long int."  ) );
			
			pol[i+2] = (long)stoi( iCoefficients[i].get_si() );
		}
		
		return pol;
	}

	vector< long int > t_POL2vector( GEN poly )
	{
		unsigned int iSize( lg(poly) );
		vector< long int > iV;
		
		for( unsigned int i(2); i < iSize; i++ )
			iV.push_back( gtolong(gel(poly,i)) );
		
		return iV;
	}
}