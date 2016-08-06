/*
Copyright (C) 2013, 2014, 2015, 2016
Rafael Guglielmetti, rafael.guglielmetti@unifr.ch
*/

/*
This file is part of CoxIter and AlVin.

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

/*!
 * \file math_tools.h
 * \author Rafael Guglielmetti
 * 
 * \brief Some mathematical functions
*/

#ifndef __MATH_TOOLS_H__
#define __MATH_TOOLS_H__

#include <string>
#include <vector>
#include <map>
#include <gmpxx.h>

using namespace std;

namespace MathTools
{
	extern unsigned int iSmallPrimeNumbers[1229];
	
	/*! \fn bIsPrime
	* 	\brief Test if a number is prime
	* 	Remark: this function could/should be optimized
	* 	\param iN( unsigned int ) Integer
	* 	\return True if the number is prime, false otherwise
	*/
	bool bIsPrime( unsigned int iN );

	/*! 	\fn iSQRT
	* 	\brief Compute the integer square root of a positive integer (the greatest integer less than or equal to the square root of the given number)
	* 
	* 	\param iN( unsigned ) Integer
	* 	\return The integer sqaure root
	*/
	template <typename Type>
	typename std::enable_if<std::is_unsigned<Type>::value, Type>::type iSQRT( Type iN )
	{
		Type place = (Type)1 << ( sizeof (Type) * 8 - 2 ); // calculated by precompiler = same runtime as: place = 0x40000000  
		while( place > iN )  
			place /= 4; // optimized by complier as place >>= 2  

		Type root = 0; 
		
		while( place )  
		{  
			if( iN >= root+place )  
			{  
				iN -= root+place;  
				root += place * 2;  
			} 
			
			root /= 2;  
			place /= 4;  
		}  
		
		return root;  
	}

	/*! 	\fn iSQRTsup
	* 	\brief Compute the sup integer square root of a positive integer
	* 
	* 	\param iN( unsigned ) Integer
	* 	\return ceil( sqrt( iN ) )
	*/
	template <typename Type>
	typename std::enable_if<std::is_unsigned<Type>::value, Type>::type iSQRTsup( Type iN )
	{
		if( iN < 2 ) // 0 or 1
			return iN;
		
		Type place = (Type)1 << ( sizeof (Type) * 8 - 2 ); // calculated by precompiler = same runtime as: place = 0x40000000  
		iN--;
		while( place > iN )  
			place /= 4; // optimized by complier as place >>= 2  

		Type root = 0; 
		
		while( place )  
		{  
			if( iN >= root+place )  
			{  
				iN -= root+place;  
				root += place * 2;  
			} 
			
			root /= 2;  
			place /= 4;  
		}  
		
		return (root+1);  
	}

	/*! 	\fn iListDivisors
	* 	\brief Return the list of the divisors of a (small) number
	* 	Rermark: this function could/should be optimized
	* 	\param iN( const Type& ) Integer
	* 	\param bNonTrivialOnly( const bool& ) If true, does not return 1 and iN
	* 	\return The list of divisors
	*/
	template <typename Type>
	vector< typename std::enable_if<std::is_unsigned<Type>::value, Type>::type > iListDivisors( const Type& iN, const bool& bNonTrivialOnly = false )
	{
		static vector< vector< Type > > iDivisors_ = vector< vector< Type > >( {vector<Type>({1}),vector<Type>({1,2}),vector<Type>({1,3}),vector<Type>({1,2,4}),vector<Type>({1,5}),vector<Type>({1,2,3,6}),vector<Type>({1,7}),vector<Type>({1,2,4,8}),vector<Type>({1,3,9}),vector<Type>({1,2,5,10}),vector<Type>({1,11}),vector<Type>({1,2,3,4,6,12}),vector<Type>({1,13}),vector<Type>({1,2,7,14}),vector<Type>({1,3,5,15}),vector<Type>({1,2,4,8,16}),vector<Type>({1,17}),vector<Type>({1,2,3,6,9,18}),vector<Type>({1,19}),vector<Type>({1,2,4,5,10,20}),vector<Type>({1,3,7,21}),vector<Type>({1,2,11,22}),vector<Type>({1,23}),vector<Type>({1,2,3,4,6,8,12,24}),vector<Type>({1,5,25}),vector<Type>({1,2,13,26}),vector<Type>({1,3,9,27}),vector<Type>({1,2,4,7,14,28}),vector<Type>({1,29}),vector<Type>({1,2,3,5,6,10,15,30}),vector<Type>({1,31}),vector<Type>({1,2,4,8,16,32}),vector<Type>({1,3,11,33}),vector<Type>({1,2,17,34}),vector<Type>({1,5,7,35}),vector<Type>({1,2,3,4,6,9,12,18,36}),vector<Type>({1,37}),vector<Type>({1,2,19,38}),vector<Type>({1,3,13,39}),vector<Type>({1,2,4,5,8,10,20,40}),vector<Type>({1,41}),vector<Type>({1,2,3,6,7,14,21,42}),vector<Type>({1,43}),vector<Type>({1,2,4,11,22,44}),vector<Type>({1,3,5,9,15,45}),vector<Type>({1,2,23,46}),vector<Type>({1,47}),vector<Type>({1,2,3,4,6,8,12,16,24,48}),vector<Type>({1,7,49}),vector<Type>({1,2,5,10,25,50}),vector<Type>({1,3,17,51}),vector<Type>({1,2,4,13,26,52}),vector<Type>({1,53}),vector<Type>({1,2,3,6,9,18,27,54}),vector<Type>({1,5,11,55}),vector<Type>({1,2,4,7,8,14,28,56}),vector<Type>({1,3,19,57}),vector<Type>({1,2,29,58}),vector<Type>({1,59}),vector<Type>({1,2,3,4,5,6,10,12,15,20,30,60})} );
		
		static vector< vector< Type > > iDivisors_bNonTrivialOnly = vector< vector< Type > >( {vector<Type>(0),vector<Type>(0),vector<Type>(0),vector<Type>({2}),vector<Type>(0),vector<Type>({2,3}),vector<Type>(0),vector<Type>({2,4}),vector<Type>({3}),vector<Type>({2,5}),vector<Type>(0),vector<Type>({2,3,4,6}),vector<Type>(0),vector<Type>({2,7}),vector<Type>({3,5}),vector<Type>({2,4,8}),vector<Type>(0),vector<Type>({2,3,6,9}),vector<Type>(0),vector<Type>({2,4,5,10}),vector<Type>({3,7}),vector<Type>({2,11}),vector<Type>(0),vector<Type>({2,3,4,6,8,12}),vector<Type>({5}),vector<Type>({2,13}),vector<Type>({3,9}),vector<Type>({2,4,7,14}),vector<Type>(0),vector<Type>({2,3,5,6,10,15}),vector<Type>(0),vector<Type>({2,4,8,16}),vector<Type>({3,11}),vector<Type>({2,17}),vector<Type>({5,7}),vector<Type>({2,3,4,6,9,12,18}),vector<Type>(0),vector<Type>({2,19}),vector<Type>({3,13}),vector<Type>({2,4,5,8,10,20}),vector<Type>(0),vector<Type>({2,3,6,7,14,21}),vector<Type>(0),vector<Type>({2,4,11,22}),vector<Type>({3,5,9,15}),vector<Type>({2,23}),vector<Type>(0),vector<Type>({2,3,4,6,8,12,16,24}),vector<Type>({7}),vector<Type>({2,5,10,25}),vector<Type>({3,17}),vector<Type>({2,4,13,26}),vector<Type>(0),vector<Type>({2,3,6,9,18,27}),vector<Type>({5,11}),vector<Type>({2,4,7,8,14,28}),vector<Type>({3,19}),vector<Type>({2,29}),vector<Type>(0),vector<Type>({2,3,4,5,6,10,12,15,20,30})} );
		
		if( iN <= 60 )
		{
			if( bNonTrivialOnly )
				return iDivisors_bNonTrivialOnly[ iN - 1 ];
			else
				return iDivisors_[ iN - 1 ];
		}
		
		vector< Type > iDivisors;
		
		if( !bNonTrivialOnly )
		{
			iDivisors.push_back( 1 );
			iDivisors.push_back( iN );
		}
		
		Type iMax( iSQRT( iN ) ), iTemp;
		for( Type i(2); i <= iMax; i++ )
		{
			if( !( iN % i ) )
			{
				iTemp = iN / i;
				iDivisors.push_back( i );
				if( i != iTemp )
					iDivisors.push_back( iTemp );
			}
		}
		
		return iDivisors;
	}

	/*! 	\fn ugcd
	* 	\brief Compute the gcd of two positive integers
	* 
	* 	\param u( integer ) First integer
	* 	\param v( inter ) Second integer
	* 	\return The gcd
	*/
	template <typename Type>
	typename std::enable_if<std::is_unsigned<Type>::value, Type>::type ugcd( Type u, Type v )
	{
		while ( v != 0) 
		{
			unsigned int r = u % v;
			u = v;
			v = r;
		}
		
		return u;
	}

	/*! 	\fn iJacobiSymbol
	* 	\brief Compute the Jacobi symbol of two integers
	* 	
	* 	\param iA( int ) Integer
	*	\param iB( int ) Integer
	* 	\return The symbol
	*/
	int iJacobiSymbol( int iA, unsigned int iB );

	/*! 	\fn iPrimeFactorsWithoutSquares
	* 	\brief Return the prime factors dividing the integer
	* 	
	* 	\param iN( unsigned int ) Integer
	* 	\return Prime factors
	*/
	vector< unsigned int > iPrimeFactorsWithoutSquares( unsigned int iN );

	/*! 	\fn iPrimeFactors
	* 	\brief Return the prime factors dividing the integer
	* 	
	* 	\param iN Integer
	* 	\return Prime factors
	*/
	template <typename Type>
	vector<Type> iPrimeFactors( const Type& iN )
	{
		Type iWork( iN < 0 ? -iN : iN );
		vector< Type > iPrimes;
	
		if( iWork == 1 )
			return vector< Type >( );
		
		if( !( iWork % 2 ) )
		{
			iPrimes.push_back(2);
			iWork /= 2;
			
			while( !( iWork % 2 ) )
				iWork/= 2;
		}
		
		for( unsigned int i(1); i < 1229 && iWork > 1; i++ )
		{
			if( !( iWork % iSmallPrimeNumbers[i] ) )
			{
				iPrimes.push_back( iSmallPrimeNumbers[i] );
				iWork/= iSmallPrimeNumbers[i];
				
				while( !( iWork % iSmallPrimeNumbers[i] ) )
					iWork/= iSmallPrimeNumbers[i];
			}
		}
		
		Type iDivisor( 10007 );
		while( iWork > 1 )
		{
			if( !( iWork % iDivisor ) )
			{
				iPrimes.push_back( iDivisor );
				iWork/= iDivisor;
				
				while( !( iWork % iDivisor ) )
					iWork/= iDivisor;
			}
			
			iDivisor += 2;
		}
		
		return iPrimes;
	}

	/*! 	\fn iPrimeDecomposition
	* 	\brief Return the prime decomposition of thee integer
	* 	
	* 	\param iN( unsigned int ) Integer
	* 	\return Prime decomposition: [ prime ] = power
	*/
	template <typename Type>
	typename std::enable_if<std::is_unsigned<Type>::value, map< Type, unsigned int >>::type iPrimeDecomposition( Type iN )
	{
		map< Type, unsigned int > iPrimes;
	
		if( iN == 1 )
			return map< Type, unsigned int >( );
		
		if( !( iN % 2 ) )
		{
			iPrimes[2] = 1;
			iN /= 2;
			
			while( !( iN % 2 ) )
			{
				iN/= 2;
				iPrimes[2]++;
			}
		}
		
		for( unsigned int i(1); i < 1229 && iN > 1; i++ )
		{
			if( !( iN % iSmallPrimeNumbers[i] ) )
			{
				iPrimes[iSmallPrimeNumbers[i]] = 1;
				iN/= iSmallPrimeNumbers[i];
				
				while( !( iN % iSmallPrimeNumbers[i] ) )
				{
					iN/= iSmallPrimeNumbers[i];
					iPrimes[iSmallPrimeNumbers[i]]++;
				}
			}
		}
		
		Type iDivisor( 10007 );
		while( iN > 1 )
		{
			if( !( iN % iDivisor ) )
			{
				iPrimes[iDivisor] = 1;
				iN/= iDivisor;
				
				while( !( iN % iDivisor ) )
				{
					iN/= iDivisor;
					iPrimes[iDivisor]++;
				}
			}
			
			iDivisor += 2;
		}
		
		return iPrimes;
	}

	/*! 	\fn iRemoveSquareFactors
	* 	\brief Remove square factors of a positive integer
	* 	
	* 	\param iN( unsigned int ) Integer
	* 	\return iN divided by all its squared factors
	*/
	template <typename Type>
	Type iRemoveSquareFactors( Type iN )
	{
		Type iWork( iN > 0 ? iN : -iN ), iRes( 1 ), iSquare;
		
		if( iWork < 3 )
			return iN;
		
		for( unsigned int i(0); i < 1229 && iWork > 1; i++ )
		{
			iSquare = iSmallPrimeNumbers[i] * iSmallPrimeNumbers[i];
			
			while( iWork % iSquare == 0 )
				iWork /= iSquare;
			
			if( iWork % iSmallPrimeNumbers[i] == 0 )
			{
				iWork /= iSmallPrimeNumbers[i];
				iRes *= iSmallPrimeNumbers[i];
			}
		}
		
		Type iDivisor( 10007 );
		while( iWork > 1 )
		{
			iSquare = iDivisor * iDivisor;
			
			while( iWork % iSquare == 0 )
				iWork /= iSquare;
			
			if( iWork % iDivisor == 0 )
			{
				iWork /= iDivisor;
				iRes *= iDivisor;
			}
				
			iDivisor += 2;
		}
		
		return ( iN < 0 ? -iRes : iRes );
	}

	/*! 	\fn iCeilQuotient
	* 	\brief Compute the ceiling of a quotient
	* 	
	* 	\param iNumerator( const unsigned int& ) Numerator
	* 	\param iDenominator( const unsigned int& ) Denominator
	* 	\return ceil( iNumerator / iDenominator )
	*/
	template <typename Type>
	inline Type iCeilQuotient( const Type& iNumerator, const Type& iDenominator )
	{
		return ( (iNumerator % iDenominator) ? iNumerator / iDenominator + 1 : iNumerator / iDenominator );
	}

	/*! 	\fn iSQRTQuotient
	* 	\brief Compute the greatest integer less than or equal to the square root of the given rational number
	* 	
	* 	\param iNumerator( const unsigned & ) Numerator
	* 	\param iDenominator( const unsigned & ) Denominator
	* 	\return ceil( sqrt( iNumerator / iDenominator ) )
	*/
	template <typename Type>
	typename std::enable_if<std::is_unsigned<Type>::value, Type>::type iSQRTQuotient( const Type& iNumerator, const Type& iDenominator )
	{
		Type tRes( iSQRT( iNumerator / iDenominator ) );
		
		while( iDenominator * ( tRes + 1 ) * ( tRes + 1 ) <= iNumerator )
			tRes++;
		
		return tRes;
	}

	mpz_class iSQRTQuotient( const mpz_class& iNumerator, const mpz_class& iDenominator );
	mpz_class iSQRTsupQuotient( const mpz_class& iNumerator, const mpz_class& iDenominator );

	/*! 	\fn iSQRTsupQuotient
	* 	\brief Compute the smalles integer greater than or equal to the square root of the given rational number
	* 	
	* 	\param iNumerator( const unsigned & ) Numerator
	* 	\param iDenominator( const unsigned & ) Denominator
	* 	\return ceil( sqrt( iNumerator / iDenominator ) )
	*/
	template <typename Type>
	typename std::enable_if<std::is_unsigned<Type>::value, Type>::type iSQRTsupQuotient( const Type& iNumerator, const Type& iDenominator )
	{
		if( !iNumerator )
			return 0;
		
		Type tRes( (iNumerator % iDenominator) != 0 ? iNumerator / iDenominator + 1 : iNumerator / iDenominator );
		tRes = iSQRTsup( tRes );

		while( iNumerator <= iDenominator * ( tRes - 1 ) * ( tRes - 1 ) )
			tRes--;
		
		return tRes;
	}
}
#endif
