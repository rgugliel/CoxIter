/*
Copyright (C) 2013-2017
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

#include <map>
#include <string>
#include <vector>
#ifdef _USE_LOCAL_GMP_
#include "gmpxx.h"
#else
#include <gmpxx.h>
#endif

using namespace std;

namespace MathTools {
extern unsigned int iSmallPrimeNumbers[1229];

/*! \fn isPrime
 * 	\brief Test if a number is prime
 * 	Remark: this function could/should be optimized
 * 	\param n(unsigned int) Integer
 * 	\return True if the number is prime, false otherwise
 */
bool isPrime(unsigned int n);

/*! 	\fn integerSqrt
 * 	\brief Compute the integer square root of a positive integer (the
 * greatest integer less than or equal to the square root of the given number)
 *
 * 	\param iN(unsigned) Integer
 * 	\return The integer sqaure root
 */
template <typename Type>
typename std::enable_if<std::is_unsigned<Type>::value, Type>::type
integerSqrt(Type n) {
  Type place =
      (Type)1
      << (sizeof(Type) * 8 -
          2); // calculated by precompiler = same runtime as: place = 0x40000000
  while (place > n)
    place /= 4; // optimized by complier as place >>= 2

  Type root = 0;

  while (place) {
    if (n >= root + place) {
      n -= root + place;
      root += place * 2;
    }

    root /= 2;
    place /= 4;
  }

  return root;
}

/*! 	\fn sqrtSup
 * 	\brief Compute the sup integer square root of a positive integer
 *
 * 	\param n(unsigned) Integer
 * 	\return ceil(sqrt(n))
 */
template <typename Type>
typename std::enable_if<std::is_unsigned<Type>::value, Type>::type
sqrtSup(Type n) {
  if (n < 2) // 0 or 1
    return n;

  Type place =
      (Type)1
      << (sizeof(Type) * 8 -
          2); // calculated by precompiler = same runtime as: place = 0x40000000
  n--;
  while (place > n)
    place /= 4; // optimized by complier as place >>= 2

  Type root = 0;

  while (place) {
    if (n >= root + place) {
      n -= root + place;
      root += place * 2;
    }

    root /= 2;
    place /= 4;
  }

  return (root + 1);
}

/*! 	\fn iListDivisors
 * 	\brief Return the list of the divisors of a (small) number
 * 	Rermark: this function could/should be optimized
 * 	\param iN(const Type&) Integer
 * 	\param nonTrivialOnly(const bool&) If true, does not return 1 and iN
 * 	\return The list of divisors
 */
template <typename Type>
vector<typename std::enable_if<std::is_unsigned<Type>::value, Type>::type>
listDivisors(const Type &n, const bool &nonTrivialOnly = false) {
#ifdef _MSC_VER
  static vector<vector<Type>> divisors_(
      vector<vector<Type>>(60, vector<Type>(0)));
  static vector<vector<Type>> divisors_nonTrivialOnly_(
      vector<vector<Type>>(60, vector<Type>(0)));

  if (n <= 60) {
    if (nonTrivialOnly && divisors_nonTrivialOnly_[n - 1].size())
      return divisors_nonTrivialOnly_[n - 1];

    if (!nonTrivialOnly && divisors_[n - 1].size())
      return divisors_[n - 1];
  }
#else // Remark: the following two initialisers don't work on Visual C++ 2013
      // and 2015
  static vector<vector<Type>> divisors_ = vector<vector<Type>>(
      {vector<Type>({1}),
       vector<Type>({1, 2}),
       vector<Type>({1, 3}),
       vector<Type>({1, 2, 4}),
       vector<Type>({1, 5}),
       vector<Type>({1, 2, 3, 6}),
       vector<Type>({1, 7}),
       vector<Type>({1, 2, 4, 8}),
       vector<Type>({1, 3, 9}),
       vector<Type>({1, 2, 5, 10}),
       vector<Type>({1, 11}),
       vector<Type>({1, 2, 3, 4, 6, 12}),
       vector<Type>({1, 13}),
       vector<Type>({1, 2, 7, 14}),
       vector<Type>({1, 3, 5, 15}),
       vector<Type>({1, 2, 4, 8, 16}),
       vector<Type>({1, 17}),
       vector<Type>({1, 2, 3, 6, 9, 18}),
       vector<Type>({1, 19}),
       vector<Type>({1, 2, 4, 5, 10, 20}),
       vector<Type>({1, 3, 7, 21}),
       vector<Type>({1, 2, 11, 22}),
       vector<Type>({1, 23}),
       vector<Type>({1, 2, 3, 4, 6, 8, 12, 24}),
       vector<Type>({1, 5, 25}),
       vector<Type>({1, 2, 13, 26}),
       vector<Type>({1, 3, 9, 27}),
       vector<Type>({1, 2, 4, 7, 14, 28}),
       vector<Type>({1, 29}),
       vector<Type>({1, 2, 3, 5, 6, 10, 15, 30}),
       vector<Type>({1, 31}),
       vector<Type>({1, 2, 4, 8, 16, 32}),
       vector<Type>({1, 3, 11, 33}),
       vector<Type>({1, 2, 17, 34}),
       vector<Type>({1, 5, 7, 35}),
       vector<Type>({1, 2, 3, 4, 6, 9, 12, 18, 36}),
       vector<Type>({1, 37}),
       vector<Type>({1, 2, 19, 38}),
       vector<Type>({1, 3, 13, 39}),
       vector<Type>({1, 2, 4, 5, 8, 10, 20, 40}),
       vector<Type>({1, 41}),
       vector<Type>({1, 2, 3, 6, 7, 14, 21, 42}),
       vector<Type>({1, 43}),
       vector<Type>({1, 2, 4, 11, 22, 44}),
       vector<Type>({1, 3, 5, 9, 15, 45}),
       vector<Type>({1, 2, 23, 46}),
       vector<Type>({1, 47}),
       vector<Type>({1, 2, 3, 4, 6, 8, 12, 16, 24, 48}),
       vector<Type>({1, 7, 49}),
       vector<Type>({1, 2, 5, 10, 25, 50}),
       vector<Type>({1, 3, 17, 51}),
       vector<Type>({1, 2, 4, 13, 26, 52}),
       vector<Type>({1, 53}),
       vector<Type>({1, 2, 3, 6, 9, 18, 27, 54}),
       vector<Type>({1, 5, 11, 55}),
       vector<Type>({1, 2, 4, 7, 8, 14, 28, 56}),
       vector<Type>({1, 3, 19, 57}),
       vector<Type>({1, 2, 29, 58}),
       vector<Type>({1, 59}),
       vector<Type>({1, 2, 3, 4, 5, 6, 10, 12, 15, 20, 30, 60})});
  static vector<vector<Type>> divisors_nonTrivialOnly_ =
      vector<vector<Type>>({vector<Type>(0),
                            vector<Type>(0),
                            vector<Type>(0),
                            vector<Type>({2}),
                            vector<Type>(0),
                            vector<Type>({2, 3}),
                            vector<Type>(0),
                            vector<Type>({2, 4}),
                            vector<Type>({3}),
                            vector<Type>({2, 5}),
                            vector<Type>(0),
                            vector<Type>({2, 3, 4, 6}),
                            vector<Type>(0),
                            vector<Type>({2, 7}),
                            vector<Type>({3, 5}),
                            vector<Type>({2, 4, 8}),
                            vector<Type>(0),
                            vector<Type>({2, 3, 6, 9}),
                            vector<Type>(0),
                            vector<Type>({2, 4, 5, 10}),
                            vector<Type>({3, 7}),
                            vector<Type>({2, 11}),
                            vector<Type>(0),
                            vector<Type>({2, 3, 4, 6, 8, 12}),
                            vector<Type>({5}),
                            vector<Type>({2, 13}),
                            vector<Type>({3, 9}),
                            vector<Type>({2, 4, 7, 14}),
                            vector<Type>(0),
                            vector<Type>({2, 3, 5, 6, 10, 15}),
                            vector<Type>(0),
                            vector<Type>({2, 4, 8, 16}),
                            vector<Type>({3, 11}),
                            vector<Type>({2, 17}),
                            vector<Type>({5, 7}),
                            vector<Type>({2, 3, 4, 6, 9, 12, 18}),
                            vector<Type>(0),
                            vector<Type>({2, 19}),
                            vector<Type>({3, 13}),
                            vector<Type>({2, 4, 5, 8, 10, 20}),
                            vector<Type>(0),
                            vector<Type>({2, 3, 6, 7, 14, 21}),
                            vector<Type>(0),
                            vector<Type>({2, 4, 11, 22}),
                            vector<Type>({3, 5, 9, 15}),
                            vector<Type>({2, 23}),
                            vector<Type>(0),
                            vector<Type>({2, 3, 4, 6, 8, 12, 16, 24}),
                            vector<Type>({7}),
                            vector<Type>({2, 5, 10, 25}),
                            vector<Type>({3, 17}),
                            vector<Type>({2, 4, 13, 26}),
                            vector<Type>(0),
                            vector<Type>({2, 3, 6, 9, 18, 27}),
                            vector<Type>({5, 11}),
                            vector<Type>({2, 4, 7, 8, 14, 28}),
                            vector<Type>({3, 19}),
                            vector<Type>({2, 29}),
                            vector<Type>(0),
                            vector<Type>({2, 3, 4, 5, 6, 10, 12, 15, 20, 30})});

  if (n <= 60) {
    if (nonTrivialOnly)
      return divisors_nonTrivialOnly_[n - 1];
    else
      return divisors_[n - 1];
  }
#endif

  vector<Type> divisors;

  if (!nonTrivialOnly) {
    divisors.push_back(1);
    divisors.push_back(n);
  }

  Type max(integerSqrt(n)), temp;
  for (Type i(2); i <= max; i++) {
    if (!(n % i)) {
      temp = n / i;
      divisors.push_back(i);
      if (i != temp)
        divisors.push_back(temp);
    }
  }

#ifdef _MSC_VER
  if (nonTrivialOnly)
    divisors_nonTrivialOnly_[n - 1] = iDivisors;
  else
    divisors_[n - 1] = iDivisors;
#endif

  return divisors;
}

/*! 	\fn ugcd
 * 	\brief Compute the gcd of two positive integers
 *
 * 	\param u(integer) First integer
 * 	\param v(inter) Second integer
 * 	\return The gcd
 */
template <typename Type>
typename std::enable_if<std::is_unsigned<Type>::value, Type>::type
ugcd(Type u, Type v) {
  while (v != 0) {
    unsigned int r = u % v;
    u = v;
    v = r;
  }

  return u;
}

/*! 	\fn jacobiSymbol
 * 	\brief Compute the Jacobi symbol of two integers
 *
 * 	\param a(int) Integer
 *	\param b(int) Integer
 * 	\return The symbol
 */
int jacobiSymbol(int a, unsigned int b);

/*! 	\fn primeFactorsWithoutSquares
 * 	\brief Return the prime factors dividing the integer
 *
 * 	\param n(unsigned int) Integer
 * 	\return Prime factors
 */
vector<unsigned int> primeFactorsWithoutSquares(unsigned int n);

/*! 	\fn primeFactors
 * 	\brief Return the prime factors dividing the integer
 *
 * 	\param n Integer
 * 	\return Prime factors
 */
template <typename Type> vector<Type> primeFactors(const Type &n) {
  Type nWork(n < 0 ? -n : n);
  vector<Type> primes;

  if (nWork == 1)
    return vector<Type>();

  if (!(nWork % 2)) {
    primes.push_back(2);
    nWork /= 2;

    while (!(nWork % 2))
      nWork /= 2;
  }

  for (unsigned int i(1); i < 1229 && nWork > 1; i++) {
    if (!(nWork % iSmallPrimeNumbers[i])) {
      primes.push_back(iSmallPrimeNumbers[i]);
      nWork /= iSmallPrimeNumbers[i];

      while (!(nWork % iSmallPrimeNumbers[i]))
        nWork /= iSmallPrimeNumbers[i];
    }
  }

  Type iDivisor(10007);
  while (nWork > 1) {
    if (!(nWork % iDivisor)) {
      primes.push_back(iDivisor);
      nWork /= iDivisor;

      while (!(nWork % iDivisor))
        nWork /= iDivisor;
    }

    iDivisor += 2;
  }

  return primes;
}

/*! 	\fn primeDecomposition
 * 	\brief Return the prime decomposition of thee integer
 *
 * 	\param n(unsigned int) Integer
 * 	\return Prime decomposition: [ prime ] = power
 */
template <typename Type>
typename std::enable_if<std::is_unsigned<Type>::value,
                        map<Type, unsigned int>>::type
primeDecomposition(Type n) {
  map<Type, unsigned int> primes;

  if (n == 1)
    return map<Type, unsigned int>();

  if (!(n % 2)) {
    primes[2] = 1;
    n /= 2;

    while (!(n % 2)) {
      n /= 2;
      primes[2]++;
    }
  }

  for (unsigned int i(1); i < 1229 && n > 1; i++) {
    if (!(n % iSmallPrimeNumbers[i])) {
      primes[iSmallPrimeNumbers[i]] = 1;
      n /= iSmallPrimeNumbers[i];

      while (!(n % iSmallPrimeNumbers[i])) {
        n /= iSmallPrimeNumbers[i];
        primes[iSmallPrimeNumbers[i]]++;
      }
    }
  }

  Type iDivisor(10007);
  while (n > 1) {
    if (!(n % iDivisor)) {
      primes[iDivisor] = 1;
      n /= iDivisor;

      while (!(n % iDivisor)) {
        n /= iDivisor;
        primes[iDivisor]++;
      }
    }

    iDivisor += 2;
  }

  return primes;
}

/*! 	\fn iRemoveSquareFactors
 * 	\brief Remove square factors of a positive integer
 *
 * 	\param iN(unsigned int) Integer
 * 	\return iN divided by all its squared factors
 */
template <typename Type> Type iRemoveSquareFactors(Type n) {
  Type nWork(n > 0 ? n : -n), res(1), square;

  if (nWork < 3)
    return n;

  for (unsigned int i(0); i < 1229 && nWork > 1; i++) {
    square = iSmallPrimeNumbers[i] * iSmallPrimeNumbers[i];

    while (nWork % square == 0)
      nWork /= square;

    if (nWork % iSmallPrimeNumbers[i] == 0) {
      nWork /= iSmallPrimeNumbers[i];
      res *= iSmallPrimeNumbers[i];
    }
  }

  Type iDivisor(10007);
  while (nWork > 1) {
    square = iDivisor * iDivisor;

    while (nWork % square == 0)
      nWork /= square;

    if (nWork % iDivisor == 0) {
      nWork /= iDivisor;
      res *= iDivisor;
    }

    iDivisor += 2;
  }

  return (n < 0 ? -res : res);
}

/*! 	\fn ceilQuotient
 * 	\brief Compute the ceiling of a quotient
 *
 * 	\param numerator(const unsigned int&) Numerator
 * 	\param denominator(const unsigned int&) Denominator
 * 	\return ceil(numerator / denominator)
 */
template <typename Type>
inline Type ceilQuotient(const Type &numerator, const Type &denominator) {
  return ((numerator % denominator) ? numerator / denominator + 1
                                    : numerator / denominator);
}

/*! 	\fn sqrtQuotient
 * 	\brief Compute the greatest integer less than or equal to the square
 * root of the given rational number
 *
 * 	\param numerator(const unsigned &) Numerator
 * 	\param denominator(const unsigned &) Denominator
 * 	\return ceil(sqrt(numerator / denominator))
 */
template <typename Type>
typename std::enable_if<std::is_unsigned<Type>::value, Type>::type
sqrtQuotient(const Type &numerator, const Type &denominator) {
  Type tRes(integerSqrt(numerator / denominator));

  while (denominator * (tRes + 1) * (tRes + 1) <= numerator)
    tRes++;

  return tRes;
}

mpz_class sqrtQuotient(const mpz_class &numerator,
                       const mpz_class &denominator);
mpz_class sqrtSupQuotient(const mpz_class &numerator,
                          const mpz_class &denominator);

/*! 	\fn sqrtSupQuotient
 * 	\brief Compute the smalles integer greater than or equal to the square
 * root of the given rational number
 *
 * 	\param numerator(const unsigned &) Numerator
 * 	\param denominator(const unsigned &) Denominator
 * 	\return ceil(sqrt(numerator / denominator))
 */
template <typename Type>
typename std::enable_if<std::is_unsigned<Type>::value, Type>::type
sqrtSupQuotient(const Type &numerator, const Type &denominator) {
  if (!numerator)
    return 0;

  Type tRes((numerator % denominator) != 0 ? numerator / denominator + 1
                                           : numerator / denominator);
  tRes = sqrtSup(tRes);

  while (numerator <= denominator * (tRes - 1) * (tRes - 1))
    tRes--;

  return tRes;
}
} // namespace MathTools
#endif
