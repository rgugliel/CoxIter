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
 * \file polynomials.h
 * \author Rafael Guglielmetti
 *
 * \brief Some mathematical functions regarding polynomials
 */

#ifndef __POLYNOMIALS_H__
#define __POLYNOMIALS_H__

#include <cmath>
#include <iostream>
#include <iterator>
#include <map>
#include <string>
#include <vector>

#ifdef _USE_LOCAL_GMP_
#include "gmpxx.h"
#else
#include <gmpxx.h>
#endif

using namespace std;

namespace Polynomials {
/*! 	\fn polynomialDisplay
 * 	\brief Display a polynomial
 * 	\param polynomial(const vector< int >& polynomial) Integer
 */
template <typename Type>
void polynomialDisplay(const vector<Type> &polynomial) {
  bool isFirst(true);
  unsigned int iSize(polynomial.size());

  for (unsigned int i(0); i < iSize; i++) {
    if (polynomial[i] != 0) {
      if (isFirst) {
        cout << polynomial[i]
             << (i ? " * x" + string(i > 1 ? "^" + to_string(i) : "") : "");
        isFirst = false;
      } else {
        if ((polynomial[i] != 1 && polynomial[i] != -1) || !i)
          cout << (polynomial[i] > 0 ? " + " : " - ") << abs(polynomial[i])
               << (i ? " * x" + string(i > 1 ? "^" + to_string(i) : "") : "");
        else
          cout << (polynomial[i] > 0 ? " + " : " - ")
               << "x" + string(i > 1 ? "^" + to_string(i) : "");
      }
    }
  }
}

/*! 	\fn symbolDisplay
 * 	\brief Display a symbol
 * 	\param symbol(const vector< int >& symbol) Integer
 */
template <typename Type> void symbolDisplay(const vector<Type> &symbol) {
  bool isFirst(true);
  unsigned int iSize(symbol.size());

  cout << "[";
  for (unsigned int i(0); i < iSize; i++) {
    if (symbol[i]) {
      for (unsigned int j(0); j < symbol[i]; j++) {
        cout << (isFirst ? "" : ",") << i;
        isFirst = false;
      }
    }
  }
  cout << "]";
}

/*! 	\fn polynomialDotSymbol
 * 	\brief Multiply a polynomial by a symbol
 * 	\param polynomial(const vector< Type >& polynomial) The polynomial
 * 	\param symbol(const unsigned int&) The symbol
 */
template <typename Type>
void polynomialDotSymbol(vector<Type> &polynomial, const unsigned int &symbol) {
  vector<Type> polynomialBackup(polynomial);
  unsigned int polynomialDegree(polynomial.size() - 1);

  for (unsigned int i(1); i < symbol - 1; i++)
    polynomial.push_back(0);

  polynomial.push_back(polynomial[polynomialDegree]);

  if (polynomialDegree < symbol) {
    for (unsigned int i(1); i <= polynomialDegree; i++)
      polynomial[i] = polynomial[i - 1] + polynomialBackup[i];

    fill(polynomial.begin() + polynomialDegree + 1,
         polynomial.end() - polynomialDegree, polynomial[polynomialDegree]);

    for (unsigned int i(1); i <= polynomialDegree;
         i++) // TODO: gérer degré plus petit que symbole
      polynomial[polynomialDegree + symbol - i - 1] =
          polynomial[polynomialDegree + symbol - i] +
          polynomialBackup[polynomialDegree - i];
  } else {
    for (unsigned int i(1); i < symbol; i++)
      polynomial[i] = polynomial[i - 1] + polynomialBackup[i];

    for (unsigned int i(symbol); i < polynomialDegree; i++)
      polynomial[i] = polynomial[i - 1] - polynomialBackup[i - symbol] +
                      polynomialBackup[i];

    for (unsigned int i(1); i < symbol && i <= polynomialDegree; i++)
      polynomial[polynomialDegree + symbol - i - 1] =
          polynomial[polynomialDegree + symbol - i] +
          polynomialBackup[polynomialDegree - i];
  }
}

template <typename Type>
bool dividePolynomialBySymbol(vector<Type> &polynomial,
                              const unsigned int &symbol) {
  unsigned int polynomialDegree(polynomial.size() - 1);

  // Removing eventual 0
  while (polynomial[polynomialDegree] == 0)
    polynomialDegree--;

  vector<Type> iWorking(polynomial.begin(),
                        polynomial.begin() + polynomialDegree + 1);
  vector<Type> quotient;

  unsigned int i;
  Type temp;

  if (polynomialDegree < symbol - 1)
    return false;

  while (polynomialDegree >= symbol) {
    temp = iWorking[polynomialDegree];
    quotient.insert(quotient.begin(), temp);

    for (i = 0; i < symbol; i++)
      iWorking[polynomialDegree - i] -= temp;

    polynomialDegree--;

    while (iWorking[polynomialDegree] == 0 && polynomialDegree >= 1) {
      quotient.insert(quotient.begin(), 0);
      polynomialDegree--;
    }
  }

  if (polynomialDegree < symbol - 1) {
    for (i = 0; i <= polynomialDegree; i++) {
      if (iWorking[i] != 0)
        return false;
    }
  }

  temp = iWorking[polynomialDegree];
  quotient.insert(quotient.begin(), temp);

  for (i = 0; i < polynomialDegree; i++) {
    if (iWorking[i] != temp)
      return false;
  }

  polynomial = quotient;

  return true;
}

/*!	\fn dividePolynomialByPolynomial
 * 	\brief Try to make a division
 *
 * 	\param numerator(vector< Type >&) The dividend ; updated if the
 * remainder is 0 \param denominator(const vector< Type >) The divisor \return
 * bool: True if numerator is divisible by denominator. In this case,
 * numerator is updated to the quotient
 */
template <typename Type>
bool dividePolynomialByPolynomial(vector<Type> &numerator,
                                  const vector<Type> &denominator) {
  unsigned int numDeg(numerator.size() - 1), denomDeg(denominator.size() - 1);

  if (numDeg < denomDeg || (denomDeg == 0 && denominator[0] != 0))
    return false;

  vector<Type> working(numerator), quotient;

  unsigned int i;
  Type temp;

  while (numDeg >= denomDeg) {
    if (working[numDeg] % denominator[denomDeg] != 0)
      return false;

    temp = working[numDeg] / denominator[denomDeg];
    quotient.insert(quotient.begin(), temp);

    for (i = 0; i <= denomDeg; i++)
      working[numDeg - i] -= temp * denominator[denomDeg - i];

    numDeg--;

    while (working[numDeg] == 0 && numDeg >= 1 && numDeg >= denomDeg) {
      quotient.insert(quotient.begin(), 0);
      numDeg--;
    }
  }

  for (i = 0; i <= numDeg; i++) {
    if (working[i] != 0)
      return false;
  }

  numerator = quotient;

  return true;
}

extern vector<vector<mpz_class>>
    cyclotomicPolynomials; ///< List of some cyclotomic polynomials (we want to
                           ///< be able to multiply/divide with the growth
                           ///< series so we use here BigInteger instead of
                           ///< int)
} // namespace Polynomials

#endif
