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
 * 	\param iPolynomial(const vector< int >& iPolynomial) Integer
 */
template <typename Type>
void polynomialDisplay(const vector<Type> &iPolynomial) {
  bool bFirst(true);
  unsigned int iSize(iPolynomial.size());

  for (unsigned int i(0); i < iSize; i++) {
    if (iPolynomial[i] != 0) {
      if (bFirst) {
        cout << iPolynomial[i]
             << (i ? " * x" + string(i > 1 ? "^" + to_string(i) : "") : "");
        bFirst = false;
      } else {
        if ((iPolynomial[i] != 1 && iPolynomial[i] != -1) || !i)
          cout << (iPolynomial[i] > 0 ? " + " : " - ") << abs(iPolynomial[i])
               << (i ? " * x" + string(i > 1 ? "^" + to_string(i) : "") : "");
        else
          cout << (iPolynomial[i] > 0 ? " + " : " - ")
               << "x" + string(i > 1 ? "^" + to_string(i) : "");
      }
    }
  }
}

/*! 	\fn symbolDisplay
 * 	\brief Display a symbol
 * 	\param iPolynomial(const vector< int >& iPolynomial) Integer
 */
template <typename Type> void symbolDisplay(const vector<Type> &iSymbol) {
  bool bFirst(true);
  unsigned int iSize(iSymbol.size());

  cout << "[";
  for (unsigned int i(0); i < iSize; i++) {
    if (iSymbol[i]) {
      for (unsigned int j(0); j < iSymbol[i]; j++) {
        cout << (bFirst ? "" : ",") << i;
        bFirst = false;
      }
    }
  }
  cout << "]";
}

/*! 	\fn polynomialDotSymbol
 * 	\brief Multiply a polynomial by a symbol
 * 	\param iPolynomial(const vector< Type >& iPolynomial) The polynomial
 * 	\param iSymbol(const unsigned int&) The symbol
 */
template <typename Type>
void polynomialDotSymbol(vector<Type> &iPolynomial,
                         const unsigned int &iSymbol) {
  vector<Type> iPolynomialBackup(iPolynomial);
  unsigned int iPolynomialDegree(iPolynomial.size() - 1);

  for (unsigned int i(1); i < iSymbol - 1; i++)
    iPolynomial.push_back(0);
  iPolynomial.push_back(iPolynomial[iPolynomialDegree]);

  if (iPolynomialDegree < iSymbol) {
    for (unsigned int i(1); i <= iPolynomialDegree; i++)
      iPolynomial[i] = iPolynomial[i - 1] + iPolynomialBackup[i];

    fill(iPolynomial.begin() + iPolynomialDegree + 1,
         iPolynomial.end() - iPolynomialDegree, iPolynomial[iPolynomialDegree]);

    for (unsigned int i(1); i <= iPolynomialDegree;
         i++) // TODO: gérer degré plus petit que symbole
      iPolynomial[iPolynomialDegree + iSymbol - i - 1] =
          iPolynomial[iPolynomialDegree + iSymbol - i] +
          iPolynomialBackup[iPolynomialDegree - i];
  } else {
    for (unsigned int i(1); i < iSymbol; i++)
      iPolynomial[i] = iPolynomial[i - 1] + iPolynomialBackup[i];

    for (unsigned int i(iSymbol); i < iPolynomialDegree; i++)
      iPolynomial[i] = iPolynomial[i - 1] - iPolynomialBackup[i - iSymbol] +
                       iPolynomialBackup[i];

    for (unsigned int i(1); i < iSymbol && i <= iPolynomialDegree; i++)
      iPolynomial[iPolynomialDegree + iSymbol - i - 1] =
          iPolynomial[iPolynomialDegree + iSymbol - i] +
          iPolynomialBackup[iPolynomialDegree - i];
  }
}

template <typename Type>
bool dividePolynomialBySymbol(vector<Type> &iPolynomial,
                              const unsigned int &iSymbol) {
  unsigned int iPolynomialDegree(iPolynomial.size() - 1);

  // Removing eventual 0
  while (iPolynomial[iPolynomialDegree] == 0)
    iPolynomialDegree--;

  vector<Type> iWorking(iPolynomial.begin(),
                        iPolynomial.begin() + iPolynomialDegree + 1);
  vector<Type> iQuotient;

  unsigned int i;
  Type iTemp;

  if (iPolynomialDegree < iSymbol - 1)
    return false;

  while (iPolynomialDegree >= iSymbol) {
    iTemp = iWorking[iPolynomialDegree];
    iQuotient.insert(iQuotient.begin(), iTemp);

    for (i = 0; i < iSymbol; i++)
      iWorking[iPolynomialDegree - i] -= iTemp;

    iPolynomialDegree--;

    while (iWorking[iPolynomialDegree] == 0 && iPolynomialDegree >= 1) {
      iQuotient.insert(iQuotient.begin(), 0);
      iPolynomialDegree--;
    }
  }

  if (iPolynomialDegree < iSymbol - 1) {
    for (i = 0; i <= iPolynomialDegree; i++) {
      if (iWorking[i] != 0)
        return false;
    }
  }

  iTemp = iWorking[iPolynomialDegree];
  iQuotient.insert(iQuotient.begin(), iTemp);

  for (i = 0; i < iPolynomialDegree; i++) {
    if (iWorking[i] != iTemp)
      return false;
  }

  iPolynomial = iQuotient;

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
  Type iTemp;

  while (numDeg >= denomDeg) {
    if (working[numDeg] % denominator[denomDeg] != 0)
      return false;

    iTemp = working[numDeg] / denominator[denomDeg];
    quotient.insert(quotient.begin(), iTemp);

    for (i = 0; i <= denomDeg; i++)
      working[numDeg - i] -= iTemp * denominator[denomDeg - i];

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
    iCyclotomicPolynomials; ///< List of some cyclotomic polynomials (we want to
                            ///< be able to multiply/divide with the growth
                            ///< series so we use here BigInteger instead of
                            ///< int)
} // namespace Polynomials

#endif
