/*
Copyright (C) 2013, 2014, 2015, 2016, 2017
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

#ifndef __TESTS_H__
#define __TESTS_H__

#include "../arithmeticity.h"
#include "../coxiter.h"
#include "../growthrate.h"
#include "../lib/numbers/mpz_rational.h"
#include "../lib/regexp.h"
#include "../signature.h"

#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <omp.h>
#include <string>
#include <vector>

using namespace std;

/*!
 * \file tests.h
 * \author Rafael Guglielmetti
 */

/*! \struct Test
 * \brief Contains information to perform a test
 */
struct Test {
  string szFile; ///< Path to the .coxiter file

  bool bTestEuler;       ///< If we have to check the Euler characteristic
  MPZ_rational brResult; ///< Theoretical Euler characteristic

  bool bTestFVector;             ///< If we have to check the f-vector
  vector<unsigned int> iFVector; ///< Theoretical value

  bool bTestCompacity;     ///< True if we test the cocompacity
  bool bIsCompact;         ///< True if the group is cocompact, false otherwise
  bool bIsFiniteVolume;    ///< True if finite covolume
  bool bTestArithmeticity; ///< True if we test the arithmeticity
  bool bIsArithmetic;      ///< True if arithmetic

  bool bTestGrowthSeries;
  vector<mpz_class>
      growthSeries_iPolynomialDenominator; ///< (i-1)th term contains the
                                           ///< coefficient of x^i
  vector<unsigned int>
      growthSeries_iCyclotomicNumerator; ///< Contains a list oif cyclotomic
                                         ///< polynomials
  string strGrowthRate;                  ///< Empty or the growth rate
};

/*! \class Tests
 * \brief Class to perform the tests
 */
class Tests {
private:
  string strError;         ///< Eventually, error code
  string strInputFilename; ///< File which contains the tests
  vector<Test> tests;      ///< Tests toi be performed

  ofstream of; ///< Output to file

  map<string, array<unsigned int, 2>>
      iTestsSucceded; ///< For each test, number of success, failure
  map<string, string> strTestDescription; ///< For each test, description
  unsigned int iTestsUnknownErrors;       ///< Number of other errors

public:
  Tests();

  /*!
   * 	\fn readGraphsFile
   * 	\brief Read the graph which contains the tests
   *
   * 	\param szInputFilename(string) File to read
   * 	\return True if success
   */
  bool readGraphsFile(string szInputFilename);

  /*!
   * 	\fn bRunTests
   * 	\brief Perform the tests
   */
  bool bRunTests();

  string get_strError() const;
  vector<Test> get_tests() const;
  unsigned int get_iTestsCount() const;

private:
  bool bRunTests_computations(const unsigned int &iTestIndex, CoxIter *ci);
  void bRunTests_arithmeticity(const unsigned int &iTestIndex, CoxIter *ci);
  void bRunTests_cocompactness_cofiniteness(const unsigned int &iTestIndex,
                                            CoxIter *ci);
  void bRunTests_growth(const unsigned int &iTestIndex, CoxIter *ci);
  void bRunTests_signature(const unsigned int &iTestIndex, CoxIter *ci,
                           const unsigned int &iDim);
  void bRunTests_Euler(const unsigned int &iTestIndex, CoxIter *ci);
  void bRunTests_FVector(const unsigned int &iTestIndex, CoxIter *ci);

  void bRunTests_error(const unsigned int &iTestIndex, const string &strTest,
                       const string &strExpected, const string &strComputed);

  void bRunTests_init();
  void bRunTests_displayInfo();

  string strIntToString(const int &i); ///< [*, 0, 1] => "?", "no", "yes"
};

#endif
