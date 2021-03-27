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
  string filename; ///< Path to the .coxiter file

  bool testEuler;        ///< If we have to check the Euler characteristic
  MPZ_rational brResult; ///< Theoretical Euler characteristic

  bool testFVector;             ///< If we have to check the f-vector
  vector<unsigned int> fVector; ///< Theoretical value

  bool testCompacity;     ///< True if we test the cocompacity
  bool isCompact;         ///< True if the group is cocompact, false otherwise
  bool isFiniteVolume;    ///< True if finite covolume
  bool testArithmeticity; ///< True if we test the arithmeticity
  bool isArithmetic;      ///< True if arithmetic

  bool testGrowthSeries;
  vector<mpz_class>
      growthSeries_polynomialDenominator; ///< (i-1)th term contains the
                                          ///< coefficient of x^i
  vector<unsigned int>
      growthSeries_cyclotomicNumerator; ///< Contains a list oif cyclotomic
                                        ///< polynomials
  string growthRate;                    ///< Empty or the growth rate
};

/*! \class Tests
 * \brief Class to perform the tests
 */
class Tests {
private:
  string error;         ///< Eventually, error code
  string inputFilename; ///< File which contains the tests
  vector<Test> tests;   ///< Tests toi be performed

  ofstream of; ///< Output to file

  map<string, array<unsigned int, 2>>
      testsSucceded; ///< For each test, number of success, failure
  map<string, string> testDescription; ///< For each test, description
  unsigned int testsUnknownErrors;     ///< Number of other errors

public:
  Tests();

  /*!
   * 	\fn readGraphsFile
   * 	\brief Read the graph which contains the tests
   *
   * 	\param inputFilename(string) File to read
   * 	\return True if success
   */
  bool readGraphsFile(string inputFilename);

  /*!
   * 	\fn runTests
   * 	\brief Perform the tests
   */
  bool runTests();

  string get_error() const;
  vector<Test> get_tests() const;
  unsigned int get_testsCount() const;

private:
  bool runTests_computations(const unsigned int &testIndex, CoxIter *ci);
  void runTests_arithmeticity(const unsigned int &testIndex, CoxIter *ci);
  void runTests_cocompactness_cofiniteness(const unsigned int &testIndex,
                                           CoxIter *ci);
  void runTests_growth(const unsigned int &testIndex, CoxIter *ci);
  void runTests_signature(const unsigned int &testIndex, CoxIter *ci,
                          const unsigned int &dim);
  void runTests_euler(const unsigned int &testIndex, CoxIter *ci);
  void runTests_fVector(const unsigned int &testIndex, CoxIter *ci);

  void runTestsError(const unsigned int &testIndex, const string &test,
                     const string &expected, const string &computed);

  void runTests_init();
  void runTests_displayInfo();

  string strIntToString(const int &i); ///< [*, 0, 1] => "?", "no", "yes"
};

#endif
