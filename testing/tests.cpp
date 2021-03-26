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

#include "tests.h"

Tests::Tests() {}

bool Tests::readGraphsFile(string input) {
  // ------------------------------------------------------
  // ouverture fichier
  inputFilename = input;
  ifstream fileIn(inputFilename.c_str());
  if (fileIn.fail())
    throw(string("READ_INPUT_FILE"));

  // ------------------------------------------------------
  // lecture
  PCRERegexp regexp;
  PCREResult regexpRes;
  string line, lineWork, number, lineWithoutSpaces;

  size_t findPos;

  Test test;

  while (getline(fileIn, line)) // tant que....
  {
    lineWithoutSpaces = line;
    str_replace(lineWithoutSpaces, " ", "");
    if (lineWithoutSpaces == "")
      continue;

    test.testEuler = false;
    test.testFVector = false;
    test.testCompacity = false;
    test.isFiniteVolume = true;
    test.testArithmeticity = false;
    test.testGrowthSeries = false;
    test.growthRate = "";

    if ((findPos = line.find("#")) != string::npos) // Forget what is after #
      line = line.substr(0, findPos);

    // --------------------------------------------------------
    // nom du fichier à tester
    regexpRes.clear();
    if (regexp.preg_match_all("\"([^\"]+)\"", line, regexpRes) ==
        0) // If we cannot read the file
    {
      cout << "Line not read: " << line << endl;
      continue;
    }

    test.filename = regexpRes[1][0];
    lineWork = line.substr(test.filename.length() + 2); // We drop the name
    str_replace(lineWork, " ", "");

    // --------------------------------------------------------
    // Growth rate
    regexpRes.clear();
    if (regexp.preg_match_all("tau\\=([[:digit:].]+)[;]?", lineWork,
                              regexpRes) != 0) {
      str_replace(lineWork, regexpRes[0][0], " ");
      test.growthRate = regexpRes[1][0];
    }

    // --------------------------------------------------------
    // Growth series
    regexpRes.clear();
    if (regexp.preg_match_all("f\\(x\\)=[C\\(]*([0-9,]*)[\\)]?/"
                              "\\(([0-9\\-\\+\\*\\^x]+)\\)[;|\n]{1,1}",
                              lineWork, regexpRes)) {
      test.testGrowthSeries = true;
      str_replace(lineWork, regexpRes[0][0], "");

      unsigned int iPower;
      test.growthSeries_polynomialDenominator = vector<mpz_class>(1, 0);
      test.growthSeries_cyclotomicNumerator.clear();

      if (regexpRes[1][0] != "") // Cyclotomic factors
        explode(",", regexpRes[1][0], test.growthSeries_cyclotomicNumerator);

      string polynomial(regexpRes[2][0]);

      regexpRes.clear();
      int iResCount(regexp.preg_match_all(
          "(\\+|\\-?)([0-9]*)[\\*]?x[\\^]?([0-9]*)", polynomial, regexpRes));

      for (int i(iResCount - 1); i >= 0; i--) {
        iPower =
            regexpRes[3][i] == "" ? 1 : abs(stoi(regexpRes[3][i])); // power

        if (iPower + 1 > test.growthSeries_polynomialDenominator.size())
          test.growthSeries_polynomialDenominator.insert(
              test.growthSeries_polynomialDenominator.end(),
              iPower - test.growthSeries_polynomialDenominator.size() + 1,
              mpz_class(0));

        test.growthSeries_polynomialDenominator[iPower] =
            (regexpRes[1][i] == "-" ? -1 : 1) *
            (regexpRes[2][i] == "" ? 1 : stoi(regexpRes[2][i]));

        string strRep(
            regexpRes[1][i] + regexpRes[2][i] +
            (regexpRes[2][i] != "" && stoi(regexpRes[2][i]) > 1 ? "*" : "") +
            "x" + (iPower > 1 ? "^" : "") + regexpRes[3][i]);

        str_replace(polynomial, strRep, "");
      }

      if (polynomial != "")
        test.growthSeries_polynomialDenominator[0] = stoi(polynomial);
    }

    // --------------------------------------------------------
    // f-vecteur
    regexpRes.clear();
    if (regexp.preg_match_all("\\(([[:digit:] ,]+)\\)", lineWork, regexpRes) !=
        0) {
      str_replace(lineWork, "(" + regexpRes[1][0] + ")", " ");
      vector<string> fVector(explode(",", regexpRes[1][0]));
      test.fVector = vector<unsigned int>(0);
      test.testFVector = true;

      for (const auto &element : fVector)
        test.fVector.push_back(stoi(element));
    }

    // --------------------------------------------------------
    // caractéristique d'Euler
    regexpRes.clear();
    if (regexp.preg_match_all("([-]{0,1})([[:digit:]/]+)", lineWork,
                              regexpRes) != 0) {
      number = regexpRes[1][0] + regexpRes[2][0];
      test.testEuler = true;

      if (regexp.preg_match_all("([[:digit:]]+)", number, regexpRes) !=
          0) // If there is at least one digit
      {
        try {
          test.brResult = MPZ_rational(number);
        } catch (int iCode) {
          cout << "Error\t " << line << endl;
          cout << "\t\tError while reading Euler characteristic" << endl;
          cout << "\t\t" << number << " is not a valid number" << endl;
          test.testEuler = false;
        }
      } else
        test.testEuler = false;
    }

    // --------------------------------------------------------
    // Cocompactness
    if (lineWork.find("compact") != string::npos ||
        lineWork.find("cocompact") != string::npos) {
      test.testCompacity = true;
      test.isCompact = true;
    }
    if (lineWork.find("non-compact") != string::npos ||
        lineWork.find("non-cocompact") != string::npos) {
      test.testCompacity = true;
      test.isCompact = false;
    }

    // --------------------------------------------------------
    // Arithmeticity
    if (lineWork.find("arithmetic") != string::npos) {
      test.testCompacity = true;
      test.testArithmeticity = true;
      test.isArithmetic = true;
    }
    if (lineWork.find("non-arithmetic") != string::npos) {
      test.testCompacity = true;
      test.testArithmeticity = true;
      test.isArithmetic = false;
    }

    // --------------------------------------------------------
    // Cofiniteness
    if (lineWork.find("non-fv") != string::npos) {
      test.isFiniteVolume = false;
    }

    tests.push_back(test);
  }

  fileIn.close();

  return true;
}

bool Tests::runTests() {
  string str_tempFilename(inputFilename + ".output");
  of = ofstream(str_tempFilename.c_str());

  if (!of.is_open()) {
    error = "Cannot open: " + str_tempFilename;
    return false;
  }

  runTests_init();

  // ------------------------------------------------------
  // Let's compute
  CoxIter ci;

  int i, testsCount(tests.size());
  unsigned int iDim, iBarWidth(70);
  double progressStep(.02), progressNext(progressStep);

  for (unsigned int j(0); j < iBarWidth; ++j) {
    if (!j)
      std::cout << ">";
    else
      cout << " ";
  }
  cout << "] " << 0 << " %\r";
  cout.flush();

  for (i = 0; i < testsCount; i++) {
    // ---------------------------------------------
    // Progression bar
    if ((double)(i + 1) / testsCount > progressNext) {
      double dProgress((double)(i + 1) / testsCount);

      unsigned int iPos(iBarWidth * dProgress);
      for (unsigned int j(0); j < iBarWidth; ++j) {
        if (j < iPos)
          cout << "=";
        else if (j == iPos)
          cout << ">";
        else
          cout << " ";
      }
      cout << "] " << int(dProgress * 100.0) << " %\r";
      cout.flush();

      progressNext += progressStep;
    }

    // ---------------------------------------------
    // Computations
    ci = CoxIter();

    ci.set_checkCocompactness(tests[i].testCompacity);
    ci.set_checkCofiniteness(true);

    if (!ci.readGraphFromFile("../../graphs/" + tests[i].filename)) {
      of << "Error\t " << tests[i].filename << endl;
      of << "\t\tError when reading file(" << ci.get_error() << ")" << endl;
      testsSucceded["readingGraph"][1]++;
      continue;
    }

    testsSucceded["readingGraph"][0]++;

    iDim = ci.get_dimension();
    ci.set_dimension(0);

    // ---------------------------------------------
    // Computations
    if (!runTests_computations(i, &ci))
      continue;

    // ---------------------------------------------
    // Tests
    runTests_growth(i, &ci);

    runTests_signature(i, &ci, iDim);

    runTests_cocompactness_cofiniteness(i, &ci);

    if (tests[i].testArithmeticity)
      runTests_arithmeticity(i, &ci);

    runTests_euler(i, &ci);
    runTests_fVector(i, &ci);
  }

  cout << "\n" << endl;

  of.close();

  runTests_displayInfo();

  return true;
}

vector<Test> Tests::get_tests() const { return tests; }

unsigned int Tests::get_testsCount() const { return tests.size(); }

string Tests::get_error() const { return error; }

void Tests::runTests_init() {
  testsSucceded.clear();

  testsSucceded["arithmeticity"] = array<unsigned int, 2>{0, 0};
  testsSucceded["cocompactness"] = array<unsigned int, 2>{0, 0};
  testsSucceded["cofiniteness"] = array<unsigned int, 2>{0, 0};
  testsSucceded["cofinitenessPartial"] = array<unsigned int, 2>{0, 0};
  testsSucceded["dimensionGuess"] = array<unsigned int, 2>{0, 0};
  testsSucceded["euler"] = array<unsigned int, 2>{0, 0};
  testsSucceded["fv"] = array<unsigned int, 2>{0, 0};
  testsSucceded["fvAlt"] = array<unsigned int, 2>{0, 0};
  testsSucceded["growthRate"] = array<unsigned int, 2>{0, 0};
  testsSucceded["growthSeries"] = array<unsigned int, 2>{0, 0};
  testsSucceded["growthSeriesDenomDimOdd"] = array<unsigned int, 2>{0, 0};
  testsSucceded["growthSeriesEuler"] = array<unsigned int, 2>{0, 0};
  testsSucceded["readingGraph"] = array<unsigned int, 2>{0, 0};
  testsSucceded["signature"] = array<unsigned int, 2>{0, 0};

  testsUnknownErrors = 0;

  testDescription["arithmeticity"] = "Arithmeticity";
  testDescription["cocompactness"] = "Cocompactness";
  testDescription["cofiniteness"] = "Cofiniteness";
  testDescription["cofinitenessPartial"] = "Cofiniteness (partial)";
  testDescription["dimensionGuess"] = "Dimension guess";
  testDescription["euler"] = "Euler charactereistic";
  testDescription["fv"] = "f-vector";
  testDescription["fvAlt"] = "Alt. sum comp. f-vector";
  testDescription["growthRate"] = "Growth rate";
  testDescription["growthSeries"] = "Growth series";
  testDescription["growthSeriesDenomDimOdd"] =
      "Denom. growth series vanish at 1";
  testDescription["growthSeriesEuler"] = "Growth series <-> Euler char.";
  testDescription["readingGraph"] = "Reading graph";
  testDescription["signature"] = "Signature";
}

bool Tests::runTests_computations(const unsigned int &iTestIndex, CoxIter *ci) {
  ci->exploreGraph();
  ci->computeGraphsProducts();
  if (!ci->computeEulerCharacteristicFVector()) {
    of << "Error\t " << tests[iTestIndex].filename << endl;
    of << "\t\tPlease check the graph; maybe it is not hyperbolic" << endl;
    testsUnknownErrors++;
    return false;
  }

  if (tests[iTestIndex].testCompacity)
    ci->isGraphCocompact();

  ci->checkCovolumeFiniteness();
  ci->canBeFiniteCovolume();
  ci->growthSeries();

  return true;
}

void Tests::runTests_growth(const unsigned int &testIndex, CoxIter *ci) {
  string growthRate;
  GrowthRate *gr(new GrowthRate());
  GrowthRate_Result grr(
      gr->grrComputations(ci->get_growthSeries_denominator()));
  delete gr;
  growthRate = grr.growthRate;

  if (tests[testIndex].growthRate != "") {
    if (tests[testIndex].growthRate == growthRate) {
      testsSucceded["growthRate"][0]++;
      of << "OK\tGrowth rate\t\t" << tests[testIndex].filename << endl;
    } else {
      testsSucceded["growthRate"][1]++;
      runTestsError(testIndex, "growth rate", tests[testIndex].growthRate,
                    growthRate);
    }
  }

  if (grr.perron != 1)
    cout << "INFO: Not a Perron number in " << tests[testIndex].filename
         << endl;

  // A small test of the growth series
  if (ci->get_isFiniteCovolume() > 0) {
    vector<mpz_class> denom;
    vector<unsigned int> cyclotomic;
    bool isReduced;

    ci->get_growthSeries(cyclotomic, denom, isReduced);

    if (tests[testIndex].testGrowthSeries) {
      if (tests[testIndex].growthSeries_cyclotomicNumerator == cyclotomic &&
          tests[testIndex].growthSeries_polynomialDenominator == denom) {
        testsSucceded["growthSeries"][0]++;
        of << "OK\tGrowth series\t\t" << tests[testIndex].filename << endl;
      } else {
        testsSucceded["growthSeries"][1]++;
        of << "Error\t " << tests[testIndex].filename << endl;
        of << "\t\tGrowth series" << endl;
      }
    }

    mpz_class iTotalDenom(0);
    for (unsigned int j(0); j < denom.size(); j++)
      iTotalDenom += denom[j];

    if (ci->get_dimension() % 2) // n is odd, the denominator should vanish in 1
    {
      if (iTotalDenom == 0) {
        testsSucceded["growthSeriesDenomDimOdd"][0]++;
        of << "OK\tGrowth series test\t" << tests[testIndex].filename << endl;
      } else {
        testsSucceded["growthSeriesDenomDimOdd"][1]++;
        of << "Error\t " << tests[testIndex].filename << endl;
        of << "\t\tGrowth series test: the denominator should vanish at 1"
           << endl;
      }
    } else // n even, we should find the euler characteristics
    {
      mpz_class iTotalNum(1), iSum;

      for (auto cyclo : cyclotomic) {
        iSum = 0;
        for (auto coeff : Polynomials::cyclotomicPolynomials[cyclo])
          iSum += coeff;

        iTotalNum *= iSum;
      }

      if (MPZ_rational(iTotalDenom, iTotalNum) ==
          ci->get_brEulerCaracteristic()) {
        testsSucceded["growthSeriesEuler"][0]++;
        of << "OK\tGrowth series test\t" << tests[testIndex].filename << endl;
      } else {
        testsSucceded["growthSeriesEuler"][1]++;
        of << "Error\t " << tests[testIndex].filename << endl;
        of << "\t\tGrowth series test: different from Euler characteristic"
           << endl;
      }
    }
  }
}

void Tests::runTests_signature(const unsigned int &testIndex, CoxIter *ci,
                               const unsigned int &dim) {
  int signatureComputed(0);

  array<unsigned int, 3> iSignature;

  try {
    if (ci->get_hasDottedLineWithoutWeight() == 0) {
      Signature s;
      iSignature = s.computeSignature(ci->get_gramMatrix_PARI());
      signatureComputed = 1;
    } else
      signatureComputed = 0;
  } catch (string ex) {
    signatureComputed = -1;
  }

  if (signatureComputed == 1 && ci->get_dimension()) {
    if (iSignature[0] != ci->get_dimension() || iSignature[1] != 1) {
      testsSucceded["signature"][1]++;
      of << "Error\t " << tests[testIndex].filename << endl;
      of << "\t\tSignature computed: (" << iSignature[0] << "," << iSignature[1]
         << "," << iSignature[2] << ")" << endl;
    } else {
      testsSucceded["signature"][0]++;
      of << "OK\tSignature\t\t" << tests[testIndex].filename << endl;
    }
  } else if (signatureComputed == -1) {
    testsSucceded["signature"][1]++;
    of << "Error\t " << tests[testIndex].filename << endl;
    of << "\t\tSignature: format" << endl;
  } else if (signatureComputed == 0 &&
             tests[testIndex].filename.find("roberts", 0) != string::npos) {
    testsSucceded["signature"][1]++;
    of << "Error\t " << tests[testIndex].filename << endl;
    of << "\t\tWeights unkown" << endl;
  }

  // ---------------------------------------
  // Guessing the dimension
  if (ci->get_dimensionGuessed()) {
    if (dim == ci->get_dimension()) {
      testsSucceded["dimensionGuess"][0]++;
      of << "OK\tDimension guessed\t" << tests[testIndex].filename << endl;
    } else {
      testsSucceded["dimensionGuess"][1]++;
      of << "Error\tDimension guessed\t" << tests[testIndex].filename << endl;
    }
  }
}

void Tests::runTests_arithmeticity(const unsigned int &testIndex, CoxIter *ci) {
  int isArithmetic;

  Arithmeticity arithmeticity;
  arithmeticity.test(*ci, false);

  isArithmetic = ci->get_isArithmetic();
  if ((isArithmetic == 1 && tests[testIndex].isArithmetic) ||
      (isArithmetic == 0 && !tests[testIndex].isArithmetic)) {
    testsSucceded["arithmeticity"][0]++;
    of << "OK\tArithmeticity\t\t" << tests[testIndex].filename << endl;
  } else {
    testsSucceded["arithmeticity"][1]++;
    runTestsError(testIndex, "arithmeticity",
                  tests[testIndex].isArithmetic ? "yes" : "no",
                  strIntToString(isArithmetic));
  }
}

void Tests::runTests_cocompactness_cofiniteness(const unsigned int &testIndex,
                                                CoxIter *ci) {
  int compacity, finiteVolume(ci->checkCovolumeFiniteness());
  bool canBeFiniteCovolume(ci->canBeFiniteCovolume());

  if (finiteVolume == 1 && !canBeFiniteCovolume) {
    testsSucceded["cofinitenessPartial"][1]++;
    of << "Error\t " << tests[testIndex].filename << endl;
    of << "\t\tTest CanBeFiniteCovolume" << endl;
  } else if ((finiteVolume == 0 && !canBeFiniteCovolume) ||
             (finiteVolume == 1 && canBeFiniteCovolume)) {
    testsSucceded["cofinitenessPartial"][0]++;
    of << "OK\tCanBeFiniteCovolume\t" << tests[testIndex].filename << endl;
  }

  if (!((finiteVolume == 1 && tests[testIndex].isFiniteVolume) ||
        (finiteVolume == 0 && !tests[testIndex].isFiniteVolume))) {
    testsSucceded["cofiniteness"][1]++;
    runTestsError(testIndex, "cofiniteness test",
                  tests[testIndex].isFiniteVolume ? "yes" : "no",
                  strIntToString(finiteVolume));
  } else {
    testsSucceded["cofiniteness"][0]++;
    of << "OK\tCofiniteness test\t" << tests[testIndex].filename << endl;
  }

  if (tests[testIndex].testCompacity) {
    compacity = ci->get_isCocompact();
    if ((compacity == 1 && tests[testIndex].isCompact) ||
        (compacity == 0 && !tests[testIndex].isCompact)) {
      testsSucceded["cocompactness"][0]++;
      of << "OK\tCompacity\t\t" << tests[testIndex].filename << endl;
    } else {
      testsSucceded["cocompactness"][1]++;
      runTestsError(testIndex, "cocompactness test",
                    tests[testIndex].isCompact ? "yes" : "no",
                    strIntToString(compacity));
    }
  }
}

void Tests::runTests_euler(const unsigned int &iTestIndex, CoxIter *ci) {
  if (tests[iTestIndex].testEuler) {
    if (ci->get_brEulerCaracteristic() != tests[iTestIndex].brResult) {
      testsSucceded["euler"][1]++;
      runTestsError(iTestIndex, "Euler characteristic",
                    tests[iTestIndex].brResult.to_string(),
                    ci->get_brEulerCaracteristic().to_string());
    } else {
      testsSucceded["euler"][0]++;
      of << "OK\tEuler\t\t\t" << tests[iTestIndex].filename << endl;
    }
  } else if (ci->get_dimension() % 2 && tests[iTestIndex].isFiniteVolume) {
    if (ci->get_brEulerCaracteristic() != 0) {
      testsSucceded["euler"][1]++;
      runTestsError(iTestIndex, "Euler characteristic", "0",
                    ci->get_brEulerCaracteristic().to_string());
    } else {
      testsSucceded["euler"][0]++;
      of << "OK\tEuler\t\t\t" << tests[iTestIndex].filename << endl;
    }
  }
}

void Tests::runTests_fVector(const unsigned int &iTestIndex, CoxIter *ci) {
  if (ci->get_isFiniteCovolume() > 0) {
    if (ci->get_fVectorAlternateSum() == (ci->get_dimension() % 2 ? 2 : 0)) {
      testsSucceded["fvAlt"][0]++;
      of << "OK\tAlt. sum\t\t" << tests[iTestIndex].filename << endl;
    } else {
      testsSucceded["fvAlt"][1]++;
      runTestsError(iTestIndex, "alt. sum of components of f-vector",
                    to_string(ci->get_dimension() % 2 ? 2 : 0),
                    to_string(ci->get_fVectorAlternateSum()));
    }
  }

  if (tests[iTestIndex].testFVector) {
    if (ci->get_fVector() == tests[iTestIndex].fVector) {
      testsSucceded["fv"][0]++;
      of << "OK\tf-vector\t\t" << tests[iTestIndex].filename << endl;
    } else {
      testsSucceded["fv"][1]++;
      runTestsError(iTestIndex, "f-vector",
                    implode(",", tests[iTestIndex].fVector),
                    implode(",", ci->get_fVector()));
    }
  }
}

void Tests::runTestsError(const unsigned int &testIndex, const string &test,
                          const string &expected, const string &computed) {
  of << "Error: " << test << "\t in " << tests[testIndex].filename << endl;
  of << "\t\tExpected: " << expected << endl;
  of << "\t\tComputed: " << computed << endl;
}

string Tests::strIntToString(const int &i) {
  if (i == 0)
    return "no";
  else if (i == 1)
    return "yes";
  else
    return "?";
}

void Tests::runTests_displayInfo() {
  array<unsigned int, 2> iTotal({0, 0});

  cout << "------------------------------------------------------" << endl;
  cout << setw(32) << left << "Test"
       << " | #success | #error |" << endl;
  cout << "------------------------------------------------------" << endl;

  for (map<string, string>::const_iterator it(testDescription.begin());
       it != testDescription.end(); it++) {
    iTotal[0] += testsSucceded[it->first][0];
    iTotal[1] += testsSucceded[it->first][1];

    cout << setw(32) << left << it->second << " | ";
    cout << setw(8) << left << testsSucceded[it->first][0] << " | ";
    cout << setw(6) << left << testsSucceded[it->first][1] << " | ";
    cout << endl;
  }

  cout << "------------------------------------------------------" << endl;
  cout << setw(32) << left << "Total"
       << " | ";
  cout << setw(8) << left << iTotal[0] << " | ";
  cout << setw(6) << left << iTotal[1] << " | ";
  cout << endl;
  cout << "------------------------------------------------------\n" << endl;

  if (testsUnknownErrors)
    cout << "Other errors: " << testsUnknownErrors << endl;

  cout << "For more information, see:\n\t" << inputFilename << ".output"
       << "\n"
       << endl;
}
