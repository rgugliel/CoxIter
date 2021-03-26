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

bool Tests::readGraphsFile(string strInput) {
  // ------------------------------------------------------
  // ouverture fichier
  strInputFilename = strInput;
  ifstream fileIn(strInputFilename.c_str());
  if (fileIn.fail())
    throw(string("READ_INPUT_FILE"));

  // ------------------------------------------------------
  // lecture
  PCRERegexp regexp;
  PCREResult regexpRes;
  string szLine, szLineWork, szNumber, strEmpty;

  size_t iFindPos;

  Test test;

  while (getline(fileIn, szLine)) // tant que....
  {
    strEmpty = szLine;
    str_replace(strEmpty, " ", "");
    if (strEmpty == "")
      continue;

    test.bTestEuler = false;
    test.bTestFVector = false;
    test.bTestCompacity = false;
    test.bIsFiniteVolume = true;
    test.bTestArithmeticity = false;
    test.bTestGrowthSeries = false;
    test.strGrowthRate = "";

    if ((iFindPos = szLine.find("#")) != string::npos) // Forget what is after #
      szLine = szLine.substr(0, iFindPos);

    // --------------------------------------------------------
    // nom du fichier à tester
    regexpRes.clear();
    if (regexp.preg_match_all("\"([^\"]+)\"", szLine, regexpRes) ==
        0) // If we cannot read the file
    {
      cout << "Line not read: " << szLine << endl;
      continue;
    }

    test.szFile = regexpRes[1][0];
    szLineWork = szLine.substr(test.szFile.length() + 2); // We drop the name
    str_replace(szLineWork, " ", "");

    // --------------------------------------------------------
    // Growth rate
    regexpRes.clear();
    if (regexp.preg_match_all("tau\\=([[:digit:].]+)[;]?", szLineWork,
                              regexpRes) != 0) {
      str_replace(szLineWork, regexpRes[0][0], " ");
      test.strGrowthRate = regexpRes[1][0];
    }

    // --------------------------------------------------------
    // Growth series
    regexpRes.clear();
    if (regexp.preg_match_all("f\\(x\\)=[C\\(]*([0-9,]*)[\\)]?/"
                              "\\(([0-9\\-\\+\\*\\^x]+)\\)[;|\n]{1,1}",
                              szLineWork, regexpRes)) {
      test.bTestGrowthSeries = true;
      str_replace(szLineWork, regexpRes[0][0], "");

      unsigned int iPower;
      test.growthSeries_polynomialDenominator = vector<mpz_class>(1, 0);
      test.growthSeries_cyclotomicNumerator.clear();

      vector<unsigned int> iTemp;

      if (regexpRes[1][0] != "") // Cyclotomic factors
        explode(",", regexpRes[1][0], test.growthSeries_cyclotomicNumerator);

      string strPol(regexpRes[2][0]);

      regexpRes.clear();
      int iResCount(regexp.preg_match_all(
          "(\\+|\\-?)([0-9]*)[\\*]?x[\\^]?([0-9]*)", strPol, regexpRes));

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

        str_replace(strPol, strRep, "");
      }

      if (strPol != "")
        test.growthSeries_polynomialDenominator[0] = stoi(strPol);
    }

    // --------------------------------------------------------
    // f-vecteur
    regexpRes.clear();
    if (regexp.preg_match_all("\\(([[:digit:] ,]+)\\)", szLineWork,
                              regexpRes) != 0) {
      str_replace(szLineWork, "(" + regexpRes[1][0] + ")", " ");
      vector<string> strFVector(explode(",", regexpRes[1][0]));
      test.iFVector = vector<unsigned int>(0);
      test.bTestFVector = true;

      for (vector<string>::const_iterator strIt(strFVector.begin());
           strIt != strFVector.end(); ++strIt)
        test.iFVector.push_back(stoi(*strIt));
    }

    // --------------------------------------------------------
    // caractéristique d'Euler
    regexpRes.clear();
    if (regexp.preg_match_all("([-]{0,1})([[:digit:]/]+)", szLineWork,
                              regexpRes) != 0) {
      szNumber = regexpRes[1][0] + regexpRes[2][0];
      test.bTestEuler = true;

      if (regexp.preg_match_all("([[:digit:]]+)", szNumber, regexpRes) !=
          0) // If there is at least one digit
      {
        try {
          test.brResult = MPZ_rational(szNumber);
        } catch (int iCode) {
          cout << "Error\t " << szLine << endl;
          cout << "\t\tError while reading Euler characteristic" << endl;
          cout << "\t\t" << szNumber << " is not a valid number" << endl;
          test.bTestEuler = false;
        }
      } else
        test.bTestEuler = false;
    }

    // --------------------------------------------------------
    // Cocompactness
    if (szLineWork.find("compact") != string::npos ||
        szLineWork.find("cocompact") != string::npos) {
      test.bTestCompacity = true;
      test.bIsCompact = true;
    }
    if (szLineWork.find("non-compact") != string::npos ||
        szLineWork.find("non-cocompact") != string::npos) {
      test.bTestCompacity = true;
      test.bIsCompact = false;
    }

    // --------------------------------------------------------
    // Arithmeticity
    if (szLineWork.find("arithmetic") != string::npos) {
      test.bTestCompacity = true;
      test.bTestArithmeticity = true;
      test.bIsArithmetic = true;
    }
    if (szLineWork.find("non-arithmetic") != string::npos) {
      test.bTestCompacity = true;
      test.bTestArithmeticity = true;
      test.bIsArithmetic = false;
    }

    // --------------------------------------------------------
    // Cofiniteness
    if (szLineWork.find("non-fv") != string::npos) {
      test.bIsFiniteVolume = false;
    }

    tests.push_back(test);
  }

  fileIn.close();

  return true;
}

bool Tests::bRunTests() {
  string str_tempFilename(strInputFilename + ".output");
  of = ofstream(str_tempFilename.c_str());

  if (!of.is_open()) {
    strError = "Cannot open: " + str_tempFilename;
    return false;
  }

  bRunTests_init();

  // ------------------------------------------------------
  // Let's compute
  CoxIter ci;

  int i, iMax(tests.size());
  unsigned int iDim, iBarWidth(70);
  double dProgressStep(.02), dProgressNext(dProgressStep);

  for (unsigned int j(0); j < iBarWidth; ++j) {
    if (!j)
      std::cout << ">";
    else
      cout << " ";
  }
  cout << "] " << 0 << " %\r";
  cout.flush();

  for (i = 0; i < iMax; i++) {
    // ---------------------------------------------
    // Progression bar
    if ((double)(i + 1) / iMax > dProgressNext) {
      double dProgress((double)(i + 1) / iMax);

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

      dProgressNext += dProgressStep;
    }

    // ---------------------------------------------
    // Computations
    ci = CoxIter();

    ci.set_bCheckCocompactness(tests[i].bTestCompacity);
    ci.set_bCheckCofiniteness(true);

    if (!ci.bReadGraphFromFile("../../graphs/" + tests[i].szFile)) {
      of << "Error\t " << tests[i].szFile << endl;
      of << "\t\tError when reading file(" << ci.get_strError() << ")" << endl;
      iTestsSucceded["readingGraph"][1]++;
      continue;
    }

    iTestsSucceded["readingGraph"][0]++;

    iDim = ci.get_dimension();
    ci.set_dimension(0);

    // ---------------------------------------------
    // Computations
    if (!bRunTests_computations(i, &ci))
      continue;

    // ---------------------------------------------
    // Tests
    bRunTests_growth(i, &ci);

    bRunTests_signature(i, &ci, iDim);

    bRunTests_cocompactness_cofiniteness(i, &ci);

    if (tests[i].bTestArithmeticity)
      bRunTests_arithmeticity(i, &ci);

    bRunTests_Euler(i, &ci);
    bRunTests_FVector(i, &ci);
  }

  cout << "\n" << endl;

  of.close();

  bRunTests_displayInfo();

  return true;
}

vector<Test> Tests::get_tests() const { return tests; }

unsigned int Tests::get_iTestsCount() const { return tests.size(); }

string Tests::get_strError() const { return strError; }

void Tests::bRunTests_init() {
  iTestsSucceded.clear();

  iTestsSucceded["arithmeticity"] = array<unsigned int, 2>{0, 0};
  iTestsSucceded["cocompactness"] = array<unsigned int, 2>{0, 0};
  iTestsSucceded["cofiniteness"] = array<unsigned int, 2>{0, 0};
  iTestsSucceded["cofinitenessPartial"] = array<unsigned int, 2>{0, 0};
  iTestsSucceded["dimensionGuess"] = array<unsigned int, 2>{0, 0};
  iTestsSucceded["euler"] = array<unsigned int, 2>{0, 0};
  iTestsSucceded["fv"] = array<unsigned int, 2>{0, 0};
  iTestsSucceded["fvAlt"] = array<unsigned int, 2>{0, 0};
  iTestsSucceded["growthRate"] = array<unsigned int, 2>{0, 0};
  iTestsSucceded["growthSeries"] = array<unsigned int, 2>{0, 0};
  iTestsSucceded["growthSeriesDenomDimOdd"] = array<unsigned int, 2>{0, 0};
  iTestsSucceded["growthSeriesEuler"] = array<unsigned int, 2>{0, 0};
  iTestsSucceded["readingGraph"] = array<unsigned int, 2>{0, 0};
  iTestsSucceded["signature"] = array<unsigned int, 2>{0, 0};

  iTestsUnknownErrors = 0;

  strTestDescription["arithmeticity"] = "Arithmeticity";
  strTestDescription["cocompactness"] = "Cocompactness";
  strTestDescription["cofiniteness"] = "Cofiniteness";
  strTestDescription["cofinitenessPartial"] = "Cofiniteness (partial)";
  strTestDescription["dimensionGuess"] = "Dimension guess";
  strTestDescription["euler"] = "Euler charactereistic";
  strTestDescription["fv"] = "f-vector";
  strTestDescription["fvAlt"] = "Alt. sum comp. f-vector";
  strTestDescription["growthRate"] = "Growth rate";
  strTestDescription["growthSeries"] = "Growth series";
  strTestDescription["growthSeriesDenomDimOdd"] =
      "Denom. growth series vanish at 1";
  strTestDescription["growthSeriesEuler"] = "Growth series <-> Euler char.";
  strTestDescription["readingGraph"] = "Reading graph";
  strTestDescription["signature"] = "Signature";
}

bool Tests::bRunTests_computations(const unsigned int &iTestIndex,
                                   CoxIter *ci) {
  ci->exploreGraph();
  ci->computeGraphsProducts();
  if (!ci->bEulerCharacteristicFVector()) {
    of << "Error\t " << tests[iTestIndex].szFile << endl;
    of << "\t\tPlease check the graph; maybe it is not hyperbolic" << endl;
    iTestsUnknownErrors++;
    return false;
  }

  if (tests[iTestIndex].bTestCompacity)
    ci->isGraphCocompact();

  ci->checkCovolumeFiniteness();
  ci->bCanBeFiniteCovolume();
  ci->growthSeries();

  return true;
}

void Tests::bRunTests_growth(const unsigned int &iTestIndex, CoxIter *ci) {
  string strGrowthRate;
  GrowthRate *gr(new GrowthRate());
  GrowthRate_Result grr(
      gr->grrComputations(ci->get_growthSeries_denominator()));
  delete gr;
  strGrowthRate = grr.strGrowthRate;

  if (tests[iTestIndex].strGrowthRate != "") {
    if (tests[iTestIndex].strGrowthRate == strGrowthRate) {
      iTestsSucceded["growthRate"][0]++;
      of << "OK\tGrowth rate\t\t" << tests[iTestIndex].szFile << endl;
    } else {
      iTestsSucceded["growthRate"][1]++;
      bRunTests_error(iTestIndex, "growth rate",
                      tests[iTestIndex].strGrowthRate, strGrowthRate);
    }
  }

  if (grr.iPerron != 1)
    cout << "INFO: Not a Perron number in " << tests[iTestIndex].szFile << endl;

  // A small test of the growth series
  if (ci->get_isFiniteCovolume() > 0) {
    vector<mpz_class> denom;
    vector<unsigned int> cyclotomic;
    bool bReduced;

    ci->get_growthSeries(cyclotomic, denom, bReduced);

    if (tests[iTestIndex].bTestGrowthSeries) {
      if (tests[iTestIndex].growthSeries_cyclotomicNumerator == cyclotomic &&
          tests[iTestIndex].growthSeries_polynomialDenominator == denom) {
        iTestsSucceded["growthSeries"][0]++;
        of << "OK\tGrowth series\t\t" << tests[iTestIndex].szFile << endl;
      } else {
        iTestsSucceded["growthSeries"][1]++;
        of << "Error\t " << tests[iTestIndex].szFile << endl;
        of << "\t\tGrowth series" << endl;
      }
    }

    mpz_class iTotalDenom(0);
    for (unsigned int j(0); j < denom.size(); j++)
      iTotalDenom += denom[j];

    if (ci->get_dimension() % 2) // n is odd, the denominator should vanish in 1
    {
      if (iTotalDenom == 0) {
        iTestsSucceded["growthSeriesDenomDimOdd"][0]++;
        of << "OK\tGrowth series test\t" << tests[iTestIndex].szFile << endl;
      } else {
        iTestsSucceded["growthSeriesDenomDimOdd"][1]++;
        of << "Error\t " << tests[iTestIndex].szFile << endl;
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
        iTestsSucceded["growthSeriesEuler"][0]++;
        of << "OK\tGrowth series test\t" << tests[iTestIndex].szFile << endl;
      } else {
        iTestsSucceded["growthSeriesEuler"][1]++;
        of << "Error\t " << tests[iTestIndex].szFile << endl;
        of << "\t\tGrowth series test: different from Euler characteristic"
           << endl;
      }
    }
  }
}

void Tests::bRunTests_signature(const unsigned int &iTestIndex, CoxIter *ci,
                                const unsigned int &iDim) {
  int iSignatureComputed(0);

  array<unsigned int, 3> iSignature;

  try {
    if (ci->get_iHasDottedLineWithoutWeight() == 0) {
      Signature s;
      iSignature = s.iComputeSignature(ci->get_strGramMatrix_PARI());
      iSignatureComputed = 1;
    } else
      iSignatureComputed = 0;
  } catch (string strE) {
    iSignatureComputed = -1;
  }

  if (iSignatureComputed == 1 && ci->get_dimension()) {
    if (iSignature[0] != ci->get_dimension() || iSignature[1] != 1) {
      iTestsSucceded["signature"][1]++;
      of << "Error\t " << tests[iTestIndex].szFile << endl;
      of << "\t\tSignature computed: (" << iSignature[0] << "," << iSignature[1]
         << "," << iSignature[2] << ")" << endl;
    } else {
      iTestsSucceded["signature"][0]++;
      of << "OK\tSignature\t\t" << tests[iTestIndex].szFile << endl;
    }
  } else if (iSignatureComputed == -1) {
    iTestsSucceded["signature"][1]++;
    of << "Error\t " << tests[iTestIndex].szFile << endl;
    of << "\t\tSignature: format" << endl;
  } else if (iSignatureComputed == 0 &&
             tests[iTestIndex].szFile.find("roberts", 0) != string::npos) {
    iTestsSucceded["signature"][1]++;
    of << "Error\t " << tests[iTestIndex].szFile << endl;
    of << "\t\tWeights unkown" << endl;
  }

  // ---------------------------------------
  // Guessing the dimension
  if (ci->get_bDimensionGuessed()) {
    if (iDim == ci->get_dimension()) {
      iTestsSucceded["dimensionGuess"][0]++;
      of << "OK\tDimension guessed\t" << tests[iTestIndex].szFile << endl;
    } else {
      iTestsSucceded["dimensionGuess"][1]++;
      of << "Error\tDimension guessed\t" << tests[iTestIndex].szFile << endl;
    }
  }
}

void Tests::bRunTests_arithmeticity(const unsigned int &iTestIndex,
                                    CoxIter *ci) {
  int iArithmeticity;

  Arithmeticity arithmeticity;
  arithmeticity.test(*ci, false);

  iArithmeticity = ci->get_isArithmetic();
  if ((iArithmeticity == 1 && tests[iTestIndex].bIsArithmetic) ||
      (iArithmeticity == 0 && !tests[iTestIndex].bIsArithmetic)) {
    iTestsSucceded["arithmeticity"][0]++;
    of << "OK\tArithmeticity\t\t" << tests[iTestIndex].szFile << endl;
  } else {
    iTestsSucceded["arithmeticity"][1]++;
    bRunTests_error(iTestIndex, "arithmeticity",
                    tests[iTestIndex].bIsArithmetic ? "yes" : "no",
                    strIntToString(iArithmeticity));
  }
}

void Tests::bRunTests_cocompactness_cofiniteness(const unsigned int &iTestIndex,
                                                 CoxIter *ci) {
  int iCompacity, iFiniteVolume(ci->checkCovolumeFiniteness());
  bool bCanBeFiniteCovolume(ci->bCanBeFiniteCovolume());

  if (iFiniteVolume == 1 && !bCanBeFiniteCovolume) {
    iTestsSucceded["cofinitenessPartial"][1]++;
    of << "Error\t " << tests[iTestIndex].szFile << endl;
    of << "\t\tTest CanBeFiniteCovolume" << endl;
  } else if ((iFiniteVolume == 0 && !bCanBeFiniteCovolume) ||
             (iFiniteVolume == 1 && bCanBeFiniteCovolume)) {
    iTestsSucceded["cofinitenessPartial"][0]++;
    of << "OK\tCanBeFiniteCovolume\t" << tests[iTestIndex].szFile << endl;
  }

  if (!((iFiniteVolume == 1 && tests[iTestIndex].bIsFiniteVolume) ||
        (iFiniteVolume == 0 && !tests[iTestIndex].bIsFiniteVolume))) {
    iTestsSucceded["cofiniteness"][1]++;
    bRunTests_error(iTestIndex, "cofiniteness test",
                    tests[iTestIndex].bIsFiniteVolume ? "yes" : "no",
                    strIntToString(iFiniteVolume));
  } else {
    iTestsSucceded["cofiniteness"][0]++;
    of << "OK\tCofiniteness test\t" << tests[iTestIndex].szFile << endl;
  }

  if (tests[iTestIndex].bTestCompacity) {
    iCompacity = ci->get_isCocompact();
    if ((iCompacity == 1 && tests[iTestIndex].bIsCompact) ||
        (iCompacity == 0 && !tests[iTestIndex].bIsCompact)) {
      iTestsSucceded["cocompactness"][0]++;
      of << "OK\tCompacity\t\t" << tests[iTestIndex].szFile << endl;
    } else {
      iTestsSucceded["cocompactness"][1]++;
      bRunTests_error(iTestIndex, "cocompactness test",
                      tests[iTestIndex].bIsCompact ? "yes" : "no",
                      strIntToString(iCompacity));
    }
  }
}

void Tests::bRunTests_Euler(const unsigned int &iTestIndex, CoxIter *ci) {
  if (tests[iTestIndex].bTestEuler) {
    if (ci->get_brEulerCaracteristic() != tests[iTestIndex].brResult) {
      iTestsSucceded["euler"][1]++;
      bRunTests_error(iTestIndex, "Euler characteristic",
                      tests[iTestIndex].brResult.to_string(),
                      ci->get_brEulerCaracteristic().to_string());
    } else {
      iTestsSucceded["euler"][0]++;
      of << "OK\tEuler\t\t\t" << tests[iTestIndex].szFile << endl;
    }
  } else if (ci->get_dimension() % 2 && tests[iTestIndex].bIsFiniteVolume) {
    if (ci->get_brEulerCaracteristic() != 0) {
      iTestsSucceded["euler"][1]++;
      bRunTests_error(iTestIndex, "Euler characteristic", "0",
                      ci->get_brEulerCaracteristic().to_string());
    } else {
      iTestsSucceded["euler"][0]++;
      of << "OK\tEuler\t\t\t" << tests[iTestIndex].szFile << endl;
    }
  }
}

void Tests::bRunTests_FVector(const unsigned int &iTestIndex, CoxIter *ci) {
  if (ci->get_isFiniteCovolume() > 0) {
    if (ci->get_iFVectorAlternateSum() == (ci->get_dimension() % 2 ? 2 : 0)) {
      iTestsSucceded["fvAlt"][0]++;
      of << "OK\tAlt. sum\t\t" << tests[iTestIndex].szFile << endl;
    } else {
      iTestsSucceded["fvAlt"][1]++;
      bRunTests_error(iTestIndex, "alt. sum of components of f-vector",
                      to_string(ci->get_dimension() % 2 ? 2 : 0),
                      to_string(ci->get_iFVectorAlternateSum()));
    }
  }

  if (tests[iTestIndex].bTestFVector) {
    if (ci->get_iFVector() == tests[iTestIndex].iFVector) {
      iTestsSucceded["fv"][0]++;
      of << "OK\tf-vector\t\t" << tests[iTestIndex].szFile << endl;
    } else {
      iTestsSucceded["fv"][1]++;
      bRunTests_error(iTestIndex, "f-vector",
                      implode(",", tests[iTestIndex].iFVector),
                      implode(",", ci->get_iFVector()));
    }
  }
}

void Tests::bRunTests_error(const unsigned int &iTestIndex,
                            const string &strTest, const string &strExpected,
                            const string &strComputed) {
  of << "Error: " << strTest << "\t in " << tests[iTestIndex].szFile << endl;
  of << "\t\tExpected: " << strExpected << endl;
  of << "\t\tComputed: " << strComputed << endl;
}

string Tests::strIntToString(const int &i) {
  if (i == 0)
    return "no";
  else if (i == 1)
    return "yes";
  else
    return "?";
}

void Tests::bRunTests_displayInfo() {
  array<unsigned int, 2> iTotal({0, 0});

  cout << "------------------------------------------------------" << endl;
  cout << setw(32) << left << "Test"
       << " | #success | #error |" << endl;
  cout << "------------------------------------------------------" << endl;

  for (map<string, string>::const_iterator it(strTestDescription.begin());
       it != strTestDescription.end(); it++) {
    iTotal[0] += iTestsSucceded[it->first][0];
    iTotal[1] += iTestsSucceded[it->first][1];

    cout << setw(32) << left << it->second << " | ";
    cout << setw(8) << left << iTestsSucceded[it->first][0] << " | ";
    cout << setw(6) << left << iTestsSucceded[it->first][1] << " | ";
    cout << endl;
  }

  cout << "------------------------------------------------------" << endl;
  cout << setw(32) << left << "Total"
       << " | ";
  cout << setw(8) << left << iTotal[0] << " | ";
  cout << setw(6) << left << iTotal[1] << " | ";
  cout << endl;
  cout << "------------------------------------------------------\n" << endl;

  if (iTestsUnknownErrors)
    cout << "Other errors: " << iTestsUnknownErrors << endl;

  cout << "For more information, see:\n\t" << strInputFilename << ".output"
       << "\n"
       << endl;
}
