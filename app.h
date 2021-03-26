/*
Copyright (C) 2013-2017
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

/*!
 * \file app.h
 * \author Rafael Guglielmetti
 *
 * \class App
 * \brief Main class for the application
 */

#ifndef APP_H
#define APP_H

#include <algorithm>
#include <chrono>
#include <string>
#include <vector>

#ifdef WIN32
#include <io.h>
#else
#include <unistd.h>
#endif

#ifdef _COMPILE_WITH_PARI_
#include "growthrate.h"
#include "signature.h"
#endif

using namespace std;

#include "arithmeticity.h"
#include "coxiter.h"
#include "index2.h"

class App {
private:
  bool checkArithmeticity; ///< If we want to whether the group is arithmetic
                           ///< or not
  bool checkCanBeFiniteCovolume; ///< If we want to check whether the group can
                                 ///< be of finite volume or not
  bool checkCocompacity;         ///< If we want to check whether the group is
                                 ///< cocompact or not
  bool checkFiniteCovolume;      ///< If we want to check whether the group has
                                 ///< finite covolume or not
  bool computeGrowthSeries;      ///< If we want compute the growth series
  bool computeGrowthRate;        ///< If we want to compute the growth rate
  bool computeEuler;       ///< If we want to compute the Euler characteristic
  bool computeSignature;   ///< If we want to compute the signature
  bool debug;              ///< Display additional information
  bool bIndex2;            ///< Trying to extract an index two subroup?
  bool useOpenMP;          ///< Use OpenMP
  bool printCoxeterGraph;  ///< Print the Coxeter graph?
  bool printCoxeterMatrix; ///< Print the Coxeter matrix?
  bool printGramMatrix;    ///< Print the Gram matrix?
  bool bPrintHelp; ///< If we want to print help (option or by default depending
                   ///< on the error)
  string ouputMathematicalFormat;  ///< Format of output: generic,
                                   ///< mathematica, pari
  vector<string> verticesToRemove; ///< The vertices we want to remove
  vector<string> vertices;         ///< If we specify the vertices
  string index2vertex_t0; ///< Reflexion and glueing with respect to the
                          ///< hyperplane corresponding to t0
  string index2vertex_s0; ///< If we want to compute the f-vector of the
                          ///< corresponding infinite sequence

public:
  bool bCoutFile;          ///< If the output is redirected to a file
  bool bOutputGraphToDraw; ///< If we write the graph in a file, to use graphviz
  bool bOutputGraph; ///< If we write the graph in a file, to use with CoxIter
  string inFilename; ///< Path to the graph
  string outFilenameBasis; ///< Path to the file for the output (+ .output,
                           ///< .graphviz, .coxiter)

public:
  App();

  bool readMainParameters(int argc, char **argv);
  void run();
  void printHelp() const;

private:
  void extractIndex2Subgroup(CoxIter &ci);
};

#endif // APP_H
