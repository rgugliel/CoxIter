/*
Copyright (C) 2013, 2014, 2015, 2016
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

#include <string>
#include <vector>
#include <chrono>
#include <algorithm>

#ifdef _COMPILE_WITH_PARI_
#include "growthrate.h"
#include "signature.h"
#endif

// cmake -DCOMPILE_WITH_PARI=1 ../ TODO

using namespace std;

#include "coxiter.h"
#include "gbd.h"
#include "arithmeticity.h"

class App
{
	private:
		bool bCheckArithmeticity; ///< If we want to whether the group is arithmetic or not
		bool bCheckCanBeFiniteCovolume; ///< If we want to check whether the group can be of finite volume or not
		bool bCheckCocompacity; ///< If we want to check whether the group is cocompact or not
		bool bCheckFiniteCovolume; ///< If we want to check whether the group has finite covolume or not
		bool bComputeGrowthSeries; ///< If we want compute the growth series
		bool bComputeGrowthRate; ///< If we want to compute the growth rate
		bool bComputeEuler; ///< If we want to compute the Euler characteristic
		bool bComputeSignature; ///< If we want to compute the signature
		bool bDoComputations; ///< If we use the big integers GMPlib) library or not
		bool bDebug; ///< Display additional information
		bool bGBD;
		bool bOpenMP; ///< Use OpenMP
		bool bPrintCoxeterMatrix; ///< Print the Coxeter matrix?
		bool bPrintGramMatrix; ///< Print the Gram matrix?
		string strOuputMathematicalFormat; ///< Format of output: generic, mathematica, pari
		vector< string > strVerticesRemove; ///< The vertices we want to remove
		vector< string > strVertices; ///< If we specify the vertices
		string strGBDvertex;
		
	public:
		bool bCoutFile; ///< If the output is redirected to a file
		bool bOutputGraphToDraw; ///< If we write the graph in a file, to use graphviz
		bool bOutputGraph; ///< If we write the graph in a file, to use with CoxIter
		string strInFilename; ///< Path to the graph
		string strOutFilenameBasis; ///< Path to the file for the output (+ .output, .graphviz, .coxiter)
		
	public:
		App( );
		
		bool bReadMainParameters( int argc, char **argv );
		void run( );
};

#endif // APP_H
