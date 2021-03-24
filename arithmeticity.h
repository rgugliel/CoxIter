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
 * \file arithmeticity.h
 * \author Rafael Guglielmetti
 *
 * \class Arithmeticity
 * \brief This class tests the arithmeticity of a graph which has no dotted edge
 * and which is non-cocompact. It uses Vinberg's criteria.
 */

#ifndef ARITHMETICITY_H
#define ARITHMETICITY_H

#include "coxiter.h"

class Arithmeticity {
private:
  string strError; ///< If an error occured, small text.

  CoxIter *ci;                 ///< Pointer to the CoxIter object
  unsigned int verticesCount; ///< Number of generators of the group
  vector<vector<unsigned int>> coxeterMatrix; ///< Coxeter matrix of the group
  vector<unsigned int> iReferencesToLabels;    ///< Correspondence for the new
                                               ///< indices to the old ones

  // For the DFS
  vector<vector<bool>> bEdgesVisited; ///< Traversed edges
  vector<bool> bVerticesVisited;      ///<  Taversed vertices
  vector<unsigned int> path;         ///< Current path

  bool notArithmetic; ///< True if not arithmetic (i.e. we have to quit the
                       ///< algorithm)

  bool listCycles; ///< If true, will list the cycles to be manually tested
  vector<string> strListCycles; ///< The list

public:
  /*! \fn Arithmeticity()
   * 	\brief Basic constructor
   */
  Arithmeticity();

  /*! 	\fn ~Arithmeticity()
   * 	\brief Destructor
   */
  ~Arithmeticity();

  /*! 	\fn test
   *	\brief Test the arithmeticity of a graph
   *
   * 	\param ci(CoxIter&) The graph
   * 	\param bListCycles_(const bool&) If true, will list the cycles to be
   *manually tested \return True if success, false otherwise. Then, use
   *ci.get_isArithmetic()
   */
  void test(CoxIter &ci, const bool &bListCycles_);

  /*!	\fn get_strListCycles
   *	\brief Return the list of cycles
   *
   *	\return List (vector<string>)
   */
  vector<string> get_strListCycles();

  /*!	\fn get_strError
   *	\brief Return the error code
   *
   *	\return Error code (string)
   */
  string get_strError();

private:
  /*!	\fn collapseQueues
   * 	\brief Try to collapse queues of the graph
   *
   * 	\return Number of vertices removed
   */
  unsigned int collapseQueues();

  /*! 	\fn testCycles
   * 	\brief Test the cycles
   */
  void testCycles();

  /*! 	\fn findCycles
   * 	\brief Look for cycles
   *
   * 	Update the vector path to find cycles
   *
   * 	\param iRoot(unsigned int&) Starting vertex
   * 	\param iFrom(unsigned int&) Previoud vertex (if recursive call); iRoot
   * otherwise
   */
  void findCycles(const unsigned int &root, const unsigned int &from);

  /*! 	\fn testCycle
   * 	\brief Test the cycle in path
   *
   * 	This function is called by findCycles. Eventually, set bNotArithmetic to
   * true
   */
  void testCycle();
};

struct CycleElement {
  unsigned int first;
  unsigned int second;
};

struct Cycles {};

#endif // ARITHMETICITY_H
