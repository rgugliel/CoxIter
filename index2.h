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
 * \file index2.h
 * \author Rafael Guglielmetti
 *
 * \class Index2
 * \brief Extract an index two subroup (based on Gal, Bonnafé-Dyer)
 */

#ifndef __INDEX2_H__
#define __INDEX2_H__

#include <string>

using namespace std;

#include "coxiter.h"

struct NewVertex {
  string label;
  unsigned int index;
  unsigned int originVertex;
};

class Index2 {
private:
  CoxIter *ci;
  vector<vector<unsigned int>> coxeterMatrix; ///< Coxeter matrix
  unsigned int verticesCount; ///< Number of vertices of the starting graph

  unsigned int vertex; ///< Index of the vertex

  vector<NewVertex> newVertices;
  vector<vector<unsigned int>> iNewCox; ///< New Coxeter matrix
  unsigned int iNewVerticesCount;

  string error;

public:
  Index2(CoxIter *ci);

  bool isVertexAdmissible(const string &vertexLabel);
  bool removeVertex(const string &vertexLabel);
  void printMatrix(vector<vector<unsigned int>> *iMatrix);

  string get_error() const;
};

#endif // __INDEX2_H__
