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

#include "index2.h"

Index2::Index2(CoxIter *ci)
    : ci(ci), verticesCount(ci->get_verticesCount()), iNewVerticesCount(0),
      coxeterMatrix(ci->get_coxeterMatrix()) {}

bool Index2::isVertexAdmissible(const string &vertexLabel) {

  if (!ci->isVertexValid(vertexLabel)) {
    error = "Invalid vertex name: " + vertexLabel;
    return false;
  }

  unsigned int index(ci->get_vertexIndex(vertexLabel));

  for (unsigned int i(0); i < verticesCount; i++) {
    if (coxeterMatrix[index][i] != 0 && coxeterMatrix[index][i] != 1 &&
        (coxeterMatrix[index][i] % 2)) {
      error =
          "m(" + vertexLabel + "," + ci->get_vertexLabel(i) + ") is not even";
      return false;
    }
  }

  return true;
}

bool Index2::removeVertex(const string &vertexName) {
  unsigned int weight;

  // -------------------------------------------------------
  // some verifications, firsts constructions
  if (!isVertexAdmissible(vertexName))
    return false;

  vertex = ci->get_vertexIndex(vertexName);

  // Remark: we don't merge the previous and the next loop. Thus, if there is a
  // problem, we do not change the CoxIter object.

  // -------------------------------------------------------
  // new vertices
  for (unsigned int i(0); i < verticesCount; i++) {
    // If they don't commute, this will add a new vertex
    if (coxeterMatrix[vertex][i] != 2) {
      NewVertex nv;
      nv.index = verticesCount + iNewVerticesCount - 1;

      nv.originVertex = i;
      nv.label = ci->get_vertexLabel(vertex) + "_" + ci->get_vertexLabel(i);

      newVertices.push_back(nv);
      ci->map_vertices_labels_addReference(nv.label);

      iNewVerticesCount++;
    }
  }

  ci->map_vertices_labels_removeReference(vertex);

  // -------------------------------------------------------
  // new adjacency matrix
  for (unsigned int i(0); i <= (verticesCount + iNewVerticesCount - 1);
       i++) // the -1 is for the vertex we removed
  {
    if (i == vertex)
      continue;

    iNewCox.push_back(
        vector<unsigned int>(verticesCount + iNewVerticesCount - 1,
                             2)); // TODO: tout construire en une fois
  }

  // we fill the old values (upper left square)
  for (unsigned int i(0); i < verticesCount; i++) {
    if (i == vertex)
      continue;

    for (unsigned int j(0); j <= i; j++)
      iNewCox[i > vertex ? i - 1 : i][j > vertex ? j - 1 : j] =
          iNewCox[j > vertex ? j - 1 : j][i > vertex ? i - 1 : i] =
              coxeterMatrix[i][j];
  }

  // -------------------------------------------------------
  // weights of the new edges
  for (vector<NewVertex>::const_iterator itNew(newVertices.begin());
       itNew != newVertices.end();
       ++itNew) // Going through the set of new vertices
  {
    // A new vertex is of the form: t0 * sj * t0, where sj =
    // NewVertex.iOriginVertex

    // for every old vertex
    for (unsigned int i(0); i < verticesCount;
         i++) // Goal: find m(t0 * sj * t0, si)
    {
      if (i == vertex) // we removed this vertex
        continue;

      if ((*itNew).originVertex == i) // m(t0 * si * t0, si) = m(t0, si) / 2
        weight = coxeterMatrix[vertex][i] == 0 || coxeterMatrix[vertex][i] == 1
                     ? coxeterMatrix[vertex][i]
                     : coxeterMatrix[vertex][i] / 2;
      else {
        if (coxeterMatrix[i][vertex] ==
            2) // If si commutes with t0, then m(t0 * sj * t0, si) = m(sj, si)
          weight = coxeterMatrix[(*itNew).originVertex][i];
        else if (coxeterMatrix[(*itNew).originVertex][vertex] != 2) {
          if (coxeterMatrix[vertex][i] == 4 &&
              coxeterMatrix[vertex][(*itNew).originVertex] == 4 &&
              coxeterMatrix[i][(*itNew).originVertex] == 2)
            weight = 0;
          else
            weight = 1;
        } else
          throw(string("Index2::removeVertex: Error")); // TODO: idem
      }

      iNewCox[i > vertex ? i - 1 : i][(*itNew).index] = weight;
      iNewCox[(*itNew).index][i > vertex ? i - 1 : i] = weight;
    }

    // for every new vertex t0 * si * t0
    for (vector<NewVertex>::const_iterator itNewSub(newVertices.begin());
         itNewSub != itNew; ++itNewSub) {
      iNewCox[(*itNew).index][(*itNewSub).index] =
          coxeterMatrix[(*itNew).originVertex][(*itNewSub).originVertex];
      iNewCox[(*itNewSub).index][(*itNew).index] =
          coxeterMatrix[(*itNew).originVertex][(*itNewSub).originVertex];
    }
  }

  ci->set_coxeterMatrix(iNewCox);

  return true;
}

void Index2::printMatrix(vector<vector<unsigned int>> *iMatrix) {
  for (vector<vector<unsigned int>>::const_iterator itRow(iMatrix->begin());
       itRow != iMatrix->end(); ++itRow) {
    for (vector<unsigned int>::const_iterator itCol((*itRow).begin());
         itCol != (*itRow).end(); ++itCol)
      cout << *itCol << " ";

    cout << endl;
  }
}

string Index2::get_error() const { return error; }
