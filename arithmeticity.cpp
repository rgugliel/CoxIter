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

#include "arithmeticity.h"

Arithmeticity::Arithmeticity()
    : verticesCount(0), notArithmetic(false), ci(0), listCycles(false) {}

Arithmeticity::~Arithmeticity() {}

void Arithmeticity::test(CoxIter &ci_, const bool &listCycles_) {
  ci = &ci_;
  listCycles = listCycles_;
  allCycles.clear();
  if (ci->get_isCocompact() !=
      0) // If the graph is cocompact (1) or if we don't know (-1)
  {
    error = "GROUP COCOMPACTNESS";
    ci->set_isArithmetic(-1);
    return;
  }

  size_t i, j;

  coxeterMatrix = ci->get_coxeterMatrix();

  verticesCount = ci->get_verticesCount();
  referencesToLabels.clear();
  notArithmetic = false;

  // ------------------------------------------------------
  // Cycles consisting of two elements
  bool isSqrt2(false), isSqrt3(false);

  for (i = 0; i < verticesCount; i++) {
    referencesToLabels.push_back(i);

    for (j = 0; j < i; j++) {
      if (i == j)
        continue;

      if (coxeterMatrix[i][j] != 1 && coxeterMatrix[i][j] != 0 &&
          coxeterMatrix[i][j] != 2 && coxeterMatrix[i][j] != 3 &&
          coxeterMatrix[i][j] != 4 && coxeterMatrix[i][j] != 6) {
        if (ci->get_debug())
          cout << "\tNot arithmetic: 2*G(" << ci->get_vertexLabel(i) << ","
               << ci->get_vertexLabel(j) << ") = pi/" << coxeterMatrix[i][j]
               << endl;

        ci->set_isArithmetic(0);
        return;
      }

      if (coxeterMatrix[i][j] == 4)
        isSqrt2 = true;
      else if (coxeterMatrix[i][j] == 6)
        isSqrt3 = true;
      else if (coxeterMatrix[i][j] == 1) {
        string cycle("4 * " + string("l") + to_string(j) + "m" + to_string(i) +
                     "^2");

        auto it(lower_bound(allCycles.begin(), allCycles.end(), cycle));
        if (it == allCycles.end() || *it != cycle)
          allCycles.insert(it, cycle);
      }
    }
  }

  // ------------------------------------------------------
  // Here, we know that m_{ij} \in {2,3,4,6,infty}

  // If true, the group is arithmetic
  if (!isSqrt2 && !isSqrt3 && !ci->get_hasDottedLine()) {
    ci->set_isArithmetic(1);
    return;
  }

  while (collapseQueues() && verticesCount)
    ;

  if (verticesCount <= 2) {
    ci->set_isArithmetic(1);
    return;
  }

  testCycles();

  return;
}

unsigned int Arithmeticity::collapseQueues() {
  vector<unsigned int> verticesToRemove;
  vector<unsigned int> neighboursCount(
      verticesCount, 0); // For each vertex, number of neighbours
  unsigned int i, j, k;
  unsigned int previousVertex, currrentVertex;

  // ------------------------------------------------------
  // We count the number of neighbours of each vertex
  for (i = 0; i < verticesCount; i++) {
    for (j = 0; j < verticesCount; j++) {
      if (coxeterMatrix[i][j] != 2)
        neighboursCount[i]++;
    }
  }

  // ------------------------------------------------------
  // We determine wich vertices we have to remove
  for (i = 0; i < verticesCount; i++) {
    if (neighboursCount[i] == 1) // If this is a queue
    {
      verticesToRemove.push_back(i);

      // vertex next to the queue
      for (currrentVertex = 0; currrentVertex < verticesCount;
           currrentVertex++) {
        if (coxeterMatrix[i][currrentVertex] != 2 &&
            neighboursCount[currrentVertex] == 2) {
          verticesToRemove.push_back(currrentVertex);
          break;
        }
      }

      if (currrentVertex ==
          verticesCount) // If we cannot remove more than one vertex
        continue;

      previousVertex = i;

      // Here: i ------- iCurrentVertex

      // We continue the path
      while (true) {
        for (k = 0; k < verticesCount; k++) {
          if (coxeterMatrix[currrentVertex][k] != 2 &&
              neighboursCount[k] == 2 && k != previousVertex) {
            verticesToRemove.push_back(k);
            previousVertex = currrentVertex;
            currrentVertex = k;
            break;
          }
        }

        if (k == verticesCount) // The end of the path
          break;
      }
    }
  }

  // ------------------------------------------------------
  // We remove the vertices
  sort(verticesToRemove.begin(), verticesToRemove.end());
  verticesToRemove = vector<unsigned int>(
      verticesToRemove.begin(),
      unique(verticesToRemove.begin(), verticesToRemove.end()));

  if (verticesToRemove.size() == verticesCount) // If we remove all the vertices
  {
    coxeterMatrix = vector<vector<unsigned int>>(0, vector<unsigned int>(0));
    verticesCount = 0;

    return verticesToRemove.size();
  }

  reverse(verticesToRemove.begin(), verticesToRemove.end());
  for (vector<unsigned int>::const_iterator it(verticesToRemove.begin());
       it != verticesToRemove.end(); ++it) {
    coxeterMatrix.erase(coxeterMatrix.begin() + *it);
    referencesToLabels.erase(referencesToLabels.begin() + *it);

    for (vector<vector<unsigned int>>::iterator itRow(coxeterMatrix.begin());
         itRow != coxeterMatrix.end(); ++itRow)
      itRow->erase(itRow->begin() + *it);
  }

  verticesCount -= verticesToRemove.size();
  return verticesToRemove.size();
}

void Arithmeticity::testCycles() {
  for (unsigned int i(0); i < verticesCount; i++) {
    path.clear();
    visitedVertices = vector<bool>(verticesCount, false);
    visitedEdges =
        vector<vector<bool>>(verticesCount, vector<bool>(verticesCount, false));
    findCycles(i, i);

    if (notArithmetic) {
      ci->set_isArithmetic(0);
      return;
    }
  }

  if (!ci->get_hasDottedLine())
    ci->set_isArithmetic(1);
  else
    ci->set_isArithmetic(-1);
}

void Arithmeticity::findCycles(const unsigned int &root,
                               const unsigned int &from) {
  path.push_back(root); // We add the vertex to the path

  for (unsigned int i(path[0]); i < verticesCount; i++) {
    // If i is a neighbour and if we did not visit this edge
    if (coxeterMatrix[root][i] != 2 && !visitedEdges[root][i]) {
      if (i == path[0]) {
        if (path[1] <
            path[path.size() - 1]) // We do not to test each cycle twice
          testCycle();

        if (notArithmetic) {
          if (ci->get_debug()) {
            cout << "\tNot arithmetic\n\t\tCycle: ";
            for (vector<unsigned int>::const_iterator it(path.begin());
                 it != path.end();
                 ++it) // We display the components of the cycle
              cout << (it == path.begin() ? "" : ", ")
                   << ci->get_vertexLabel(referencesToLabels[*it]);
            cout << endl;
          }

          return;
        }
      } else if (find(path.begin(), path.end(), i) == path.end()) {
        visitedEdges[root][i] = visitedEdges[i][root] = true;
        findCycles(i, root);

        if (notArithmetic)
          return;
      }
    }
  }

  if (from != root)
    visitedEdges[root][from] = visitedEdges[from][root] = false;

  path.pop_back();
}

void Arithmeticity::testCycle() {
  unsigned int pathSize(path.size());

  if (!listCycles) {
    bool isSqrt2CountEven = coxeterMatrix[path[0]][path[pathSize - 1]] != 4;
    bool isSqrt3CountEven = coxeterMatrix[path[0]][path[pathSize - 1]] != 6;

    for (unsigned int i(1); i < pathSize; i++) {
      if (coxeterMatrix[path[i]][path[i - 1]] == 4)
        isSqrt2CountEven = isSqrt2CountEven ? false : true;
      else if (coxeterMatrix[path[i]][path[i - 1]] == 6)
        isSqrt3CountEven = isSqrt3CountEven ? false : true;
      else if (coxeterMatrix[path[i]][path[i - 1]] ==
               1) // Because of the dotted line we cannot say anything for this
                  // cycle
        return;
    }

    // Because of the dotted line we cannot say anything for this cycle
    if (coxeterMatrix[path[0]][path[pathSize - 1]] == 1)
      return;

    notArithmetic = !isSqrt2CountEven || !isSqrt3CountEven;
  } else {
    unsigned int pathSize(path.size());

    int dottedCount(coxeterMatrix[path[0]][path[pathSize - 1]] == 1 ? 1 : 0);
    int twoCount(coxeterMatrix[path[0]][path[pathSize - 1]] == 0 ? 1 : 0);
    int sqrt2Count(coxeterMatrix[path[0]][path[pathSize - 1]] == 4 ? 1 : 0);
    int sqrt3Count(coxeterMatrix[path[0]][path[pathSize - 1]] == 6 ? 1 : 0);

    for (unsigned int i(1); i < pathSize; i++) {
      if (coxeterMatrix[path[i]][path[i - 1]] == 0)
        twoCount++;
      else if (coxeterMatrix[path[i]][path[i - 1]] == 4)
        sqrt2Count++;
      else if (coxeterMatrix[path[i]][path[i - 1]] == 6)
        sqrt3Count++;
      else if (coxeterMatrix[path[i]][path[i - 1]] ==
               1) // Because of the dotted line we cannot say anything for this
                  // cycle
        dottedCount++;
    }

    if (dottedCount == 0)
      notArithmetic = (sqrt2Count % 2) || (sqrt3Count % 2);
    else {
      string temp;

      twoCount += dottedCount;

      if (sqrt2Count > 1)
        twoCount += ((sqrt2Count % 2) ? sqrt2Count - 1 : sqrt2Count) / 2;

      if (twoCount)
        temp += (temp == "" ? "" : " * ") + string("2^") + to_string(twoCount);

      if (sqrt3Count > 1)
        temp += (temp == "" ? "" : " * ") + string("3^") +
                to_string(((sqrt2Count % 2) ? sqrt2Count - 1 : sqrt2Count) / 2);

      if (sqrt2Count % 2)
        temp += (temp == "" ? "" : " * ") + string("Sqrt[2]");

      if (sqrt3Count % 2)
        temp += (temp == "" ? "" : " * ") + string("Sqrt[3]");

      for (unsigned int i(1); i < pathSize; i++) {
        if (coxeterMatrix[path[i]][path[i - 1]] == 1)
          temp += (temp == "" ? "" : " * ") + string("l") +
                  to_string(min(referencesToLabels[path[i]],
                                referencesToLabels[path[i - 1]])) +
                  "m" +
                  to_string(max(referencesToLabels[path[i]],
                                referencesToLabels[path[i - 1]]));
      }

      if (coxeterMatrix[path[0]][path[pathSize - 1]] == 1)
        temp += (temp == "" ? "" : " * ") + string("l") +
                to_string(min(referencesToLabels[path[0]],
                              referencesToLabels[path[pathSize - 1]])) +
                "m" +
                to_string(max(referencesToLabels[path[0]],
                              referencesToLabels[path[pathSize - 1]]));

      const auto it(lower_bound(allCycles.begin(), allCycles.end(), temp));
      if (it == allCycles.end() || *it != temp)
        allCycles.insert(it, temp);
    }
  }
}

vector<string> Arithmeticity::get_allCycles() { return allCycles; }

string Arithmeticity::get_error() { return error; }
