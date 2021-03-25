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

void Arithmeticity::test(CoxIter &ci_, const bool &bListCycles_) {
  ci = &ci_;
  listCycles = bListCycles_;
  strListCycles.clear();
  if (ci->get_isCocompact() !=
      0) // If the graph is cocompact (1) or if we don't know (-1)
  {
    strError = "GROUP COCOMPACTNESS";
    ci->set_isArithmetic(-1);
    return;
  }

  size_t i, j;

  coxeterMatrix = ci->get_coxeterMatrix();

  verticesCount = ci->get_verticesCount();
  iReferencesToLabels.clear();
  notArithmetic = false;

  // ------------------------------------------------------
  // Cycles consisting of two elements
  bool bSQRT2(false), bSQRT3(false);

  for (i = 0; i < verticesCount; i++) {
    iReferencesToLabels.push_back(i);

    for (j = 0; j < i; j++) {
      if (i == j)
        continue;

      if (coxeterMatrix[i][j] != 1 && coxeterMatrix[i][j] != 0 &&
          coxeterMatrix[i][j] != 2 && coxeterMatrix[i][j] != 3 &&
          coxeterMatrix[i][j] != 4 && coxeterMatrix[i][j] != 6) {
        if (ci->get_bDebug())
          cout << "\tNot arithmetic: 2*G(" << ci->get_vertexLabel(i) << ","
               << ci->get_vertexLabel(j) << ") = pi/" << coxeterMatrix[i][j]
               << endl;

        ci->set_isArithmetic(0);
        return;
      }

      if (coxeterMatrix[i][j] == 4)
        bSQRT2 = true;
      else if (coxeterMatrix[i][j] == 6)
        bSQRT3 = true;
      else if (coxeterMatrix[i][j] == 1) {
        string strC("4 * " + string("l") + to_string(j) + "m" + to_string(i) +
                    "^2");

        auto it(lower_bound(strListCycles.begin(), strListCycles.end(), strC));
        if (it == strListCycles.end() || *it != strC)
          strListCycles.insert(it, strC);
      }
    }
  }

  // ------------------------------------------------------
  // Here, we know that m_{ij} \in {2,3,4,6,infty}

  // If true, the group is arithmetic
  if (!bSQRT2 && !bSQRT3 && !ci->get_bHasDottedLine()) {
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
  unsigned int iPreviousVertex, iCurrrentVertex;

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
      for (iCurrrentVertex = 0; iCurrrentVertex < verticesCount;
           iCurrrentVertex++) {
        if (coxeterMatrix[i][iCurrrentVertex] != 2 &&
            neighboursCount[iCurrrentVertex] == 2) {
          verticesToRemove.push_back(iCurrrentVertex);
          break;
        }
      }

      if (iCurrrentVertex ==
          verticesCount) // If we cannot remove more than one vertex
        continue;

      iPreviousVertex = i;

      // Here: i ------- iCurrentVertex

      // We continue the path
      while (true) {
        for (k = 0; k < verticesCount; k++) {
          if (coxeterMatrix[iCurrrentVertex][k] != 2 &&
              neighboursCount[k] == 2 && k != iPreviousVertex) {
            verticesToRemove.push_back(k);
            iPreviousVertex = iCurrrentVertex;
            iCurrrentVertex = k;
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
    iReferencesToLabels.erase(iReferencesToLabels.begin() + *it);

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
    bVerticesVisited = vector<bool>(verticesCount, false);
    bEdgesVisited =
        vector<vector<bool>>(verticesCount, vector<bool>(verticesCount, false));
    findCycles(i, i);

    if (notArithmetic) {
      ci->set_isArithmetic(0);
      return;
    }
  }

  if (!ci->get_bHasDottedLine())
    ci->set_isArithmetic(1);
  else
    ci->set_isArithmetic(-1);
}

void Arithmeticity::findCycles(const unsigned int &iRoot,
                               const unsigned int &iFrom) {
  path.push_back(iRoot); // We add the vertex to the path

  for (unsigned int i(path[0]); i < verticesCount; i++) {
    // If i is a neighbour and if we did not visit this edge
    if (coxeterMatrix[iRoot][i] != 2 && !bEdgesVisited[iRoot][i]) {
      if (i == path[0]) {
        if (path[1] <
            path[path.size() - 1]) // We do not to test each cycle twice
          testCycle();

        if (notArithmetic) {
          if (ci->get_bDebug()) {
            cout << "\tNot arithmetic\n\t\tCycle: ";
            for (vector<unsigned int>::const_iterator it(path.begin());
                 it != path.end();
                 ++it) // We display the components of the cycle
              cout << (it == path.begin() ? "" : ", ")
                   << ci->get_vertexLabel(iReferencesToLabels[*it]);
            cout << endl;
          }

          return;
        }
      } else if (find(path.begin(), path.end(), i) == path.end()) {
        bEdgesVisited[iRoot][i] = bEdgesVisited[i][iRoot] = true;
        findCycles(i, iRoot);

        if (notArithmetic)
          return;
      }
    }
  }

  if (iFrom != iRoot)
    bEdgesVisited[iRoot][iFrom] = bEdgesVisited[iFrom][iRoot] = false;

  path.pop_back();
}

void Arithmeticity::testCycle() {
  unsigned int pathSize(path.size());

  if (!listCycles) {
    bool bNumberSQRT2Even(
        coxeterMatrix[path[0]][path[pathSize - 1]] == 4 ? false : true),
        bNumberSQRT3Even(
            coxeterMatrix[path[0]][path[pathSize - 1]] == 6 ? false : true);

    for (unsigned int i(1); i < pathSize; i++) {
      if (coxeterMatrix[path[i]][path[i - 1]] == 4)
        bNumberSQRT2Even = bNumberSQRT2Even ? false : true;
      else if (coxeterMatrix[path[i]][path[i - 1]] == 6)
        bNumberSQRT3Even = bNumberSQRT3Even ? false : true;
      else if (coxeterMatrix[path[i]][path[i - 1]] ==
               1) // Because of the dotted line we cannot say anything for this
                  // cycle
        return;
    }

    if (coxeterMatrix[path[0]][path[pathSize - 1]] ==
        1) // Because of the dotted line we cannot say anything for this cycle
      return;

    notArithmetic = !bNumberSQRT2Even || !bNumberSQRT3Even;
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
      string strTemp;

      twoCount += dottedCount;

      if (sqrt2Count > 1)
        twoCount += ((sqrt2Count % 2) ? sqrt2Count - 1 : sqrt2Count) / 2;

      if (twoCount)
        strTemp +=
            (strTemp == "" ? "" : " * ") + string("2^") + to_string(twoCount);

      if (sqrt3Count > 1)
        strTemp +=
            (strTemp == "" ? "" : " * ") + string("3^") +
            to_string(((sqrt2Count % 2) ? sqrt2Count - 1 : sqrt2Count) / 2);

      if (sqrt2Count % 2)
        strTemp += (strTemp == "" ? "" : " * ") + string("Sqrt[2]");

      if (sqrt3Count % 2)
        strTemp += (strTemp == "" ? "" : " * ") + string("Sqrt[3]");

      for (unsigned int i(1); i < pathSize; i++) {
        if (coxeterMatrix[path[i]][path[i - 1]] == 1)
          strTemp += (strTemp == "" ? "" : " * ") + string("l") +
                     to_string(min(iReferencesToLabels[path[i]],
                                   iReferencesToLabels[path[i - 1]])) +
                     "m" +
                     to_string(max(iReferencesToLabels[path[i]],
                                   iReferencesToLabels[path[i - 1]]));
      }

      if (coxeterMatrix[path[0]][path[pathSize - 1]] == 1)
        strTemp += (strTemp == "" ? "" : " * ") + string("l") +
                   to_string(min(iReferencesToLabels[path[0]],
                                 iReferencesToLabels[path[pathSize - 1]])) +
                   "m" +
                   to_string(max(iReferencesToLabels[path[0]],
                                 iReferencesToLabels[path[pathSize - 1]]));

      auto it(lower_bound(strListCycles.begin(), strListCycles.end(), strTemp));
      if (it == strListCycles.end() || *it != strTemp)
        strListCycles.insert(it, strTemp);
    }
  }
}

vector<string> Arithmeticity::get_strListCycles() { return strListCycles; }

string Arithmeticity::get_strError() { return strError; }
