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

#include "graph.h"

Graph::Graph(const vector<short unsigned int> &vertices,
             vector<string> *ptr_map_vertices_indexToLabel,
             const vector<bool> &bVerticesLinkable, const unsigned int &iType,
             const bool &isSpherical, const unsigned int &iDataSupp)
    : type(iType), vertices(vertices), bVerticesLinkable(bVerticesLinkable),
      dataSupp(iDataSupp), isSpherical(isSpherical),
      ptr_map_vertices_indexToLabel(ptr_map_vertices_indexToLabel),
      b_map_vertices_indexToLabelIsEmpty(
          !ptr_map_vertices_indexToLabel ||
          ptr_map_vertices_indexToLabel->size() == 0) {}

ostream &operator<<(ostream &o, const Graph &g) {
  unsigned int iMax(g.vertices.size()), i;

  // -----------------------------------------------------------------------
  // Nom du graphe
  o << "\t\t" << (g.isSpherical ? "" : "T") << (char)(g.type + 65)
    << (g.isSpherical ? iMax : iMax - 1) << " ; ";

  // -----------------------------------------------------------------------
  // sommets qui constituent le graphe
  for (i = 0; i < iMax; i++) {
    if (g.b_map_vertices_indexToLabelIsEmpty)
      o << (g.vertices[i] + 1) << " ";
    else
      o << ((*g.ptr_map_vertices_indexToLabel)[g.vertices[i]]) << " ";

    if (((g.type == 3 || g.type == 4) && i == (iMax - 2)) ||
        ((g.type == 1 && !g.isSpherical) &&
         (i == (iMax - 2) || i == 0))) // (Dn ou En) ou (\tilde Bn)
      o << "| ";
  }

  // -----------------------------------------------------------------------
  // Poids (pour G_2^k)
  if (g.dataSupp && g.isSpherical)
    o << " (" << g.dataSupp << ")";

  // Extended debbuging info
  /*
  o << "N={";
  for (i = 0; i < g.bVerticesLinkable.size(); i++)
  {
          if (!g.bVerticesLinkable[i])
                  o << (i + 1) << ",";
  }
  o << "}";*/

  o << endl;

  return o;
}

bool Graph::bIsSubgraphOf(const Graph *grBig) const {
  if (*this == *grBig)
    return true;

  if (isSpherical && grBig->isSpherical)
    return isSubgraphOf_spherical_spherical(grBig);
  else if (isSpherical && !grBig->isSpherical)
    return isSubgraphOf_spherical_euclidean(grBig);
  else // We do not test other cases
    throw(0);
}

bool Graph::isSubgraphOf_spherical_euclidean(const Graph *grBig) const {
  size_t verticesCount(vertices.size()),
      verticesBigCount(grBig->vertices.size());

  if (verticesCount == 1)
    return (find(grBig->vertices.begin(), grBig->vertices.end(), vertices[0]) !=
            grBig->vertices.end());

  // ---------------------------------------------------------------------
  // A_n < TA_m
  if (type == 0 && grBig->type == 0 && verticesBigCount > 2) {
    vector<short unsigned int>::const_iterator itBig(
        find(grBig->vertices.begin(), grBig->vertices.end(), vertices[0]));
    if (itBig == grBig->vertices.end())
      return false;

    int iDir;

    if ((itBig + 1) != grBig->vertices.end() && *(itBig + 1) == vertices[1])
      iDir = 1;
    else if ((itBig + 1) == grBig->vertices.end() &&
             grBig->vertices[0] == vertices[1])
      iDir = 1;
    else if ((itBig - 1) >= grBig->vertices.begin() &&
             *(itBig - 1) == vertices[1])
      iDir = -1;
    else if ((itBig - 1) < grBig->vertices.begin() &&
             grBig->vertices[verticesBigCount - 1] == vertices[1])
      iDir = -1;
    else
      return false;

    for (vector<short unsigned int>::const_iterator itSub(vertices.begin());
         itSub != vertices.end(); ++itSub) {
      if (*itSub != *itBig)
        return false;

      itBig += iDir;
      if (itBig == grBig->vertices.end())
        itBig = grBig->vertices.begin(); // move to the beginning
      else if (itBig < grBig->vertices.begin())
        itBig = grBig->vertices.end() - 1; // move to the end
    }

    return true;
  }

  // ---------------------------------------------------------------------
  // A_n < TB_m
  if (type == 0 && grBig->type == 1) {
    if (verticesBigCount == 4) {
      /* In this case, if we remove the first edge, the obtained A3 is not
       * correct */
      vector<short unsigned int> iVTemp;
      iVTemp.push_back(grBig->vertices[2]);
      iVTemp.push_back(grBig->vertices[1]);
      iVTemp.push_back(grBig->vertices[3]);

      return isAnSubAm(vertices, iVTemp);
    } else {
      Graph g(vector<short unsigned int>(grBig->vertices.begin() + 1,
                                         grBig->vertices.end()),
              0, vector<bool>(0), 3, true, 0);

      return isSubgraphOf_spherical_spherical(&g);
    }
  }

  // ---------------------------------------------------------------------
  // A_n < TC_m
  if (type == 0 && grBig->type == 2) {
    return isAnSubAm(vertices,
                     vector<short unsigned int>(grBig->vertices.begin() + 1,
                                                grBig->vertices.end() - 1));
  }

  // ---------------------------------------------------------------------
  // A_n < TD_m
  if (type == 0 && grBig->type == 3) {
    if (verticesCount <= 3) // could be in one of the ends
    {
      vector<short unsigned int> iVTemp;

      // left hand side
      iVTemp.push_back(grBig->vertices[0]);
      iVTemp.push_back(grBig->vertices[2]);
      iVTemp.push_back(grBig->vertices[1]);
      if (isAnSubAm(vertices, iVTemp))
        return true;
      iVTemp.clear();

      // right hand side
      iVTemp.push_back(grBig->vertices[verticesBigCount - 2]);
      iVTemp.push_back(grBig->vertices[verticesBigCount - 3]);
      iVTemp.push_back(grBig->vertices[verticesBigCount - 1]);
      if (isAnSubAm(vertices, iVTemp))
        return true;
    }

    vector<short unsigned int>::const_iterator itStart(
        find(grBig->vertices.begin(), grBig->vertices.end(), vertices[0]));
    vector<short unsigned int>::const_iterator itEnd(
        find(grBig->vertices.begin(), grBig->vertices.end(),
             vertices[verticesCount - 1]));

    vector<short unsigned int> iVTemp(grBig->vertices.begin() + 2,
                                      grBig->vertices.end() -
                                          2); // basis (middle)

    if (itStart == grBig->vertices.begin())
      iVTemp.insert(iVTemp.begin(), grBig->vertices[0]);
    else if (itStart == (grBig->vertices.begin() + 1))
      iVTemp.insert(iVTemp.begin(), grBig->vertices[1]);

    if (itEnd == grBig->vertices.end())
      iVTemp.push_back(grBig->vertices[verticesBigCount - 1]);
    else if (itEnd == (grBig->vertices.end() - 1))
      iVTemp.push_back(grBig->vertices[verticesBigCount - 2]);

    return isAnSubAm(vertices, iVTemp);
  }

  // ---------------------------------------------------------------------
  // A_n < TE_6
  if (type == 0 && grBig->type == 4 && verticesBigCount == 7) {
    vector<short unsigned int> iVTemp(grBig->vertices.begin(),
                                      grBig->vertices.begin() + 2);
    iVTemp.push_back(grBig->vertices[6]);
    iVTemp.push_back(grBig->vertices[4]);
    iVTemp.push_back(grBig->vertices[5]);
    if (isAnSubAm(vertices, iVTemp))
      return true;

    iVTemp = vector<short unsigned int>(grBig->vertices.begin(),
                                        grBig->vertices.begin() + 2);
    iVTemp.push_back(grBig->vertices[6]);
    iVTemp.push_back(grBig->vertices[2]);
    iVTemp.push_back(grBig->vertices[3]);
    if (isAnSubAm(vertices, iVTemp))
      return true;

    iVTemp.clear();
    iVTemp.push_back(grBig->vertices[5]);
    iVTemp.push_back(grBig->vertices[4]);
    iVTemp.push_back(grBig->vertices[6]);
    iVTemp.push_back(grBig->vertices[2]);
    iVTemp.push_back(grBig->vertices[3]);
    if (isAnSubAm(vertices, iVTemp))
      return true;
  }

  // ---------------------------------------------------------------------
  // A_n < TE_m
  if (type == 0 && grBig->type == 4) {
    vector<short unsigned int> iVTemp;

    //-------------------------------------
    // basis
    iVTemp = vector<short unsigned int>(grBig->vertices.begin(),
                                        grBig->vertices.end() - 1);

    if (isAnSubAm(vertices, iVTemp))
      return true;

    iVTemp.clear();

    // -------------------------------------
    // left & queue or right & queue
    if (verticesBigCount == 8) // TE7
    {
      iVTemp = vector<short unsigned int>(grBig->vertices.begin(),
                                          grBig->vertices.begin() + 4);
      iVTemp.push_back(grBig->vertices[7]);
      if (isAnSubAm(vertices, iVTemp))
        return true;

      iVTemp.clear();
      iVTemp = vector<short unsigned int>(grBig->vertices.begin() + 3,
                                          grBig->vertices.end() - 1);
      iVTemp.insert(iVTemp.begin(), grBig->vertices[7]);

      return isAnSubAm(vertices, iVTemp);
    } else if (verticesBigCount == 9) // TE8
    {
      iVTemp = vector<short unsigned int>(grBig->vertices.begin(),
                                          grBig->vertices.begin() + 3);
      iVTemp.push_back(grBig->vertices[8]);
      if (isAnSubAm(vertices, iVTemp))
        return true;

      iVTemp.clear();
      iVTemp = vector<short unsigned int>(grBig->vertices.begin() + 2,
                                          grBig->vertices.end() - 1);
      iVTemp.insert(iVTemp.begin(), grBig->vertices[8]);
      return isAnSubAm(vertices, iVTemp);
    }
  }

  // ---------------------------------------------------------------------
  // A_n < TF_4
  if (type == 0 && grBig->type == 5) {
    if (verticesCount > 3)
      return false;

    if (isAnSubAm(vertices,
                  vector<short unsigned int>(grBig->vertices.begin(),
                                             grBig->vertices.begin() + 3)))
      return true;

    if (verticesCount == 2)
      return isAnSubAm(vertices,
                       vector<short unsigned int>(grBig->vertices.end() - 2,
                                                  grBig->vertices.end()));

    return false;
  }

  // ---------------------------------------------------------------------
  // A_n < TG_2
  if (type == 0 && grBig->type == 6) {
    if (verticesCount > 2)
      return false;

    return isAnSubAm(vertices,
                     vector<short unsigned int>(grBig->vertices.begin(),
                                                grBig->vertices.begin() + 2));
  }

  // ---------------------------------------------------------------------
  // B_n < TB_m
  if (type == 1 && grBig->type == 1) {
    vector<short unsigned int> iVTemp(vertices);
    reverse(iVTemp.begin(), iVTemp.end());

    vector<short unsigned int>::const_iterator itSearch(
        find(grBig->vertices.begin(), grBig->vertices.end(), vertices[0]));

    if (itSearch == grBig->vertices.end())
      return false;

    // if we don't use any of the two ends (right ends)
    if (itSearch + 1 != grBig->vertices.end() &&
        itSearch + 2 != grBig->vertices.end())
      return (iVTemp == vector<short unsigned int>(grBig->vertices.begin(),
                                                   itSearch + 1));
    else if (itSearch + 1 == grBig->vertices.end()) {

      vector<short unsigned int> iVTemp2(
          vector<short unsigned int>(grBig->vertices.begin(), itSearch - 1));
      iVTemp2.push_back(grBig->vertices[verticesBigCount - 1]);

      return (iVTemp == iVTemp2);

    } else // itSearch + 2 == grBig->vertices.end()
    {
      vector<short unsigned int> iVTemp2(
          vector<short unsigned int>(grBig->vertices.begin(), itSearch));
      iVTemp2.push_back(grBig->vertices[verticesBigCount - 2]);

      return (iVTemp == iVTemp2);
    }

    return false;
  }

  // ---------------------------------------------------------------------
  // B_n < TC_m
  if (type == 1 && grBig->type == 2) {
    vector<short unsigned int>::const_iterator itSearch(
        find(grBig->vertices.begin(), grBig->vertices.end(), vertices[0]));

    if (itSearch == grBig->vertices.end())
      return false;

    // [3, ..., 3, 4] < [3, ..., 3, 4]
    if (vertices[verticesCount - 1] == grBig->vertices[verticesBigCount - 1])
      return (vertices ==
              vector<short unsigned int>(itSearch, grBig->vertices.end()));

    // [3, ..., 3, 4] < [4, 3, ..., 3]
    if (vertices[verticesCount - 1] == grBig->vertices[0]) {
      vector<short unsigned int> iVTemp(grBig->vertices.begin(), itSearch + 1);
      reverse(iVTemp.begin(), iVTemp.end());

      return (iVTemp == vertices);
    }
  }

  // ---------------------------------------------------------------------
  // B_n < TF_4
  if (type == 1 && grBig->type == 5) {
    if (verticesCount > 4)
      return false;

    // ----------------------------------------
    // ([4,3] OR [4,3,3]) < ([3,3,4])
    vector<short unsigned int> iVTemp(grBig->vertices.begin() +
                                          (verticesCount == 4 ? 0 : 1),
                                      grBig->vertices.end() - 1);
    if (vertices == iVTemp)
      return true;

    // ----------------------------------------
    // [4,3] < ([4,3])
    iVTemp = vector<short unsigned int>(grBig->vertices.begin() + 2,
                                        grBig->vertices.end());
    reverse(iVTemp.begin(), iVTemp.end());
    return (vertices == iVTemp);
  }

  // ---------------------------------------------------------------------
  // D_n < TB_m
  if (type == 3 && grBig->type == 1) {
    if (verticesBigCount < 5)
      return false;

    // "central" node
    if (vertices[verticesCount - 3] != grBig->vertices[verticesBigCount - 3])
      return false;

    Graph g(vector<short unsigned int>(grBig->vertices.begin() + 1,
                                       grBig->vertices.end()),
            0, vector<bool>(false), 3, true, 0);

    return isSubgraphOf_spherical_spherical(&g);
  }

  // ---------------------------------------------------------------------
  // D_4 < TD_4
  if (type == 3 && verticesCount == 4 && grBig->type == 3 &&
      verticesBigCount == 5) {
    if (vertices[1] != grBig->vertices[4]) // comparison of the central node
      return false;

    /*
     * Remark:
     * We should use find(grBig->vertices.begin(), grBig->vertices.end() - 1,
     * -), instead of find(grBig->vertices.begin(), grBig->vertices.end(), -)
     * but it doesn't really matter since the elements of grBig->vertices are
     * all distincts.
     */
    if (find(grBig->vertices.begin(), grBig->vertices.end(), vertices[0]) ==
        grBig->vertices.end())
      return false;
    if (find(grBig->vertices.begin(), grBig->vertices.end(), vertices[2]) ==
        grBig->vertices.end())
      return false;
    if (find(grBig->vertices.begin(), grBig->vertices.end(), vertices[3]) ==
        grBig->vertices.end())
      return false;

    return true;
  }

  // ---------------------------------------------------------------------
  // D_n < TD_m
  if (type == 3 && grBig->type == 3) {
    // Test 1
    Graph g(vector<short unsigned int>(grBig->vertices.begin() + 1,
                                       grBig->vertices.end()),
            0, vector<bool>(false), 3, true, 0);
    if (isSubgraphOf_spherical_spherical(&g))
      return true;

    // Test 2
    vector<short unsigned int> iVTemp(vector<short unsigned int>(
        grBig->vertices.begin() + 2, grBig->vertices.end()));
    iVTemp.insert(iVTemp.begin(), grBig->vertices[0]);
    g = Graph(iVTemp, 0, vector<bool>(false), 3, true, 0);
    if (isSubgraphOf_spherical_spherical(&g))
      return true;

    // Test 3
    iVTemp = vector<short unsigned int>(vector<short unsigned int>(
        grBig->vertices.begin(), grBig->vertices.end() - 1));
    iVTemp[0] = max(grBig->vertices[0], grBig->vertices[1]);
    iVTemp[1] = min(grBig->vertices[0], grBig->vertices[1]);
    reverse(iVTemp.begin(), iVTemp.end());
    g = Graph(iVTemp, 0, vector<bool>(false), 3, true, 0);
    if (isSubgraphOf_spherical_spherical(&g))
      return true;

    // Test 4
    iVTemp = vector<short unsigned int>(vector<short unsigned int>(
        grBig->vertices.begin(), grBig->vertices.end() - 2));
    iVTemp[0] = max(grBig->vertices[0], grBig->vertices[1]);
    iVTemp[1] = min(grBig->vertices[0], grBig->vertices[1]);
    iVTemp.push_back(grBig->vertices[verticesBigCount - 1]);
    reverse(iVTemp.begin(), iVTemp.end());
    g = Graph(iVTemp, 0, vector<bool>(false), 3, true, 0);
    if (isSubgraphOf_spherical_spherical(&g))
      return true;
  }

  // ---------------------------------------------------------------------
  // D_n < TE_6
  if (type == 3 && grBig->type == 4 && verticesBigCount == 7) {
    // -------------------------------------------------
    // first test
    vector<short unsigned int> iVTemp(grBig->vertices.begin(),
                                      grBig->vertices.begin() + 2);
    iVTemp.push_back(grBig->vertices[6]);
    iVTemp.push_back(min(grBig->vertices[2], grBig->vertices[4]));
    iVTemp.push_back(max(grBig->vertices[2], grBig->vertices[4]));
    Graph g(iVTemp, 0, vector<bool>(false), 3, true, 0);

    if (isSubgraphOf_spherical_spherical(&g))
      return true;

    // -------------------------------------------------
    // second test
    iVTemp.clear();
    iVTemp.push_back(grBig->vertices[5]);
    iVTemp.push_back(grBig->vertices[4]);
    iVTemp.push_back(grBig->vertices[6]);
    iVTemp.push_back(min(grBig->vertices[2], grBig->vertices[1]));
    iVTemp.push_back(max(grBig->vertices[2], grBig->vertices[1]));

    g = Graph(iVTemp, 0, vector<bool>(false), 3, true, 0);

    if (isSubgraphOf_spherical_spherical(&g))
      return true;

    // -------------------------------------------------
    // third test
    iVTemp.clear();
    iVTemp.push_back(grBig->vertices[3]);
    iVTemp.push_back(grBig->vertices[2]);
    iVTemp.push_back(grBig->vertices[6]);
    iVTemp.push_back(min(grBig->vertices[4], grBig->vertices[1]));
    iVTemp.push_back(max(grBig->vertices[4], grBig->vertices[1]));

    g = Graph(iVTemp, 0, vector<bool>(false), 3, true, 0);

    if (isSubgraphOf_spherical_spherical(&g))
      return true;

    return false;
  }

  // ---------------------------------------------------------------------
  // D_n < TE_m
  if (type == 3 && grBig->type == 4) {
    unsigned int iBigQueueIndex(verticesBigCount == 8 ? 3 : 2);

    // -------------------------------------------------
    // first test
    vector<short unsigned int> iVTemp(
        grBig->vertices.begin(), grBig->vertices.begin() + iBigQueueIndex + 1);
    iVTemp.push_back(min(grBig->vertices[iBigQueueIndex + 1],
                         grBig->vertices[verticesBigCount - 1]));
    iVTemp.push_back(max(grBig->vertices[iBigQueueIndex + 1],
                         grBig->vertices[verticesBigCount - 1]));
    Graph g(iVTemp, 0, vector<bool>(false), 3, true, 0);

    if (isSubgraphOf_spherical_spherical(&g))
      return true;

    // -------------------------------------------------
    // second test
    iVTemp = vector<short unsigned int>(
        grBig->vertices.begin() + iBigQueueIndex, grBig->vertices.end() - 1);
    reverse(iVTemp.begin(), iVTemp.end());
    iVTemp.push_back(min(grBig->vertices[iBigQueueIndex - 1],
                         grBig->vertices[verticesBigCount - 1]));
    iVTemp.push_back(max(grBig->vertices[iBigQueueIndex - 1],
                         grBig->vertices[verticesBigCount - 1]));
    g = Graph(iVTemp, 0, vector<bool>(false), 3, true, 0);

    if (isSubgraphOf_spherical_spherical(&g))
      return true;

    return false;
  }

  // ---------------------------------------------------------------------
  // E_6 < TE_6
  if (type == 4 && verticesCount == 6 && grBig->type == 4 &&
      verticesBigCount == 7) {
    vector<short unsigned int> iVBigBasis;

    // We must find out what the base is
    if (vertices[5] == grBig->vertices[4]) {
      iVBigBasis.push_back(grBig->vertices[0]);
      iVBigBasis.push_back(grBig->vertices[1]);
      iVBigBasis.push_back(grBig->vertices[6]);
      iVBigBasis.push_back(grBig->vertices[2]);
      iVBigBasis.push_back(grBig->vertices[3]);
    } else if (vertices[5] == grBig->vertices[1]) {
      iVBigBasis.push_back(grBig->vertices[5]);
      iVBigBasis.push_back(grBig->vertices[4]);
      iVBigBasis.push_back(grBig->vertices[6]);
      iVBigBasis.push_back(grBig->vertices[2]);
      iVBigBasis.push_back(grBig->vertices[3]);
    } else if (vertices[5] == grBig->vertices[2]) {
      iVBigBasis.push_back(grBig->vertices[0]);
      iVBigBasis.push_back(grBig->vertices[1]);
      iVBigBasis.push_back(grBig->vertices[6]);
      iVBigBasis.push_back(grBig->vertices[4]);
      iVBigBasis.push_back(grBig->vertices[5]);
    } else
      return false;

    if (iVBigBasis[0] > iVBigBasis[4])
      reverse(iVBigBasis.begin(), iVBigBasis.end());

    return (iVBigBasis ==
            vector<short unsigned int>(vertices.begin(), vertices.end() - 1));
  }

  // ---------------------------------------------------------------------
  // E_n < TE_7, TE_8
  if (type == 4 && grBig->type == 4) {
    if (vertices[verticesCount - 1] != grBig->vertices[verticesBigCount - 1])
      return false;

    // E_6, E_7 < TE7
    if (verticesCount <= 7 && verticesBigCount == 8) {
      vector<short unsigned int> iVBigBasis(
          grBig->vertices.begin() + 1, grBig->vertices.begin() + verticesCount);
      if (isAnSubAm(
              vector<short unsigned int>(vertices.begin(), vertices.end() - 1),
              iVBigBasis))
        return true;

      if (verticesCount == 7 && verticesBigCount == 8) // E_7 < TE_7
      {
        iVBigBasis = vector<short unsigned int>(grBig->vertices.begin(),
                                                grBig->vertices.begin() + 6);
        reverse(iVBigBasis.begin(), iVBigBasis.end());

        return (vector<short unsigned int>(vertices.begin(),
                                           vertices.end() - 1) == iVBigBasis);
      }

      return false;
    } else if (verticesBigCount == 9) // < TE_8
    {
      vector<short unsigned int> iVBigBasis(
          grBig->vertices.begin(), grBig->vertices.begin() + verticesCount - 1);

      if (vector<short unsigned int>(vertices.begin(), vertices.end() - 1) ==
          iVBigBasis)
        return true;

      if (verticesCount == 6) // E_5 is symmetric
      {
        reverse(iVBigBasis.begin(), iVBigBasis.end());

        return ((vector<short unsigned int>(vertices.begin(),
                                            vertices.end() - 1) == iVBigBasis));
      }

      return false;
    } else
      return false;
  }

  // ---------------------------------------------------------------------
  // F_4 < TF_4
  if (type == 5 && grBig->type == 5) {
    return isAnSubAm(vertices,
                     vector<short unsigned int>(grBig->vertices.begin() + 1,
                                                grBig->vertices.end()));
  }

  // ---------------------------------------------------------------------
  // G_4 < TB_m
  if (type == 6 && dataSupp == 4 && grBig->type == 1) {
    return ((vertices[0] == grBig->vertices[0] &&
             vertices[1] == grBig->vertices[1]) ||
            (vertices[0] == grBig->vertices[1] &&
             vertices[1] == grBig->vertices[0]));
  }

  // ---------------------------------------------------------------------
  // G_4 < TC_m
  if (type == 6 && dataSupp == 4 && grBig->type == 2) {
    if ((vertices[0] == grBig->vertices[0] &&
         vertices[1] == grBig->vertices[1]) ||
        (vertices[0] == grBig->vertices[1] &&
         vertices[1] == grBig->vertices[0]))
      return true;

    return ((vertices[0] == grBig->vertices[verticesBigCount - 1] &&
             vertices[1] == grBig->vertices[verticesBigCount - 2]) ||
            (vertices[0] == grBig->vertices[verticesBigCount - 2] &&
             vertices[1] == grBig->vertices[verticesBigCount - 1]));
  }

  // ---------------------------------------------------------------------
  // G_4 < TF_4
  if (type == 6 && dataSupp == 4 && grBig->type == 5) {
    return ((vertices[0] == grBig->vertices[2] &&
             vertices[1] == grBig->vertices[3]) ||
            (vertices[0] == grBig->vertices[3] &&
             vertices[1] == grBig->vertices[2]));
  }

  // ---------------------------------------------------------------------
  // G_6 < TG_2
  if (type == 6 && dataSupp == 6 && grBig->type == 6) {
    return ((vertices[0] == grBig->vertices[1] &&
             vertices[1] == grBig->vertices[2]) ||
            (vertices[1] == grBig->vertices[1] &&
             vertices[0] == grBig->vertices[2]));
  }

  return false;
}

bool Graph::isSubgraphOf_spherical_spherical(const Graph *grBig) const {
  size_t verticesCount = vertices.size();
  size_t verticesBigCount = grBig->vertices.size();

  if (verticesBigCount < verticesCount)
    return false;

  if (verticesCount == 1)
    return find(grBig->vertices.begin(), grBig->vertices.end(), vertices[0]) !=
           grBig->vertices.end();

  // ---------------------------------------------------------------------
  // A_n < A_m
  if (type == 0 && grBig->type == 0)
    return isAnSubAm(vertices, grBig->vertices);

  // ---------------------------------------------------------------------
  // A_n < B_m
  if (type == 0 && grBig->type == 1) {
    return isAnSubAm(vertices,
                     vector<short unsigned int>(grBig->vertices.begin(),
                                                grBig->vertices.end() - 1));
  }

  // ---------------------------------------------------------------------
  // A_n < D_m
  if (type == 0 && grBig->type == 3) {
    // base plus une des deux extrémités
    vector<short unsigned int> iVTemp(grBig->vertices.begin(),
                                      grBig->vertices.end() - 1);
    if (isAnSubAm(vertices, iVTemp))
      return true;

    // base plus l'autre extrémité
    iVTemp = vector<short unsigned int>(grBig->vertices.begin(),
                                        grBig->vertices.end() - 2);
    iVTemp.push_back(grBig->vertices[verticesBigCount - 1]);
    if (isAnSubAm(vertices, iVTemp))
      return true;

    // A_3 dans le "bout du Y"
    if (verticesCount <= 3 && verticesBigCount > 3) {
      iVTemp.clear();
      iVTemp.push_back(min(grBig->vertices[verticesBigCount - 1],
                           grBig->vertices[verticesBigCount - 2]));
      iVTemp.push_back(grBig->vertices[verticesBigCount - 3]);
      iVTemp.push_back(max(grBig->vertices[verticesBigCount - 1],
                           grBig->vertices[verticesBigCount - 2]));

      return isAnSubAm(vertices, iVTemp);
    }

    return false;
  }

  // ---------------------------------------------------------------------
  // A_n < E_m
  if (type == 0 && grBig->type == 4) {
    // A_n dans la base du E_n
    vector<short unsigned int> iVTemp(grBig->vertices.begin(),
                                      grBig->vertices.end() - 1);
    if (isAnSubAm(vertices, iVTemp))
      return true;

    // A_n contient la queue du E_m (début)
    iVTemp = vector<short unsigned int>(grBig->vertices.begin(),
                                        grBig->vertices.begin() + 3);
    iVTemp.push_back(grBig->vertices[verticesBigCount - 1]);
    if (isAnSubAm(vertices, iVTemp))
      return true;

    // A_n contient la queue du E_m (fin)
    iVTemp = vector<short unsigned int>(grBig->vertices.begin() + 2,
                                        grBig->vertices.end() - 1);
    iVTemp.insert(iVTemp.begin(), grBig->vertices[verticesBigCount - 1]);

    return isAnSubAm(vertices, iVTemp);
  }

  // ---------------------------------------------------------------------
  // A_2 < F_4
  if (type == 0 && verticesCount == 2 && grBig->type == 5) {
    return ((vertices[0] == grBig->vertices[0] &&
             vertices[1] == grBig->vertices[1]) ||
            (vertices[1] == grBig->vertices[0] &&
             vertices[0] == grBig->vertices[1]) ||
            (vertices[0] == grBig->vertices[2] &&
             vertices[1] == grBig->vertices[3]) ||
            (vertices[1] == grBig->vertices[2] &&
             vertices[0] == grBig->vertices[3]));
  }

  // ---------------------------------------------------------------------
  // A_n < H_m
  if (type == 0 && grBig->type == 7) {
    return isAnSubAm(vertices,
                     vector<short unsigned int>(grBig->vertices.begin(),
                                                grBig->vertices.end() - 1));
  }

  // ---------------------------------------------------------------------
  // B_n < B_m
  if (type == 1 && grBig->type == 1) {
    return isAnSubAm(vertices, grBig->vertices);
  }

  // ---------------------------------------------------------------------
  // B_3 < F_4
  if (type == 1 && grBig->type == 5 && verticesCount == 3) {
    if (vertices == vector<short unsigned int>(grBig->vertices.begin(),
                                               grBig->vertices.begin() + 3))
      return true;

    vector<short unsigned int> iVTemp(grBig->vertices.begin() + 1,
                                      grBig->vertices.begin() + 4);
    reverse(iVTemp.begin(), iVTemp.end());

    return (vertices == iVTemp);
  }

  // ---------------------------------------------------------------------
  // D_n < D_m
  if (type == 3 && grBig->type == 3) {
    if (verticesCount == 4) // D_4 est très symétrique
    {
      // central node
      if (vertices[1] != grBig->vertices[verticesBigCount - 3])
        return false;

      vector<short unsigned int> iVTemp1;
      iVTemp1.push_back(vertices[0]);
      iVTemp1.push_back(vertices[2]);
      iVTemp1.push_back(vertices[3]);

      vector<short unsigned int> iVTemp2;
      iVTemp2.push_back(grBig->vertices[verticesBigCount - 1]);
      iVTemp2.push_back(grBig->vertices[verticesBigCount - 2]);
      iVTemp2.push_back(grBig->vertices[verticesBigCount - 4]);

      sort(iVTemp1.begin(), iVTemp1.end());
      sort(iVTemp2.begin(), iVTemp2.end());

      return (iVTemp1 == iVTemp2);
    }

    // Autres D_n
    if (vertices[verticesCount - 2] != grBig->vertices[verticesBigCount - 2] ||
        vertices[verticesCount - 1] != grBig->vertices[verticesBigCount - 1])
      return false;

    vector<short unsigned int> iVTemp1(vertices.begin(),
                                       vertices.begin() + verticesCount - 2);
    vector<short unsigned int> iVTemp2(
        grBig->vertices.begin() + verticesBigCount - verticesCount,
        grBig->vertices.begin() + verticesBigCount - 2);

    return (iVTemp1 == iVTemp2);
  }

  // ---------------------------------------------------------------------
  // D_n < E_m
  if (type == 3 && grBig->type == 4) {
    // -------------------------------------------------
    // first test
    vector<short unsigned int> iVTemp(grBig->vertices.begin(),
                                      grBig->vertices.begin() + 3);
    iVTemp.push_back(
        min(grBig->vertices[3], grBig->vertices[verticesBigCount - 1]));
    iVTemp.push_back(
        max(grBig->vertices[3], grBig->vertices[verticesBigCount - 1]));
    Graph g(iVTemp, 0, vector<bool>(false), 3, true, 0);

    if (isSubgraphOf_spherical_spherical(&g))
      return true;

    // -------------------------------------------------
    // second test
    iVTemp = vector<short unsigned int>(grBig->vertices.begin() + 2,
                                        grBig->vertices.end() - 1);
    reverse(iVTemp.begin(), iVTemp.end());
    iVTemp.push_back(
        min(grBig->vertices[1], grBig->vertices[verticesBigCount - 1]));
    iVTemp.push_back(
        max(grBig->vertices[1], grBig->vertices[verticesBigCount - 1]));
    g = Graph(iVTemp, 0, vector<bool>(false), 3, true, 0);

    if (isSubgraphOf_spherical_spherical(&g))
      return true;

    return false;
  }

  // ---------------------------------------------------------------------
  // E_n < E_m
  if (type == 4 && grBig->type == 4) {
    if (vertices[verticesCount - 1] != grBig->vertices[verticesBigCount - 1])
      return false;

    vector<short unsigned int> iVTemp1(vertices.begin(), vertices.end() - 1);
    vector<short unsigned int> iVTemp2(
        grBig->vertices.begin(), grBig->vertices.begin() + verticesCount - 1);

    if (iVTemp1 == iVTemp2)
      return true;

    if (verticesCount == 6) // E_5 a une base symétrique
    {
      reverse(iVTemp1.begin(), iVTemp1.end());
      if (iVTemp1 == iVTemp2)
        return true;
    }

    return false;
  }

  // ---------------------------------------------------------------------
  // G_4 < B_n
  if (type == 6 && dataSupp == 4 && grBig->type == 1) {
    return ((vertices[0] == grBig->vertices[verticesBigCount - 2] &&
             vertices[1] == grBig->vertices[verticesBigCount - 1]) ||
            (vertices[0] == grBig->vertices[verticesBigCount - 1] &&
             vertices[1] == grBig->vertices[verticesBigCount - 2]));
  }

  // ---------------------------------------------------------------------
  // G_m < G_m
  if (type == 6 && grBig->type == 6 && dataSupp == grBig->dataSupp &&
      vertices == grBig->vertices) {
    return ((vertices[0] == grBig->vertices[verticesBigCount - 2] &&
             vertices[1] == grBig->vertices[verticesBigCount - 1]) ||
            (vertices[0] == grBig->vertices[verticesBigCount - 1] &&
             vertices[1] == grBig->vertices[verticesBigCount - 2]));
  }

  // ---------------------------------------------------------------------
  // G_5 < H_m
  if ((type == 6 && dataSupp == 5) && grBig->type == 7) {
    return ((vertices[0] == grBig->vertices[verticesBigCount - 2] &&
             vertices[1] == grBig->vertices[verticesBigCount - 1]) ||
            (vertices[1] == grBig->vertices[verticesBigCount - 2] &&
             vertices[0] == grBig->vertices[verticesBigCount - 1]));
  }

  // ---------------------------------------------------------------------
  // H_3 < H_4
  if (type == 7 && grBig->type == 7) {
    return (vector<short unsigned int>(grBig->vertices.begin() + 1,
                                       grBig->vertices.begin() + verticesCount +
                                           1) == vertices);
  }

  return false;
}

bool Graph::isAnSubAm(const vector<short unsigned int> &iSubV,
                      const vector<short unsigned int> &iBigV) {
  vector<short unsigned int>::const_iterator itBig, itSub,
      it(find(iBigV.begin(), iBigV.end(), iSubV[0]));

  if (it == iBigV.end())
    return false;

  if (iSubV.size() == 1)
    return true;

  if (iSubV.size() > iBigV.size())
    return false;

  if ((it + 1) == iBigV.end() && *(it - 1) != iSubV[1])
    return false;

  // ----------------------------------------------------------------------
  // recherche en avant depuis it
  if ((it + 1) != iBigV.end() && *(it + 1) == iSubV[1]) {
    itBig = it;
    for (itSub = iSubV.begin(); itSub != iSubV.end() && itBig != iBigV.end();
         ++itSub) {
      if (*itSub != *itBig)
        break;
      ++itBig;
    }

    return (itSub == iSubV.end());
  }

  // ----------------------------------------------------------------------
  // recherche en arrière depuis it
  itBig = it;
  for (itSub = iSubV.begin(); itSub != iSubV.end() && itBig >= iBigV.begin();
       ++itSub) {
    if (*itSub != *itBig)
      return false;

    --itBig;
  }

  return (itSub == iSubV.end());
}

bool operator==(const Graph &g1, const Graph &g2) {
  return (g1.type == g2.type && g1.dataSupp == g2.dataSupp &&
          g1.vertices == g2.vertices && g1.isSpherical == g2.isSpherical);
}

bool operator<(const Graph &g1, const Graph &g2) {
  if (g1 == g2)
    return false;

  if (g1.isSpherical && !g2.isSpherical)
    return true;
  if (!g1.isSpherical && g2.isSpherical)
    return false;

  if (g1.vertices.size() < g2.vertices.size())
    return true;
  if (g1.vertices.size() > g2.vertices.size())
    return false;

  if (g1.type < g2.type)
    return true;
  if (g1.type > g2.type)
    return false;

  if (g1.type == 6 && g2.type == 6) {
    if (g1.dataSupp < g2.dataSupp)
      return true;
    if (g1.dataSupp > g2.dataSupp)
      return false;
  }

  if (g1.vertices < g2.vertices)
    return true;
  if (g1.vertices > g2.vertices)
    return false;

  throw(string("Graph::operator<: One missed case"));

  return false;
}
