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
             const vector<bool> &linkableVertices, const unsigned int &type,
             const bool &isSpherical, const unsigned int &dataSupp)
    : type(type), vertices(vertices), linkableVertices(linkableVertices),
      dataSupp(dataSupp), isSpherical(isSpherical),
      ptr_map_vertices_indexToLabel(ptr_map_vertices_indexToLabel),
      b_map_vertices_indexToLabelIsEmpty(
          !ptr_map_vertices_indexToLabel ||
          ptr_map_vertices_indexToLabel->size() == 0) {}

ostream &operator<<(ostream &o, const Graph &g) {
  unsigned int verticesCount(g.vertices.size()), i;

  // -----------------------------------------------------------------------
  // Nom du graphe
  o << "\t\t" << (g.isSpherical ? "" : "T") << (char)(g.type + 65)
    << (g.isSpherical ? verticesCount : verticesCount - 1) << " ; ";

  // -----------------------------------------------------------------------
  // sommets qui constituent le graphe
  for (i = 0; i < verticesCount; i++) {
    if (g.b_map_vertices_indexToLabelIsEmpty)
      o << (g.vertices[i] + 1) << " ";
    else
      o << ((*g.ptr_map_vertices_indexToLabel)[g.vertices[i]]) << " ";

    if (((g.type == 3 || g.type == 4) && i == (verticesCount - 2)) ||
        ((g.type == 1 && !g.isSpherical) &&
         (i == (verticesCount - 2) || i == 0))) // (Dn ou En) ou (\tilde Bn)
      o << "| ";
  }

  // -----------------------------------------------------------------------
  // Poids (pour G_2^k)
  if (g.dataSupp && g.isSpherical)
    o << " (" << g.dataSupp << ")";

  o << endl;

  return o;
}

bool Graph::isSubgraphOf(const Graph *grBig) const {
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
    auto itBig(
        find(grBig->vertices.begin(), grBig->vertices.end(), vertices[0]));
    if (itBig == grBig->vertices.end())
      return false;

    int direction;

    if ((itBig + 1) != grBig->vertices.end() && *(itBig + 1) == vertices[1])
      direction = 1;
    else if ((itBig + 1) == grBig->vertices.end() &&
             grBig->vertices[0] == vertices[1])
      direction = 1;
    else if ((itBig - 1) >= grBig->vertices.begin() &&
             *(itBig - 1) == vertices[1])
      direction = -1;
    else if ((itBig - 1) < grBig->vertices.begin() &&
             grBig->vertices[verticesBigCount - 1] == vertices[1])
      direction = -1;
    else
      return false;

    for (const auto &vertex : vertices) {
      if (vertex != *itBig)
        return false;

      itBig += direction;
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
      vector<short unsigned int> temp;
      temp.push_back(grBig->vertices[2]);
      temp.push_back(grBig->vertices[1]);
      temp.push_back(grBig->vertices[3]);

      return isAnSubAm(vertices, temp);
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
      vector<short unsigned int> temp;

      // left hand side
      temp.push_back(grBig->vertices[0]);
      temp.push_back(grBig->vertices[2]);
      temp.push_back(grBig->vertices[1]);
      if (isAnSubAm(vertices, temp))
        return true;
      temp.clear();

      // right hand side
      temp.push_back(grBig->vertices[verticesBigCount - 2]);
      temp.push_back(grBig->vertices[verticesBigCount - 3]);
      temp.push_back(grBig->vertices[verticesBigCount - 1]);
      if (isAnSubAm(vertices, temp))
        return true;
    }

    auto itStart =
        find(grBig->vertices.begin(), grBig->vertices.end(), vertices[0]);
    auto itEnd = find(grBig->vertices.begin(), grBig->vertices.end(),
                      vertices[verticesCount - 1]);

    vector<short unsigned int> temp(grBig->vertices.begin() + 2,
                                    grBig->vertices.end() -
                                        2); // basis (middle)

    if (itStart == grBig->vertices.begin())
      temp.insert(temp.begin(), grBig->vertices[0]);
    else if (itStart == (grBig->vertices.begin() + 1))
      temp.insert(temp.begin(), grBig->vertices[1]);

    if (itEnd == grBig->vertices.end())
      temp.push_back(grBig->vertices[verticesBigCount - 1]);
    else if (itEnd == (grBig->vertices.end() - 1))
      temp.push_back(grBig->vertices[verticesBigCount - 2]);

    return isAnSubAm(vertices, temp);
  }

  // ---------------------------------------------------------------------
  // A_n < TE_6
  if (type == 0 && grBig->type == 4 && verticesBigCount == 7) {
    vector<short unsigned int> temp(grBig->vertices.begin(),
                                    grBig->vertices.begin() + 2);
    temp.push_back(grBig->vertices[6]);
    temp.push_back(grBig->vertices[4]);
    temp.push_back(grBig->vertices[5]);
    if (isAnSubAm(vertices, temp))
      return true;

    temp = vector<short unsigned int>(grBig->vertices.begin(),
                                      grBig->vertices.begin() + 2);
    temp.push_back(grBig->vertices[6]);
    temp.push_back(grBig->vertices[2]);
    temp.push_back(grBig->vertices[3]);
    if (isAnSubAm(vertices, temp))
      return true;

    temp.clear();
    temp.push_back(grBig->vertices[5]);
    temp.push_back(grBig->vertices[4]);
    temp.push_back(grBig->vertices[6]);
    temp.push_back(grBig->vertices[2]);
    temp.push_back(grBig->vertices[3]);
    if (isAnSubAm(vertices, temp))
      return true;
  }

  // ---------------------------------------------------------------------
  // A_n < TE_m
  if (type == 0 && grBig->type == 4) {
    vector<short unsigned int> temp;

    //-------------------------------------
    // basis
    temp = vector<short unsigned int>(grBig->vertices.begin(),
                                      grBig->vertices.end() - 1);

    if (isAnSubAm(vertices, temp))
      return true;

    temp.clear();

    // -------------------------------------
    // left & queue or right & queue
    if (verticesBigCount == 8) // TE7
    {
      temp = vector<short unsigned int>(grBig->vertices.begin(),
                                        grBig->vertices.begin() + 4);
      temp.push_back(grBig->vertices[7]);
      if (isAnSubAm(vertices, temp))
        return true;

      temp.clear();
      temp = vector<short unsigned int>(grBig->vertices.begin() + 3,
                                        grBig->vertices.end() - 1);
      temp.insert(temp.begin(), grBig->vertices[7]);

      return isAnSubAm(vertices, temp);
    } else if (verticesBigCount == 9) // TE8
    {
      temp = vector<short unsigned int>(grBig->vertices.begin(),
                                        grBig->vertices.begin() + 3);
      temp.push_back(grBig->vertices[8]);
      if (isAnSubAm(vertices, temp))
        return true;

      temp.clear();
      temp = vector<short unsigned int>(grBig->vertices.begin() + 2,
                                        grBig->vertices.end() - 1);
      temp.insert(temp.begin(), grBig->vertices[8]);
      return isAnSubAm(vertices, temp);
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
    vector<short unsigned int> temp(vertices);
    reverse(temp.begin(), temp.end());

    auto itSearch =
        find(grBig->vertices.begin(), grBig->vertices.end(), vertices[0]);

    if (itSearch == grBig->vertices.end())
      return false;

    // if we don't use any of the two ends (right ends)
    if (itSearch + 1 != grBig->vertices.end() &&
        itSearch + 2 != grBig->vertices.end())
      return (temp == vector<short unsigned int>(grBig->vertices.begin(),
                                                 itSearch + 1));
    else if (itSearch + 1 == grBig->vertices.end()) {

      auto temp2(
          vector<short unsigned int>(grBig->vertices.begin(), itSearch - 1));
      temp2.push_back(grBig->vertices[verticesBigCount - 1]);

      return (temp == temp2);

    } else // itSearch + 2 == grBig->vertices.end()
    {
      vector<short unsigned int> temp2(
          vector<short unsigned int>(grBig->vertices.begin(), itSearch));
      temp2.push_back(grBig->vertices[verticesBigCount - 2]);

      return (temp == temp2);
    }

    return false;
  }

  // ---------------------------------------------------------------------
  // B_n < TC_m
  if (type == 1 && grBig->type == 2) {
    auto itSearch =
        find(grBig->vertices.begin(), grBig->vertices.end(), vertices[0]);

    if (itSearch == grBig->vertices.end())
      return false;

    // [3, ..., 3, 4] < [3, ..., 3, 4]
    if (vertices[verticesCount - 1] == grBig->vertices[verticesBigCount - 1])
      return (vertices ==
              vector<short unsigned int>(itSearch, grBig->vertices.end()));

    // [3, ..., 3, 4] < [4, 3, ..., 3]
    if (vertices[verticesCount - 1] == grBig->vertices[0]) {
      vector<short unsigned int> temp(grBig->vertices.begin(), itSearch + 1);
      reverse(temp.begin(), temp.end());

      return (temp == vertices);
    }
  }

  // ---------------------------------------------------------------------
  // B_n < TF_4
  if (type == 1 && grBig->type == 5) {
    if (verticesCount > 4)
      return false;

    // ----------------------------------------
    // ([4,3] OR [4,3,3]) < ([3,3,4])
    vector<short unsigned int> temp(grBig->vertices.begin() +
                                        (verticesCount == 4 ? 0 : 1),
                                    grBig->vertices.end() - 1);
    if (vertices == temp)
      return true;

    // ----------------------------------------
    // [4,3] < ([4,3])
    temp = vector<short unsigned int>(grBig->vertices.begin() + 2,
                                      grBig->vertices.end());
    reverse(temp.begin(), temp.end());
    return (vertices == temp);
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
    vector<short unsigned int> temp(vector<short unsigned int>(
        grBig->vertices.begin() + 2, grBig->vertices.end()));
    temp.insert(temp.begin(), grBig->vertices[0]);
    g = Graph(temp, 0, vector<bool>(false), 3, true, 0);
    if (isSubgraphOf_spherical_spherical(&g))
      return true;

    // Test 3
    temp = vector<short unsigned int>(vector<short unsigned int>(
        grBig->vertices.begin(), grBig->vertices.end() - 1));
    temp[0] = max(grBig->vertices[0], grBig->vertices[1]);
    temp[1] = min(grBig->vertices[0], grBig->vertices[1]);
    reverse(temp.begin(), temp.end());
    g = Graph(temp, 0, vector<bool>(false), 3, true, 0);
    if (isSubgraphOf_spherical_spherical(&g))
      return true;

    // Test 4
    temp = vector<short unsigned int>(vector<short unsigned int>(
        grBig->vertices.begin(), grBig->vertices.end() - 2));
    temp[0] = max(grBig->vertices[0], grBig->vertices[1]);
    temp[1] = min(grBig->vertices[0], grBig->vertices[1]);
    temp.push_back(grBig->vertices[verticesBigCount - 1]);
    reverse(temp.begin(), temp.end());
    g = Graph(temp, 0, vector<bool>(false), 3, true, 0);
    if (isSubgraphOf_spherical_spherical(&g))
      return true;
  }

  // ---------------------------------------------------------------------
  // D_n < TE_6
  if (type == 3 && grBig->type == 4 && verticesBigCount == 7) {
    // -------------------------------------------------
    // first test
    vector<short unsigned int> temp(grBig->vertices.begin(),
                                    grBig->vertices.begin() + 2);
    temp.push_back(grBig->vertices[6]);
    temp.push_back(min(grBig->vertices[2], grBig->vertices[4]));
    temp.push_back(max(grBig->vertices[2], grBig->vertices[4]));
    Graph g(temp, 0, vector<bool>(false), 3, true, 0);

    if (isSubgraphOf_spherical_spherical(&g))
      return true;

    // -------------------------------------------------
    // second test
    temp.clear();
    temp.push_back(grBig->vertices[5]);
    temp.push_back(grBig->vertices[4]);
    temp.push_back(grBig->vertices[6]);
    temp.push_back(min(grBig->vertices[2], grBig->vertices[1]));
    temp.push_back(max(grBig->vertices[2], grBig->vertices[1]));

    g = Graph(temp, 0, vector<bool>(false), 3, true, 0);

    if (isSubgraphOf_spherical_spherical(&g))
      return true;

    // -------------------------------------------------
    // third test
    temp.clear();
    temp.push_back(grBig->vertices[3]);
    temp.push_back(grBig->vertices[2]);
    temp.push_back(grBig->vertices[6]);
    temp.push_back(min(grBig->vertices[4], grBig->vertices[1]));
    temp.push_back(max(grBig->vertices[4], grBig->vertices[1]));

    g = Graph(temp, 0, vector<bool>(false), 3, true, 0);

    if (isSubgraphOf_spherical_spherical(&g))
      return true;

    return false;
  }

  // ---------------------------------------------------------------------
  // D_n < TE_m
  if (type == 3 && grBig->type == 4) {
    unsigned int bigQueueIndex(verticesBigCount == 8 ? 3 : 2);

    // -------------------------------------------------
    // first test
    vector<short unsigned int> temp(
        grBig->vertices.begin(), grBig->vertices.begin() + bigQueueIndex + 1);
    temp.push_back(min(grBig->vertices[bigQueueIndex + 1],
                       grBig->vertices[verticesBigCount - 1]));
    temp.push_back(max(grBig->vertices[bigQueueIndex + 1],
                       grBig->vertices[verticesBigCount - 1]));
    Graph g(temp, 0, vector<bool>(false), 3, true, 0);

    if (isSubgraphOf_spherical_spherical(&g))
      return true;

    // -------------------------------------------------
    // second test
    temp = vector<short unsigned int>(grBig->vertices.begin() + bigQueueIndex,
                                      grBig->vertices.end() - 1);
    reverse(temp.begin(), temp.end());
    temp.push_back(min(grBig->vertices[bigQueueIndex - 1],
                       grBig->vertices[verticesBigCount - 1]));
    temp.push_back(max(grBig->vertices[bigQueueIndex - 1],
                       grBig->vertices[verticesBigCount - 1]));
    g = Graph(temp, 0, vector<bool>(false), 3, true, 0);

    if (isSubgraphOf_spherical_spherical(&g))
      return true;

    return false;
  }

  // ---------------------------------------------------------------------
  // E_6 < TE_6
  if (type == 4 && verticesCount == 6 && grBig->type == 4 &&
      verticesBigCount == 7) {
    vector<short unsigned int> bigBasis;

    // We must find out what the base is
    if (vertices[5] == grBig->vertices[4]) {
      bigBasis.push_back(grBig->vertices[0]);
      bigBasis.push_back(grBig->vertices[1]);
      bigBasis.push_back(grBig->vertices[6]);
      bigBasis.push_back(grBig->vertices[2]);
      bigBasis.push_back(grBig->vertices[3]);
    } else if (vertices[5] == grBig->vertices[1]) {
      bigBasis.push_back(grBig->vertices[5]);
      bigBasis.push_back(grBig->vertices[4]);
      bigBasis.push_back(grBig->vertices[6]);
      bigBasis.push_back(grBig->vertices[2]);
      bigBasis.push_back(grBig->vertices[3]);
    } else if (vertices[5] == grBig->vertices[2]) {
      bigBasis.push_back(grBig->vertices[0]);
      bigBasis.push_back(grBig->vertices[1]);
      bigBasis.push_back(grBig->vertices[6]);
      bigBasis.push_back(grBig->vertices[4]);
      bigBasis.push_back(grBig->vertices[5]);
    } else
      return false;

    if (bigBasis[0] > bigBasis[4])
      reverse(bigBasis.begin(), bigBasis.end());

    return (bigBasis ==
            vector<short unsigned int>(vertices.begin(), vertices.end() - 1));
  }

  // ---------------------------------------------------------------------
  // E_n < TE_7, TE_8
  if (type == 4 && grBig->type == 4) {
    if (vertices[verticesCount - 1] != grBig->vertices[verticesBigCount - 1])
      return false;

    // E_6, E_7 < TE7
    if (verticesCount <= 7 && verticesBigCount == 8) {
      vector<short unsigned int> bigBasis(
          grBig->vertices.begin() + 1, grBig->vertices.begin() + verticesCount);
      if (isAnSubAm(
              vector<short unsigned int>(vertices.begin(), vertices.end() - 1),
              bigBasis))
        return true;

      if (verticesCount == 7 && verticesBigCount == 8) // E_7 < TE_7
      {
        bigBasis = vector<short unsigned int>(grBig->vertices.begin(),
                                              grBig->vertices.begin() + 6);
        reverse(bigBasis.begin(), bigBasis.end());

        return (vector<short unsigned int>(vertices.begin(),
                                           vertices.end() - 1) == bigBasis);
      }

      return false;
    } else if (verticesBigCount == 9) // < TE_8
    {
      vector<short unsigned int> bigBasis(
          grBig->vertices.begin(), grBig->vertices.begin() + verticesCount - 1);

      if (vector<short unsigned int>(vertices.begin(), vertices.end() - 1) ==
          bigBasis)
        return true;

      if (verticesCount == 6) // E_5 is symmetric
      {
        reverse(bigBasis.begin(), bigBasis.end());

        return ((vector<short unsigned int>(vertices.begin(),
                                            vertices.end() - 1) == bigBasis));
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
    vector<short unsigned int> temp(grBig->vertices.begin(),
                                    grBig->vertices.end() - 1);
    if (isAnSubAm(vertices, temp))
      return true;

    // base plus l'autre extrémité
    temp = vector<short unsigned int>(grBig->vertices.begin(),
                                      grBig->vertices.end() - 2);
    temp.push_back(grBig->vertices[verticesBigCount - 1]);
    if (isAnSubAm(vertices, temp))
      return true;

    // A_3 dans le "bout du Y"
    if (verticesCount <= 3 && verticesBigCount > 3) {
      temp.clear();
      temp.push_back(min(grBig->vertices[verticesBigCount - 1],
                         grBig->vertices[verticesBigCount - 2]));
      temp.push_back(grBig->vertices[verticesBigCount - 3]);
      temp.push_back(max(grBig->vertices[verticesBigCount - 1],
                         grBig->vertices[verticesBigCount - 2]));

      return isAnSubAm(vertices, temp);
    }

    return false;
  }

  // ---------------------------------------------------------------------
  // A_n < E_m
  if (type == 0 && grBig->type == 4) {
    // A_n dans la base du E_n
    vector<short unsigned int> temp(grBig->vertices.begin(),
                                    grBig->vertices.end() - 1);
    if (isAnSubAm(vertices, temp))
      return true;

    // A_n contient la queue du E_m (début)
    temp = vector<short unsigned int>(grBig->vertices.begin(),
                                      grBig->vertices.begin() + 3);
    temp.push_back(grBig->vertices[verticesBigCount - 1]);
    if (isAnSubAm(vertices, temp))
      return true;

    // A_n contient la queue du E_m (fin)
    temp = vector<short unsigned int>(grBig->vertices.begin() + 2,
                                      grBig->vertices.end() - 1);
    temp.insert(temp.begin(), grBig->vertices[verticesBigCount - 1]);

    return isAnSubAm(vertices, temp);
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

    vector<short unsigned int> temp(grBig->vertices.begin() + 1,
                                    grBig->vertices.begin() + 4);
    reverse(temp.begin(), temp.end());

    return (vertices == temp);
  }

  // ---------------------------------------------------------------------
  // D_n < D_m
  if (type == 3 && grBig->type == 3) {
    if (verticesCount == 4) // D_4 est très symétrique
    {
      // central node
      if (vertices[1] != grBig->vertices[verticesBigCount - 3])
        return false;

      vector<short unsigned int> temp1;
      temp1.push_back(vertices[0]);
      temp1.push_back(vertices[2]);
      temp1.push_back(vertices[3]);

      vector<short unsigned int> temp2;
      temp2.push_back(grBig->vertices[verticesBigCount - 1]);
      temp2.push_back(grBig->vertices[verticesBigCount - 2]);
      temp2.push_back(grBig->vertices[verticesBigCount - 4]);

      sort(temp1.begin(), temp1.end());
      sort(temp2.begin(), temp2.end());

      return (temp1 == temp2);
    }

    // Autres D_n
    if (vertices[verticesCount - 2] != grBig->vertices[verticesBigCount - 2] ||
        vertices[verticesCount - 1] != grBig->vertices[verticesBigCount - 1])
      return false;

    vector<short unsigned int> temp1(vertices.begin(),
                                     vertices.begin() + verticesCount - 2);
    vector<short unsigned int> temp2(
        grBig->vertices.begin() + verticesBigCount - verticesCount,
        grBig->vertices.begin() + verticesBigCount - 2);

    return (temp1 == temp2);
  }

  // ---------------------------------------------------------------------
  // D_n < E_m
  if (type == 3 && grBig->type == 4) {
    // -------------------------------------------------
    // first test
    vector<short unsigned int> temp(grBig->vertices.begin(),
                                    grBig->vertices.begin() + 3);
    temp.push_back(
        min(grBig->vertices[3], grBig->vertices[verticesBigCount - 1]));
    temp.push_back(
        max(grBig->vertices[3], grBig->vertices[verticesBigCount - 1]));
    Graph g(temp, 0, vector<bool>(false), 3, true, 0);

    if (isSubgraphOf_spherical_spherical(&g))
      return true;

    // -------------------------------------------------
    // second test
    temp = vector<short unsigned int>(grBig->vertices.begin() + 2,
                                      grBig->vertices.end() - 1);
    reverse(temp.begin(), temp.end());
    temp.push_back(
        min(grBig->vertices[1], grBig->vertices[verticesBigCount - 1]));
    temp.push_back(
        max(grBig->vertices[1], grBig->vertices[verticesBigCount - 1]));
    g = Graph(temp, 0, vector<bool>(false), 3, true, 0);

    if (isSubgraphOf_spherical_spherical(&g))
      return true;

    return false;
  }

  // ---------------------------------------------------------------------
  // E_n < E_m
  if (type == 4 && grBig->type == 4) {
    if (vertices[verticesCount - 1] != grBig->vertices[verticesBigCount - 1])
      return false;

    vector<short unsigned int> temp1(vertices.begin(), vertices.end() - 1);
    vector<short unsigned int> temp2(
        grBig->vertices.begin(), grBig->vertices.begin() + verticesCount - 1);

    if (temp1 == temp2)
      return true;

    if (verticesCount == 6) // E_5 a une base symétrique
    {
      reverse(temp1.begin(), temp1.end());
      if (temp1 == temp2)
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

bool Graph::isAnSubAm(const vector<short unsigned int> &subGraphVertices,
                      const vector<short unsigned int> &bigGraphVertices) {
  vector<short unsigned int>::const_iterator itBig, itSub;
  auto it = find(bigGraphVertices.begin(), bigGraphVertices.end(),
                 subGraphVertices[0]);

  if (it == bigGraphVertices.end())
    return false;

  if (subGraphVertices.size() == 1)
    return true;

  if (subGraphVertices.size() > bigGraphVertices.size())
    return false;

  if ((it + 1) == bigGraphVertices.end() && *(it - 1) != subGraphVertices[1])
    return false;

  // ----------------------------------------------------------------------
  // recherche en avant depuis it
  if ((it + 1) != bigGraphVertices.end() && *(it + 1) == subGraphVertices[1]) {
    itBig = it;
    for (itSub = subGraphVertices.begin();
         itSub != subGraphVertices.end() && itBig != bigGraphVertices.end();
         ++itSub) {
      if (*itSub != *itBig)
        break;
      ++itBig;
    }

    return (itSub == subGraphVertices.end());
  }

  // ----------------------------------------------------------------------
  // recherche en arrière depuis it
  itBig = it;
  for (itSub = subGraphVertices.begin();
       itSub != subGraphVertices.end() && itBig >= bigGraphVertices.begin();
       ++itSub) {
    if (*itSub != *itBig)
      return false;

    --itBig;
  }

  return (itSub == subGraphVertices.end());
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
