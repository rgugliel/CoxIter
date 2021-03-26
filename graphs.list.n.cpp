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

#include "graphs.list.n.h"

GraphsListN::GraphsListN(unsigned int verticesCount,
                         vector<string> *ptr_map_vertices_indexToLabel)
    : verticesCount(verticesCount),
      ptr_map_vertices_indexToLabel(ptr_map_vertices_indexToLabel) {}

size_t GraphsListN::size() const { return graphs.size(); }

Graph *GraphsListN::next(const size_t &graphIndex) {
  return &graphs[graphIndex];
}

Graph *GraphsListN::begin() {
  if (graphs.empty())
    return 0;

  return &graphs[0];
}

void GraphsListN::addGraph(vector<short unsigned int> vertices,
                           const vector<bool> &linkableVertices,
                           const unsigned int &type, bool isSpherical,
                           const short unsigned int &vertexSupp1,
                           const short unsigned int &vertexSupp2,
                           const unsigned int &dataSupp) {
  unsigned int short temp, i;

  if (isSpherical) {
    if (type == 0) // An
    {
      // on veut une fois les A_n (et pas une fois pour chaque sens de lecture)
      if (vertices.front() > vertices.back())
        reverse(vertices.begin(), vertices.end());
    } else if (type == 1) // Bn
    {
      vertices.push_back(vertexSupp1);
    } else if (type == 3) // Dn
    {
      if (verticesCount > 4) {
        temp = vertices[verticesCount - 2];
        vertices[verticesCount - 2] = min(temp, vertexSupp1);
        vertices.push_back(max(temp, vertexSupp1));
      } else {
        if (vertices[0] > vertices[2]) {
          temp = vertices[0];
          vertices[0] = vertices[2];
          vertices[2] = temp;
        }

        if (vertexSupp1 < vertices[0]) {
          vertices.push_back(vertices[2]);
          vertices[2] = vertices[0];
          vertices[0] = vertexSupp1;
        } else if (vertexSupp1 < vertices[2]) {
          vertices.push_back(vertices[2]);
          vertices[2] = vertexSupp1;
        } else
          vertices.push_back(vertexSupp1);
      }
    } else if (type == 4) // En
    {
      if (verticesCount == 6) {
        if (vertices.front() > vertices.back())
          reverse(vertices.begin(), vertices.end());
      }
      vertices.push_back(vertexSupp1);
    } else if (type == 5) // Fn
    {
      vertices.push_back(vertexSupp1);
      if (vertices.front() > vertices.back())
        reverse(vertices.begin(), vertices.end());
    } else if (type == 7) // Hn
    {
      vertices.push_back(vertexSupp1);
    }
  } else {
    if (type == 0 && dataSupp) {
      int iMinIndex(0);
      unsigned int iMinValue(vertices[0]);

      vertices.push_back(vertexSupp1);
      vector<short unsigned int> verticesTemp(vertices);

      // on répère le sommet avec l'indice le plus petit
      for (i = 1; i < verticesCount; i++) {
        if (vertices[i] < iMinValue) {
          iMinValue = vertices[i];
          iMinIndex = i;
        }
      }

      // et le sens de parcours
      int iDirection(
          vertices[!iMinIndex ? (verticesCount - 1) : (iMinIndex - 1)] >
                  vertices[iMinIndex == (int)(verticesCount - 1)
                               ? 0
                               : (iMinIndex + 1)]
              ? 1
              : -1);

      // on réordonne
      for (unsigned int j(0); j < verticesCount; j++) {
        vertices[j] = verticesTemp[iMinIndex];

        iMinIndex += iDirection;
        if (iMinIndex == -1)
          iMinIndex = verticesCount - 1;
        else if (iMinIndex == (int)verticesCount)
          iMinIndex = 0;
      }
    } else if (type == 1) // TBn
    {
      if (verticesCount == 4) // TB3
      {
        i = vertices[1];
        temp = max(vertices[0], vertices[2]);
        vertices[1] = min(vertices[0], vertices[2]);
        vertices[2] = temp;
        vertices[0] = i;

        vertices.insert(vertices.begin(), vertexSupp1);
      } else // autres \tilde Bn
      {
        temp = min(vertices[verticesCount - 3], vertexSupp1);
        vertices.push_back(max(vertices[verticesCount - 3], vertexSupp1));
        vertices[verticesCount - 3] = temp;

        vertices.insert(vertices.begin(), vertexSupp2);
      }
    } else if (type == 2) // \tilde Cn
    {
      vertices.insert(vertices.begin(), vertexSupp1);
      vertices.push_back(vertexSupp2);

      if (vertices[0] > vertices[verticesCount - 1])
        reverse(vertices.begin(), vertices.end());
    } else if (type == 3) // \tilde Dn
    {
      if (verticesCount >= 6) {
// le vecteur verticesBase contient la base (i.e. sans les 4 extrémités)
#ifdef _MSC_VER
        vector<short unsigned int> verticesBase;
        temp = verticesCount - 3;
        for (i = 1; i < temp; i++)
          verticesBase.push_back(vertices[i]);
#else
        vector<short unsigned int> verticesBase(
            vertices.begin() + 1, vertices.begin() + verticesCount -
                                      3); // ne marche pas sous Visual Studio
                                          // 2010 & 2012 malheureusement
#endif

        // on regarde si on doit modifier l'ordre
        if (verticesBase[0] > verticesBase[verticesCount - 5]) {
          reverse(verticesBase.begin(), verticesBase.end());

          // ajout des 4 extrémités
          verticesBase.push_back(min(vertices[0], vertexSupp1));
          verticesBase.push_back(max(vertices[0], vertexSupp1));
          verticesBase.insert(verticesBase.begin(),
                              max(vertices[verticesCount - 3], vertexSupp2));
          verticesBase.insert(verticesBase.begin(),
                              min(vertices[verticesCount - 3], vertexSupp2));
        } else {
          // ajout des 4 extrémités
          verticesBase.push_back(min(vertices[verticesCount - 3], vertexSupp2));
          verticesBase.push_back(max(vertices[verticesCount - 3], vertexSupp2));
          verticesBase.insert(verticesBase.begin(),
                              max(vertices[0], vertexSupp1));
          verticesBase.insert(verticesBase.begin(),
                              min(vertices[0], vertexSupp1));
        }

        vertices = verticesBase;
      } else // le TD4, "+" est traité différemment
      {
        vector<short unsigned int> verticesTemp(4, 0);
        verticesTemp[0] = vertices[0];
        verticesTemp[1] = vertices[2];
        verticesTemp[2] = vertexSupp1;
        verticesTemp[3] = vertexSupp2;
        sort(verticesTemp.begin(), verticesTemp.end());
        verticesTemp.push_back(vertices[1]);
        vertices = verticesTemp;
      }
    } else if (type == 4) // \tilde En
    {
      if (verticesCount == 7) // TE6
      {
        // TODO: refaire l'encodage de ce graphe et modifer la fonction
        // bIsSubgraphOf_spherical_euclidean?
        vector<short unsigned int> verticesTemp;
        unsigned int iMin(min(min(vertices[1], vertices[3]), vertexSupp1));

        if (iMin == vertices[1]) {
          verticesTemp.push_back(vertices[0]);
          verticesTemp.push_back(vertices[1]);

          if (min(vertices[3], vertexSupp1) == vertices[3]) {
            verticesTemp.push_back(vertices[3]);
            verticesTemp.push_back(vertices[4]);
            verticesTemp.push_back(vertexSupp1);
            verticesTemp.push_back(vertexSupp2);
          } else {
            verticesTemp.push_back(vertexSupp1);
            verticesTemp.push_back(vertexSupp2);
            verticesTemp.push_back(vertices[3]);
            verticesTemp.push_back(vertices[4]);
          }
        } else if (iMin == vertices[3]) {
          verticesTemp.push_back(vertices[4]);
          verticesTemp.push_back(vertices[3]);

          if (min(vertices[1], vertexSupp1) == vertices[1]) {
            verticesTemp.push_back(vertices[1]);
            verticesTemp.push_back(vertices[0]);
            verticesTemp.push_back(vertexSupp1);
            verticesTemp.push_back(vertexSupp2);
          } else {
            verticesTemp.push_back(vertexSupp1);
            verticesTemp.push_back(vertexSupp2);
            verticesTemp.push_back(vertices[1]);
            verticesTemp.push_back(vertices[0]);
          }
        } else {
          verticesTemp.push_back(vertexSupp2);
          verticesTemp.push_back(vertexSupp1);

          if (min(vertices[1], vertices[3]) == vertices[1]) {
            verticesTemp.push_back(vertices[1]);
            verticesTemp.push_back(vertices[0]);
            verticesTemp.push_back(vertices[3]);
            verticesTemp.push_back(vertices[4]);
          } else {
            verticesTemp.push_back(vertices[3]);
            verticesTemp.push_back(vertices[4]);
            verticesTemp.push_back(vertices[1]);
            verticesTemp.push_back(vertices[0]);
          }
        }

        verticesTemp.push_back(vertices[2]);
        vertices = verticesTemp;
      } else if (verticesCount == 8) // TE7
      {
        if (vertices.front() >
            vertices.back()) // la base du \tilde E7 est symétrique
          reverse(vertices.begin(), vertices.end());

        vertices.push_back(vertexSupp1);
      } else // TE8
      {
        vertices.push_back(vertexSupp1);
      }
    } else if (type == 5) {
      vertices.push_back(vertexSupp1);
    } else if (type == 6 && dataSupp) // \tilde G_2
    {
      vertices.push_back(vertexSupp1);
    }
  }

  Graph g(vertices, ptr_map_vertices_indexToLabel, linkableVertices, type,
          isSpherical, dataSupp);

  auto it(lower_bound(graphs.begin(), graphs.end(), g));
  if (it == graphs.end() || !(*it == g))
    graphs.insert(it, g);
}

bool GraphsListN::addGraphsList(const GraphsListN &gln) {
  if (verticesCount != gln.get_verticesCount())
    return false;

  vector<Graph> gr(gln.get_graphs());
  graphs.insert(graphs.end(), gr.begin(), gr.end());

  return true;
}

unsigned int GraphsListN::get_verticesCount() const { return verticesCount; }

vector<Graph> GraphsListN::get_graphs() const { return graphs; }

ostream &operator<<(ostream &o, const GraphsListN &g) {
  o << "\tGraphs of rank " << g.verticesCount << endl;

  for (const auto &graph : g.graphs)
    o << graph;

  return o;
}
