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
 * \file graph.h
 * \author Rafael Guglielmetti
 *
 * \class Graph
 * \brief This class represents one graph
 */

#ifndef GRAPH_H
#define GRAPH_H

#include <algorithm>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

using namespace std;

class Graph {
public:
  vector<short unsigned int> vertices; ///< Vertices of the graph
  vector<bool>
      linkableVertices; ///< A quels sommets on peut lier le graphe { sommets }
                        ///< \ { sommets du graphes et leurs voisins)

  unsigned int type; ///< Type of the graph: A=0, B=1, ...

  unsigned int dataSupp; ///< In the case of G_2^n, the weight

  bool isSpherical; ///< True if spherical, false if euclidean

private:
  vector<string>
      *ptr_map_vertices_indexToLabel; ///< Pointeur vers la correspondance
  bool b_map_vertices_indexToLabelIsEmpty;

public:
  /*! \fn Graph
   * 	\brief Constructeur
   * 	\param vertices(const vector<short unsigned int>&) Sommets composant le
   * graphe \param ptr_map_vertices_indexToLabel(vector< string > *) Pointeur
   * vers la correspondance index --> label des sommets \param
   * linkableVertices(const vector< bool >&) Ce à quoi le graphe est liable
   * 	\param type (const unsigned int &) Type du graphe: A=0, B=1, ...
   * 	\param isSpherical(const bool& ) True si sphérique, false si euclidien
   * 	\param dataSupp(const unsigned int &) Eventuelle information
   * supplémentaire, par exemple poids pour le G2
   */
  Graph(const vector<short unsigned int> &vertices,
        vector<string> *ptr_map_vertices_indexToLabel,
        const vector<bool> &linkableVertices, const unsigned int &type,
        const bool &isSpherical, const unsigned int &dataSupp = 0);

  /*! \fn isSubgraphOf
   * 	\brief Test if a graph (spherical) is a subgraph of another graph
   *(spherical or euclidean)
   *
   * 	\param grBig(Graph*) Bigger graph
   *	\return True if the *this is a subgraph of the graph given in parameter,
   *false otherwise
   */
  bool isSubgraphOf(const Graph *grBig) const;

  /*! \fn isSubgraphOf_spherical_spherical
   * 	\brief Test if a graph (spherical) is a subgraph of another graph
   *(spherical) graph
   *
   * 	\param grBig(Graph*) Bigger graph
   *	\return True if the *this is a subgraph of the graph given in parameter,
   *false otherwise
   */
  bool isSubgraphOf_spherical_spherical(const Graph *grBig) const;

  /*! \fn isSubgraphOf_spherical_euclidean
   * 	\brief Test if a graph (spherical) is a subgraph of another graph
   *(euclidean) graph
   *
   * 	\param grBig(Graph*) Bigger graph
   *	\return True if the *this is a subgraph of the graph given in parameter,
   *false otherwise
   */
  bool isSubgraphOf_spherical_euclidean(const Graph *grBig) const;

  friend bool operator==(const Graph &g1, const Graph &g2);
  friend bool operator<(const Graph &g1, const Graph &g2);

private:
  static bool isAnSubAm(const vector<short unsigned int> &subGraphVertices,
                        const vector<short unsigned int> &bigGraphVertices);

public:
  friend ostream &operator<<(ostream &, Graph const &);
};

struct GraphPtrComp {
  bool operator()(const Graph *g1, const Graph *g2) const {
    return (*g1 < *g2);
  }
};

#endif // GRAPH_H
