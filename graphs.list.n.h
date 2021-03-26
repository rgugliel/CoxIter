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

#ifndef GRAPHS_LIST_N_H
#define GRAPHS_LIST_N_H

#include "graph.h"

#include <algorithm>
#include <iterator>
#include <vector>

using namespace std;

/*!
 * \file graphs.list.n.h
 * \author Rafael Guglielmetti
 *
 * \class GraphsListN
 * \brief Liste des graphes d'une taille donnée
 */
class GraphsListN {
private:
  vector<Graph> graphs;       ///< Liste des graphes trouvés
  unsigned int verticesCount; ///< Nombre de sommets des graphes de la liste

  vector<string>
      *ptr_map_vertices_indexToLabel; ///< Pointeur vers la correspondance

public:
  GraphsListN(unsigned int verticesCount,
              vector<string> *ptr_map_vertices_indexToLabel);

  /*! 	\fn addGraph
   * 	\brief Ajoute un graphe à la liste
   * 	\param vertices Sommets du graphe
   * 	\param linkableVertices Sommets qui sont liables au graphes
   * 	\param type(const unsigned int &): Type du graphe (A, B, D, E, F, G, H)
   *= (0, 1, 3, 4, 5, 6, 7)
   *  \param isSpherical(bool) True si sphérique, false si euclidien
   *  \param vertexSupp1 Premier sommet supplémentaire (par exemple pour les B)
   *  \param vertexSupp2 Second sommet supplémentaire (par exemple pour les B)
   *  \param dataSupp(const unsigned int &): donnée supplémentaire (par exemple,
   *pour le G_2, le poids)
   */
  void addGraph(vector<short unsigned int> vertices,
                const vector<bool> &linkableVertices, const unsigned &type,
                bool isSpherical, const short unsigned int &vertexSupp1 = 0,
                const short unsigned int &vertexSupp2 = 0,
                const unsigned int &dataSupp = 0);

  /*!	\fn addGraphsList
   * 	\brief Concatenate another list to the current list
   *
   * 	\param gln(const GraphsListN&) Other list
   * 	\return bool (true if success)
   */
  bool addGraphsList(const GraphsListN &gln);

  /*! 	\fn size
   * 	\brief Retourne la taille de la liste de graphes
   * 	\return Taille de la liste de graphes (size_t)
   */
  size_t size() const;

  /*! \fn begin
   * \brief Renvoie un pointeur vers le premier élément de la liste
   * \return Pointeur vers l'élément (ou 0 si la liste est vide)
   */
  Graph *begin();

  /*! \fn next
   * \brief Renvoie un pointeur vers le prochain élément
   * \param graphIndex Index du graphe en cours
   * \return Pointeur vers l'élément (ou 0 si la fin de la liste est atteinte)
   */
  Graph *next(const size_t &graphIndex);

  /*!	\fn get_verticesCount
   * 	\brief Return the number of vertices of the graphs of the list
   * (this->verticesCount)
   *
   * 	\return this->verticesCount(unsigned int)
   */
  unsigned int get_verticesCount() const;

  /*!	\fn get_graphs
   * 	\brief Return the list of graphs
   *
   * 	\return this->graphs(vector< Graph >)
   */
  vector<Graph> get_graphs() const;

public:
  friend ostream &operator<<(ostream &, GraphsListN const &);
};

#endif // GRAPHS_LIST_N_H
