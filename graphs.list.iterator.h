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
 * \file graphs.list.iterator.h
 * \author Rafael Guglielmetti
 *
 * \class GraphsListIterator graphs.list.iterator.h
 * \brief Permet de parcourir une liste de graphes
 * Utilisation (pour un GraphsList *gl)<br>
 * 	GraphsListIterator iter(gl);<br>
        while (iter.ptr)<br>
        {<br>
                cout << "iter; " << *iter.ptr << endl;<br>
                ++iter;<br>
        }
*/

#ifndef GRAPHS_LIST_ITERATOR_H
#define GRAPHS_LIST_ITERATOR_H

#include "graphs.list.h"

class GraphsListIterator {
private:
  size_t verticesCount; ///< Nombre de sommets du graphe courant
  size_t graphIndex;    ///< Index du graphe courant

  GraphsList *graphsList; ///< Pointeur vers la liste de graphes

  unsigned int verticesCountMax;
  bool limitMaxVertices;

public:
  /*!
   * \brief Constructeur
   * \param gl vers la liste de graphes à considérer
   */
  GraphsListIterator(GraphsList *gl);

  GraphsListIterator(const GraphsListIterator &gl);

  /*!
   * \brief Constructeur
   * \param gl vers la liste de graphes à considérer
   * \param verticesCountMin(const unsigned int&) Nombre de sommets où l'on
   * commence \param verticesCountMax(const unsigned int&) Nombre de sommets où
   * l'on s'arrête
   */
  GraphsListIterator(GraphsList *gl, const unsigned int &verticesCountMin,
                     const unsigned int &verticesCountMax = 0);

  GraphsListIterator();

  /*!
   *	\brief Déplace ptr vers le prochain graphe
   *
   *	Si tous la fin de la liste est atteinte, ptr est mis à 0
   *
   * 	\return Pointeur vers le prochain graphe (ou 0 si la fin de la liste a
   *été atteinte)
   */
  Graph *next();

  /*!	\fn operator++()
   * 	\overload ++
   *	\brief Appel à next()
   */
  GraphsListIterator &operator++();

public:
  Graph *ptr; ///< Pointeur vers le graphe courant
};

#endif // GRAPHS_LIST_ITERATOR_H
