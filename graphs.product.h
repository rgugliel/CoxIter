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
 * \file graphs.product.h
 * \author Rafael Guglielmetti
 *
 * \class GraphsProduct
 * \brief: Un produit de graphs
 *
 * Ces produits sont utilisés de manière volatile dans le programme. Les graphes
 * sont ajoutés dans le std::vector au fur et à mesure. Au contraire,
 * GraphsProductSet est utilisé pour garder les produits de manière persistante;
 * les graphes sont stockés dans un std::set
 */

#ifndef GRAPHS_PRODUCT_H
#define GRAPHS_PRODUCT_H

#include <iterator>
#include <string>
#include <vector>

#include "graph.h"

using namespace std;

class GraphsProduct {
public:
  vector<Graph *> graphs; ///< Pointers to the graphs of the product
  unsigned int rank;      ///< Rank of the product

public:
  GraphsProduct();

  /*! \fn createFootPrint
   * 	\brief Create a footprint for the graph.
   *  This footprint is used as a key for the map CoxIter.graphsProductsCount
   *  \return The footprint of the graph
   */
  vector<vector<short unsigned int>> createFootPrint();

  friend ostream &operator<<(ostream &o, const GraphsProduct &gp);
};

#endif // GRAPHS_PRODUCT_H
