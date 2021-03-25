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

#include "graphs.product.h"

GraphsProduct::GraphsProduct() : rank(0) {}

vector<vector<short unsigned int>> GraphsProduct::createFootPrint() {
  vector<vector<short unsigned int>> graphsCountByType(8);

  size_t max, type;

  for (const auto &graph : graphs) {
    type = graph->type;

    // If of type G, we are interested in the weight
    max = (type == 6 && graph->isSpherical) ? graph->dataSupp
                                            : graph->vertices.size();

    for (size_t i = graphsCountByType[type].size(); i < max; i++)
      graphsCountByType[type].push_back(0);

    graphsCountByType[type][max - 1]++;
  }

  return graphsCountByType;
}

ostream &operator<<(ostream &o, const GraphsProduct &gp) {
  for (const auto &graph : gp.graphs)
    o << *graph;

  return o;
}
