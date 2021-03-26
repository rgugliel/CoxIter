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

#include "graphs.list.h"

GraphsList::GraphsList(size_t maxVertices,
                       vector<string> *ptr_map_vertices_indexToLabel)
    : maxVertices(maxVertices),
      graphsCount(vector<size_t>(maxVertices + 1, 0)) {
  for (unsigned int i(0); i <= maxVertices; ++i)
    graphs.push_back(GraphsListN(i, ptr_map_vertices_indexToLabel));

  totalGraphsCount = 0;
}

void GraphsList::addGraph(const std::vector<short unsigned int> &vertices,
                          const std::vector<bool> &linkableVertices,
                          const unsigned int &type, bool isSpherical,
                          const unsigned int &vertexSupp1,
                          const unsigned int &vertexSupp2,
                          const unsigned int &dataSupp) {
  size_t sizeTemp, verticesCount(vertices.size());

  if (isSpherical) {
    if (type == 1 || type == 3 || type == 4 || type == 5 || type == 7)
      verticesCount++;
  } else {
    if ((type == 1 && verticesCount >= 3 && dataSupp) || type == 2 ||
        (type == 3 && verticesCount >= 3) ||
        (type == 4 && verticesCount == 5)) // TBn, TCn, TE6
      verticesCount += 2;
    else if ((type == 0 && dataSupp) || type == 1 || type == 3 || type == 4 ||
             type == 5 ||
             (type == 6 && dataSupp)) // TAn ou TB3 ou TD4 ou TEn ou TG_2
      verticesCount++;
  }

  sizeTemp = graphs[verticesCount].size();
  graphs[verticesCount].addGraph(vertices, linkableVertices, type, isSpherical,
                                 vertexSupp1, vertexSupp2, dataSupp);
  if (sizeTemp != graphs[verticesCount].size()) {
    totalGraphsCount++;
    graphsCount[verticesCount]++;
  }
}

Graph *GraphsList::begin() {
  if (!totalGraphsCount)
    return 0;

  for (unsigned int i(1); i <= maxVertices; i++) {
    if (graphsCount[i])
      return graphs[i].begin();
  }

  return 0;
}

Graph *GraphsList::next(size_t &verticesCount, size_t &graphIndex) {
  if ((++graphIndex) < graphsCount[verticesCount])
    return graphs[verticesCount].next(graphIndex);
  else {
    for (size_t i(verticesCount + 1); i <= maxVertices; i++) {
      if (graphsCount[i]) {
        verticesCount = i;
        graphIndex = 0;
        return graphs[verticesCount].next(graphIndex);
      }
    }
  }

  return 0;
}

ostream &operator<<(ostream &o, const GraphsList &g) {
  for (const auto &graph : g.graphs)
    o << graph;

  return o;
}
