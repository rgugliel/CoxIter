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

GraphsList::GraphsList(size_t verticesCount,
                       vector<string> *ptr_map_vertices_indexToLabel)
    : verticesCount(verticesCount),
      graphsCount(vector<size_t>(verticesCount + 1, 0)) {
  for (unsigned int i(0); i <= verticesCount; ++i)
    graphs.push_back(GraphsListN(i, ptr_map_vertices_indexToLabel));

  totalGraphsCount = 0;
}

void GraphsList::addGraph(const std::vector<short unsigned int> &vertices,
                          const std::vector<bool> &bVerticesLinkable,
                          const unsigned int &type, bool isSpherical,
                          const unsigned int &iVertexSupp1,
                          const unsigned int &iVertexSupp2,
                          const unsigned int &iDataSupp) {
  size_t sizeTemp, iVCount(vertices.size());

  if (isSpherical) {
    if (type == 1 || type == 3 || type == 4 || type == 5 || type == 7)
      iVCount++;
  } else {
    if ((type == 1 && iVCount >= 3 && iDataSupp) || type == 2 ||
        (type == 3 && iVCount >= 3) ||
        (type == 4 && iVCount == 5)) // TBn, TCn, TE6
      iVCount += 2;
    else if ((type == 0 && iDataSupp) || type == 1 || type == 3 || type == 4 ||
             type == 5 ||
             (type == 6 && iDataSupp)) // TAn ou TB3 ou TD4 ou TEn ou TG_2
      iVCount++;
  }

  sizeTemp = graphs[iVCount].size();
  graphs[iVCount].addGraph(vertices, bVerticesLinkable, type, isSpherical,
                           iVertexSupp1, iVertexSupp2, iDataSupp);
  if (sizeTemp != graphs[iVCount].size()) {
    totalGraphsCount++;
    graphsCount[iVCount]++;
  }
}

Graph *GraphsList::begin() {
  if (!totalGraphsCount)
    return 0;

  for (unsigned int i(1); i <= verticesCount; i++) {
    if (graphsCount[i])
      return graphs[i].begin();
  }

  return 0;
}

Graph *GraphsList::next(size_t &iVCount, size_t &graphIndex) {
  if ((++graphIndex) < graphsCount[iVCount])
    return graphs[iVCount].next(graphIndex);
  else {
    for (size_t i(iVCount + 1); i <= verticesCount; i++) {
      if (graphsCount[i]) {
        iVCount = i;
        graphIndex = 0;
        return graphs[iVCount].next(graphIndex);
      }
    }
  }

  return 0;
}

ostream &operator<<(ostream &o, const GraphsList &g) {
  unsigned int iMax(g.graphs.size());
  for (unsigned int i(0); i < iMax; ++i)
    o << g.graphs[i];

  return o;
}
