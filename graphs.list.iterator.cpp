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

#include "graphs.list.iterator.h"

GraphsListIterator::GraphsListIterator(GraphsList *gl)
    : limitMaxVertices(false), graphsList(gl),
      verticesCountMax(
          0) // The value 0 has no effect because of limitVerticesMax(false)
{
  ptr = graphsList->begin();
  if (ptr) {
    verticesCount = ptr->vertices.size();
    graphIndex = 0;
  }
}

GraphsListIterator::GraphsListIterator(const GraphsListIterator &gl)
    : verticesCount(gl.verticesCount), graphIndex(gl.graphIndex),
      graphsList(gl.graphsList), verticesCountMax(gl.verticesCountMax),
      limitMaxVertices(gl.limitMaxVertices), ptr(gl.ptr) {}

GraphsListIterator::GraphsListIterator(GraphsList *gl,
                                       const unsigned int &verticesCountMin,
                                       const unsigned int &verticesCountMax)
    : verticesCountMax(verticesCountMax), limitMaxVertices(false) {
  if (verticesCountMax && verticesCountMin <= verticesCountMax)
    limitMaxVertices = true;

  graphsList = gl;

  if (graphsList->graphs.size() < verticesCountMin)
    throw("Graphs of this size don't exist");

  // ---------------------------------------------------
  // We look fot the first graph in the list
  unsigned int i(verticesCountMin);
  while (i < graphsList->graphs.size() && !graphsList->graphs[i].size())
    i++;

  if (i < graphsList->graphs.size()) {
    ptr = graphsList->graphs[i].begin();
    verticesCount = ptr->vertices.size();
    graphIndex = 0;
  } else
    ptr = 0;
}

GraphsListIterator::GraphsListIterator() {}

Graph *GraphsListIterator::next() {
  ptr = graphsList->next(verticesCount, graphIndex);

  if (limitMaxVertices && verticesCount > verticesCountMax)
    return 0;

  return ptr;
}

GraphsListIterator &GraphsListIterator::operator++() {
  if (!next())
    ptr = 0;

  return *this;
}
