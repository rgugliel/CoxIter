/*
Copyright (C) 2013, 2014, 2015, 2016
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

/*
 * THIS IS AN EXPERIMENTAL MODULE!!!
 * Use at your own risk!!
 * TODO: remove
 */

/*!
 * \file invariantTraceField.h
 * \author Rafael Guglielmetti
 * 
 * \class InvariantTraceField
 * \brief To compute the invariant trace field
*/

#ifndef __INVARIANT_TRACE_FIELD_H__
#define __INVARIANT_TRACE_FIELD_H__ 1

#include <pari/pari.h>

#include "coxiter.h"

class InvariantTraceField
{
	private:
		vector<vector<unsigned int>> iAdjacency; ///< Neighbours of each vertex TODO: séparer pour chaque thread?
		unsigned int iVerticesCount; ///< Number of vertices of the graph
		
		// To find cycles
		vector<vector<unsigned int>> iPath;
		vector<vector<vector<bool>>> bEdgesVisited;
		
		const unsigned int iOMPMaxThreads;
		
	public:
		InvariantTraceField(const CoxIter& ci);
		~InvariantTraceField();
		
		void compute();
		
	private:
		void initializations(const CoxIter& ci);
		
		void findCycles();
		void findCycles(const unsigned int& iThreadIndex, const unsigned int& iRoot, const unsigned int& iFrom);
};

#endif
