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

#include "invariantTraceField.h"

// TODO: tester compilation sans OpenMP

InvariantTraceField::InvariantTraceField(const CoxIter& ci)
: iOMPMaxThreads(omp_get_max_threads())
{
	pari_init(50000000, 2);
	initializations(ci);
}

InvariantTraceField::~InvariantTraceField()
{
	pari_close();
}

void InvariantTraceField::initializations(const CoxIter& ci)
{
	vector<vector<unsigned int>> iCoxeterMatrix(ci.get_iCoxeterMatrix());
	iVerticesCount = iCoxeterMatrix.size();
	iAdjacency = vector<vector<unsigned int>>(iVerticesCount, vector<unsigned int>());
	
	unsigned int j;
	
	for (unsigned int i(0); i < iVerticesCount; i++)
	{
		for (j = 0; j < iVerticesCount; j++)
		{
			if (iCoxeterMatrix[i][j] != 2 && i != j)
				iAdjacency[i].push_back(j);
		}
	}
	
	// -----------------------------------------------
	// Debug data
	cout << "Number of threads: " << iOMPMaxThreads << endl;
	cout << "#Vertices: " << iVerticesCount << endl;
}

void InvariantTraceField::compute()
{
	findCycles();
}

// TODO: plus carré des edges

void InvariantTraceField::findCycles()
{
	iPath = vector<vector<unsigned int>>(iOMPMaxThreads, vector<unsigned int>());
	bEdgesVisited = vector<vector<vector<bool>>>(iOMPMaxThreads, vector<vector<bool>>());
	
	#pragma omp parallel for
	for (unsigned int i = 0; i < iVerticesCount; i++)
	{
		unsigned int iThreadIndex(omp_get_thread_num());
		
		iPath[iThreadIndex].clear();
		bEdgesVisited[iThreadIndex] = vector<vector<bool>>(iVerticesCount, vector<bool>(iVerticesCount, false));
		
		findCycles(omp_get_thread_num(), i, i);
	}
}

void InvariantTraceField::findCycles(const unsigned int& iThreadIndex, const unsigned int& iRoot, const unsigned int& iFrom)
{
	iPath[iThreadIndex].push_back(iRoot); // We add the vertex to the path
	
	for (vector<unsigned int>::const_iterator it(iAdjacency[iRoot].begin()); it != iAdjacency[iRoot].end(); ++it)
	{
		// If we did not visit this edge
		if (!bEdgesVisited[iThreadIndex][iRoot][*it] && *it >= iPath[iThreadIndex][0])
		{
			if (*it == iPath[iThreadIndex][0]) // A cycle is found
			{
				if (iPath[iThreadIndex][1] < iPath[iThreadIndex][iPath[iThreadIndex].size() - 1]) // We do not to test each cycle twice
				{
					#pragma omp critical
					{
						for(vector<unsigned int>::const_iterator itSub(iPath[iThreadIndex].begin()); itSub != iPath[iThreadIndex].end(); ++itSub)
							cout << (*itSub + 1) << " ";
						cout << endl;
					}
				}
			}
			else if (find(iPath[iThreadIndex].begin(), iPath[iThreadIndex].end(), *it) == iPath[iThreadIndex].end())
			{
				bEdgesVisited[iThreadIndex][iRoot][*it] = bEdgesVisited[iThreadIndex][*it][iRoot] = true;
				findCycles(iThreadIndex, *it, iRoot);
			}
		}
	}
	
	if (iFrom != iRoot)
		bEdgesVisited[iThreadIndex][iRoot][iFrom] = bEdgesVisited[iThreadIndex][iFrom][iRoot] = false;
	
	iPath[iThreadIndex].pop_back();
}
