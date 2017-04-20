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

InvariantTraceField::InvariantTraceField(const CoxIter& ci)
{
	pari_init(50000000, 2);
}

InvariantTraceField::~InvariantTraceField()
{
	pari_close();
}

void InvariantTraceField::initializations(const CoxIter& ci)
{
	vector<vector<unsigned int>> iCoxeterMatrix(ci.get_iCoxeterMatrix());
	iVerticesCount = iCoxeterMatrix.size();
	iAdjacency = vector<vector<short unsigned int>>(iVerticesCount, vector<short unsigned int>());
	
	unsigned int j;
	
	for (unsigned int i(0); i < iVerticesCount; i++)
	{
		for (j = i + 1; j < iVerticesCount; j++)
		{
			if (iCoxeterMatrix[i][j] != 2)
				iAdjacency[i].push_back(j);
		}
	}
}
