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

#include "arithmeticity.h"

Arithmeticity::Arithmeticity()
: iVerticesCount(0), bNotArithmetic(false), ci(0), bListCycles(false)
{
}

Arithmeticity::~Arithmeticity()
{
}

void Arithmeticity::test(CoxIter& ci_, const bool& bListCycles_)
{
	ci = &ci_;
	bListCycles = bListCycles_;
	strListCycles.clear();
	if (ci->get_iIsCocompact() != 0) // If the graph is cocompact (1) or if we don't know (-1)
	{
		strError = "GROUP COCOMPACTNESS";
		ci->set_iIsArithmetic(-1);
		return;
	}
	
	size_t i, j;
	
	iCoxeterMatrix = ci->get_iCoxeterMatrix();
	
	iVerticesCount = ci->get_iVerticesCount();
	iReferencesToLabels.clear();
	bNotArithmetic = false;
	
	// ------------------------------------------------------
	// Cycles consisting of two elements
	bool bSQRT2(false), bSQRT3(false);
	
	for (i = 0; i < iVerticesCount; i++)
	{
		iReferencesToLabels.push_back(i);
		
		for (j = 0; j < i; j++)
		{
			if (i == j)
				continue;
			
			if (iCoxeterMatrix[i][j] != 1 && iCoxeterMatrix[i][j] != 0 && iCoxeterMatrix[i][j] != 2 && iCoxeterMatrix[i][j] != 3 &&  iCoxeterMatrix[i][j] != 4 && iCoxeterMatrix[i][j] != 6)
			{
				if (ci->get_bDebug())
					cout << "\tNot arithmetic: 2*G(" << ci->get_strVertexLabel(i) << "," << ci->get_strVertexLabel(j) << ") = pi/" << iCoxeterMatrix[i][j] << endl;
				
				ci->set_iIsArithmetic(0);
				return;
			}
			
			if (iCoxeterMatrix[i][j] == 4)
				bSQRT2 = true;
			else if (iCoxeterMatrix[i][j] == 6)
				bSQRT3 = true;
			else if (iCoxeterMatrix[i][j] == 1)
			{
				string strC("4 * " + string("l") + to_string(j) + "m" + to_string(i) + "^2");
				
				auto it(lower_bound(strListCycles.begin(), strListCycles.end(), strC));
				if (it == strListCycles.end() || *it != strC)
					strListCycles.insert(it, strC);
			}
		}
	}
	
	// ------------------------------------------------------
	// Here, we know that m_{ij} \in {2,3,4,6,infty}

	// If true, the group is arithmetic
	if (!bSQRT2 && !bSQRT3 && !ci->get_bHasDottedLine())
	{
		ci->set_iIsArithmetic(1);
		return;
	}
	
	while (collapseQueues() && iVerticesCount)
		;
	
	if (iVerticesCount <= 2)
	{
		ci->set_iIsArithmetic(1);
		return;
	}
	
	testCycles();
	
	return;
}

unsigned int Arithmeticity::collapseQueues()
{
	vector< unsigned int > iVerticesToRemove;
	vector< unsigned int > iNumberNeighbours(iVerticesCount, 0); // For each vertex, number of neighbours
	unsigned int i, j, k;
	unsigned int iPreviousVertex, iCurrrentVertex;
	
	// ------------------------------------------------------
	// We count the number of neighbours of each vertex
	for (i = 0; i < iVerticesCount; i++)
	{
		for (j = 0; j < iVerticesCount; j++)
		{
			if (iCoxeterMatrix[i][j] != 2)
				iNumberNeighbours[i]++;
		}
	}
	
	// ------------------------------------------------------
	// We determine wich vertices we have to remove
	for (i = 0; i < iVerticesCount; i++)
	{
		if (iNumberNeighbours[i] == 1) // If this is a queue
		{
			iVerticesToRemove.push_back(i);
			
			// vertex next to the queue
			for (iCurrrentVertex = 0; iCurrrentVertex < iVerticesCount; iCurrrentVertex++)
			{
				if (iCoxeterMatrix[i][iCurrrentVertex] != 2 && iNumberNeighbours[iCurrrentVertex] == 2)
				{
					iVerticesToRemove.push_back(iCurrrentVertex);
					break;
				}
			}
			
			if (iCurrrentVertex == iVerticesCount) // If we cannot remove more than one vertex
				continue;
			
			iPreviousVertex = i;
			
			// Here: i ------- iCurrentVertex
			
			// We continue the path
			while (true)
			{
				for (k = 0; k < iVerticesCount; k++)
				{
					if (iCoxeterMatrix[iCurrrentVertex][k] != 2 && iNumberNeighbours[k] == 2 && k != iPreviousVertex)
					{
						iVerticesToRemove.push_back(k);
						iPreviousVertex = iCurrrentVertex;
						iCurrrentVertex = k;
						break;
					}
				}
				
				if (k == iVerticesCount) // The end of the path
					break;
			}
		}
	}
	
	// ------------------------------------------------------
	// We remove the vertices
	sort(iVerticesToRemove.begin(), iVerticesToRemove.end());
	iVerticesToRemove = vector< unsigned int >(iVerticesToRemove.begin(), unique(iVerticesToRemove.begin(), iVerticesToRemove.end()));
	
	if (iVerticesToRemove.size() == iVerticesCount) // If we remove all the vertices
	{
		iCoxeterMatrix = vector< vector< unsigned int > >(0, vector< unsigned int >(0));
		iVerticesCount = 0;
		
		return iVerticesToRemove.size();
	}
	
	reverse(iVerticesToRemove.begin(), iVerticesToRemove.end());
	for (vector< unsigned int >::const_iterator it(iVerticesToRemove.begin()); it != iVerticesToRemove.end(); ++it)
	{
		iCoxeterMatrix.erase(iCoxeterMatrix.begin() + *it);
		iReferencesToLabels.erase(iReferencesToLabels.begin() + *it);
		
		for (vector< vector< unsigned int > >::iterator itRow(iCoxeterMatrix.begin()); itRow != iCoxeterMatrix.end(); ++itRow)
			itRow->erase(itRow->begin() + *it);
	}
	
	iVerticesCount -= iVerticesToRemove.size();
	return iVerticesToRemove.size();
}

void Arithmeticity::testCycles()
{
	for (unsigned int i(0); i < iVerticesCount; i++)
	{
		iPath.clear();
		bVerticesVisited = vector< bool >(iVerticesCount, false);
		bEdgesVisited = vector< vector<bool> >(iVerticesCount, vector<bool>(iVerticesCount, false));
		findCycles(i, i);
		
		if (bNotArithmetic)
		{
			ci->set_iIsArithmetic(0);
			return;
		}
	}
	
	if (!ci->get_bHasDottedLine()) 
		ci->set_iIsArithmetic(1);
	else
		ci->set_iIsArithmetic(-1);
}

void Arithmeticity::findCycles(const unsigned int& iRoot, const unsigned int& iFrom)
{
	iPath.push_back(iRoot); // We add the vertex to the path
	
	for (unsigned int i(iPath[0]); i < iVerticesCount; i++)
	{
		// If i is a neighbour and if we did not visit this edge
		if (iCoxeterMatrix[iRoot][i] != 2 && !bEdgesVisited[iRoot][i])
		{
			if (i == iPath[0])
			{
				if (iPath[1] < iPath[ iPath.size() - 1 ]) // We do not to test each cycle twice
					testCycle();
				
				if (bNotArithmetic)
				{
					if (ci->get_bDebug())
					{
						cout << "\tNot arithmetic\n\t\tCycle: ";
						for (vector< unsigned int >::const_iterator it(iPath.begin()); it != iPath.end(); ++it) // We display the components of the cycle
							cout << (it ==  iPath.begin() ? "" : ", ") << ci->get_strVertexLabel(iReferencesToLabels[ *it ]);
						cout << endl;
					}
					
					return;
				}
			}
			else if (find(iPath.begin(), iPath.end(), i) == iPath.end())
			{
				bEdgesVisited[iRoot][i] = bEdgesVisited[i][iRoot] = true;
				findCycles(i, iRoot);
				
				if (bNotArithmetic)
					return;
			}
		}
	}
	
	if (iFrom != iRoot)
		bEdgesVisited[iRoot][iFrom] = bEdgesVisited[iFrom][iRoot] = false;
	
	iPath.pop_back();
}

void Arithmeticity::testCycle()
{
	unsigned int iPathSize(iPath.size());
	
	if (!bListCycles)
	{
		bool bNumberSQRT2Even(iCoxeterMatrix[ iPath[0] ][ iPath[iPathSize - 1] ] == 4 ? false : true), bNumberSQRT3Even(iCoxeterMatrix[ iPath[0] ][ iPath[iPathSize - 1] ] == 6 ? false : true);
		
		for (unsigned int i(1); i < iPathSize; i++)
		{
			if (iCoxeterMatrix[ iPath[i] ][ iPath[i-1] ] == 4)
				bNumberSQRT2Even = bNumberSQRT2Even ? false : true;
			else if (iCoxeterMatrix[ iPath[i] ][ iPath[i-1] ] == 6)
				bNumberSQRT3Even = bNumberSQRT3Even ? false : true;
			else if (iCoxeterMatrix[ iPath[i] ][ iPath[i-1] ] == 1) // Because of the dotted line we cannot say anything for this cycle
				return;
		}
		
		if (iCoxeterMatrix[ iPath[0] ][ iPath[iPathSize - 1] ] == 1) // Because of the dotted line we cannot say anything for this cycle
			return;
		
		bNotArithmetic = !bNumberSQRT2Even || !bNumberSQRT3Even;
	}
	else
	{
		unsigned int iPathSize(iPath.size());
	
		int iNumberDotted(iCoxeterMatrix[ iPath[0] ][ iPath[iPathSize - 1] ] == 1 ? 1 : 0);
		int iNumber2(iCoxeterMatrix[ iPath[0] ][ iPath[iPathSize - 1] ] == 0 ? 1 : 0);
		int iNumberSQRT2(iCoxeterMatrix[ iPath[0] ][ iPath[iPathSize - 1] ] == 4 ? 1 : 0);
		int iNumberSQRT3(iCoxeterMatrix[ iPath[0] ][ iPath[iPathSize - 1] ] == 6 ? 1 : 0);
		
		for (unsigned int i(1); i < iPathSize; i++)
		{
			if (iCoxeterMatrix[ iPath[i] ][ iPath[i-1] ] == 0)
				iNumber2++;
			else if (iCoxeterMatrix[ iPath[i] ][ iPath[i-1] ] == 4)
				iNumberSQRT2++;
			else if (iCoxeterMatrix[ iPath[i] ][ iPath[i-1] ] == 6)
				iNumberSQRT3++;
			else if (iCoxeterMatrix[ iPath[i] ][ iPath[i-1] ] == 1) // Because of the dotted line we cannot say anything for this cycle
				iNumberDotted++;
		}
		
		if (iNumberDotted == 0)
			bNotArithmetic = (iNumberSQRT2%2) || (iNumberSQRT3%2);
		else
		{
			string strTemp;
			
			iNumber2 += iNumberDotted;
			
			if (iNumberSQRT2 > 1)
				iNumber2 += ((iNumberSQRT2 % 2) ? iNumberSQRT2 - 1 : iNumberSQRT2) / 2;
			
			if (iNumber2)
				strTemp += (strTemp == "" ? "" : " * ") + string("2^") + to_string(iNumber2);
			
			if (iNumberSQRT3 > 1)
				strTemp += (strTemp == "" ? "" : " * ") + string("3^") + to_string(((iNumberSQRT2 % 2) ? iNumberSQRT2 - 1 : iNumberSQRT2) / 2);
			
			if (iNumberSQRT2 % 2)
				strTemp += (strTemp == "" ? "" : " * ") + string("Sqrt[2]");
			
			if (iNumberSQRT3 % 2)
				strTemp += (strTemp == "" ? "" : " * ") + string("Sqrt[3]");
			
			for (unsigned int i(1); i < iPathSize; i++)
			{
				if (iCoxeterMatrix[ iPath[i] ][ iPath[i-1] ] == 1)
					strTemp += (strTemp == "" ? "" : " * ") + string("l") + to_string(min(iReferencesToLabels[iPath[i]], iReferencesToLabels[iPath[i-1]])) + "m" + to_string(max(iReferencesToLabels[iPath[i]], iReferencesToLabels[iPath[i-1]]));
			}
			
			if (iCoxeterMatrix[ iPath[0] ][ iPath[iPathSize - 1] ] == 1)
				strTemp += (strTemp == "" ? "" : " * ") + string("l") + to_string(min(iReferencesToLabels[iPath[0]], iReferencesToLabels[iPath[iPathSize - 1]])) + "m" + to_string(max(iReferencesToLabels[iPath[0]], iReferencesToLabels[iPath[iPathSize - 1]]));
			
			
			auto it(lower_bound(strListCycles.begin(), strListCycles.end(), strTemp));
			if (it == strListCycles.end() || *it != strTemp)
				strListCycles.insert(it, strTemp);
		}
	}
}

vector< string > Arithmeticity::get_strListCycles()
{
	return strListCycles;
}

string Arithmeticity::get_strError()
{
	return strError;
}
