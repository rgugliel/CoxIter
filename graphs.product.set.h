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

#ifndef GRAPHSPRODUCTSET_H
#define GRAPHSPRODUCTSET_H

/*!
 * \file graphs.product.set.h
 * \author Rafael Guglielmetti
 * 
 * \class GraphsProduct
 * \brief: Un produit de graphess
 * 
 * Ces produits sont utilisés pour tester la compacité
*/

#include <set>
#include <string>
#include <iterator>

#include "graph.h"
#include "graphs.product.h"

using namespace std;

class GraphsProductSet
{
	public:
		unsigned int iRank; ///< Rank of the graph
		set< Graph*, GraphPtrComp > graphs; ///< Pointeurs vers les graphes qui constituent le produit
		
	public:
		GraphsProductSet( );
		GraphsProductSet( const GraphsProduct& gp );
		
		/*!	\fn get_iVertices
		 * 	\brief Get the list of vertices of the product
		 * 
		 * 	\return vertices( vector< unsigned int > )
		 */
		vector< short unsigned int > get_iVertices( ) const;
		
		/*!	\fn b_areVerticesSubsetOf
		 * 	\brief Test if the vertices appear in another Product
		 * 	\param gp( const GraphsProductSet& ): The other product
		 * 	\return Bool
		 */
		bool b_areVerticesSubsetOf( const GraphsProductSet& gp ) const;
		
	friend ostream& operator<<( ostream &o, const GraphsProductSet& gp );
};

#endif // GRAPHSPRODUCTSET_H
