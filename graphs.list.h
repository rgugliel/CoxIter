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

/*!
 * \file graphs.list.h
 * \author Rafael Guglielmetti
 * 
 * \class GraphsList
 * \brief Liste des graphes
*/

#ifndef GRAPHS_LIST_H
#define GRAPHS_LIST_H

#include "graphs.list.n.h"

#include <vector>

using namespace std;

class GraphsList
{
	public:
		/*! 	\fn GraphsList
		 * 	\brief Crée la liste de graphes
		 * 	\param iVerticesCount Nombre de sommets du graphe lu
		 * 	\param ptr_map_vertices_indexToLabel( vector< string > * ) Pointeur vers la correspondance index --> label des sommets
		 */
		GraphsList( size_t iVerticesCount, vector< string > *ptr_map_vertices_indexToLabel );
		
		/*!
		 * 	\fn addGraph
		 * 	\brief Ajoute un graphe
		 * 	
		 * 	\param iVertices( const vector< unsigned int > & ): tableau contenant les sommets
		 * 	\param bVerticesLinkable( const vector< bool > & ): sommets qui sont liables (ou non) au graphe
		 * 	\param iType( const unsigned int & ): Type du graphe (A, B, D, E, F, G, H) = (0, 1, 3, 4, 5, 6, 7)
		 * 	\param bSpherical( bool ): true si graphe sphérique, false sinon
		 * 	\param iVertexSupp1( const unsigned int & ): éventuellement, premier sommet supplémentaire (dans le cas du Dn, par exemple) // TODO: valeur par défaut meilleure que 0?
		 * 	\param iVertexSupp2( const unsigned int & ): éventuellement, deuxième sommet supplémentaire (dans le cas du Dn, par exemple)
		 * 	\param iDataSupp( const unsigned int & ): donnée supplémentaire (par exemple, pour le G_2, le poids)
		 */
		void addGraph( const vector<short unsigned int> &iVertices, const vector< bool > &bVerticesLinkable, const unsigned int &iType, bool bSpherical, const unsigned int &iVertexSupp1 = 0, const unsigned int &iVertexSupp2 = 0, const unsigned int &iDataSupp = 0 );
		
		/*!	\fn begin
		 * 	\brief Retourne un pointeur sur le premier graphe de la liste 
		 * 
		 * 	\return Pointeur sur le graphe (Graph *) ou 0 si la liste est vide
		 */
		Graph* begin( );
		
		/*!	\fn next
		 * 	\brief Retourne un pointeur sur l'élément suivant de la liste
		 * 
		 * 	\param iVCount( size_t & ): Nombre de sommets du graphe actuel
		 * 	\param iGraphIndex( size_t & ): Index du graphe actuel (i.e. position dans la liste des graphes de taille iVCount)
		 * 
		 * 	\return Pointeur sur le graphe (Graph *) ou 0 si la fin de la liste est atteinte
		 */
		Graph* next( size_t &iVCount, size_t &iGraphIndex );
		
	public: // Remark: this is public for read-only purpose!
		vector< GraphsListN > graphs; ///< List of list of graphs (by number of vertices)
		vector< size_t > graphsCount; ///< Number of graphs (by number of vertices)
		
		size_t iGraphsCount; ///< Total number of graphs
		size_t iVerticesCount; ///< Maximum number of vertices in the graphs
		
	public:
		friend ostream& operator<<( ostream& , GraphsList const & );
};

#endif // GRAPHS_LIST_H
