/*
Copyright (C) 2013, 2014
Rafael Guglielmetti, rafael.guglielmetti@unifr.ch
*/

/*
This file is part of CoxIter and AlVin.

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
 * \file string.h
 * \brief Quelques fonctions relatives au string
 * \author Rafael Guglielmetti
*/

#ifndef __STRING_H__
#define __STRING_H__ 1

#include <string>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <iterator>

using namespace std;

/**
	* \fn void str_replace( string &str, const string &from, const string &to )
	* \brief Rechercher remplacer
	*
	* \param[out] str Chaîne dans laquelle se fait la recherche/remplacement
	* \param[in] from Chaîne à rechercher
	* \param[in] to Ce par quoi on remplace
	* \return void
*/
void str_replace( string &str, const string &from, const string &to );

/**
	* \fn vector<string> explode( string source, const string &separator );
	* \brief Équivalent de la fonction explode de PHP
	*
	* \param[in] source Chaîne de caractères sur laquelle on travaille
	* \param[in] separator Séparateur
	* \return vector<string>: Tableau contenant les résultats

	Prend une string ainsi qu'un séparateur et découpe la chaîne vers un tableau
*/
vector<string> explode( const string &separator, string source );

/**
	* \fn vector<string> explode( string source, const string &separator, vector<string> &results );
	* \brief Équivalent de la fonction explode en C++
	*
	* \param[in] source Chaîne de caractères sur laquelle on travaille
	* \param[in] separator Séparateur
	* \param[out] results Tableau auquel sont ajoutés les résultats (par référence)
	* \return vector<string>: void

	Prend une string ainsi qu'un séparateur et découpe la chaîne vers un tableau
*/
void explode( const string &separator, string source, vector<string> &results );

/**
	* \fn vector<string> explode( string source, const string &separator, vector<int> &results );
	* \brief Équivalent de la fonction explode en C++
	*
	* \param[in] source Chaîne de caractères sur laquelle on travaille
	* \param[in] separator Séparateur
	* \param[out] results Tableau auquel sont ajoutés les résultats (par référence)
	* \return vector<string>: void

	Prend une string ainsi qu'un séparateur et découpe la chaîne vers un tableau
*/
void explode( const string &separator, string source, vector<int> &results );

/**
	* \fn vector<string> explode( string source, const string &separator, vector<unsigned int> &results );
	* \brief Équivalent de la fonction explode en C++
	*
	* \param[in] source Chaîne de caractères sur laquelle on travaille
	* \param[in] separator Séparateur
	* \param[out] results Tableau auquel sont ajoutés les résultats (par référence)
	* \return vector<string>: void

	Prend une string ainsi qu'un séparateur et découpe la chaîne vers un tableau
*/
void explode( const string &separator, string source, vector<unsigned int> &results );

/*!
 * 	\fn implode
 * 	\brief Implode function (as in the PHP language)
 * 
 * 	\param strSeparator( const string & ) Separator
 * 	\param strVector( const vector< string >& ) Vector to implode
 * 	\return Imploded string
 */
string implode( const string& strSeparator, const vector< string >& strVector );

/*!
 * 	\fn implode
 * 	\brief Implode function (as in the PHP language)
 * 
 * 	\param strSeparator( const string & ) Separator
 * 	\param iVector( const vector< int >& ) Vector to implode
 * 	\return Imploded string
 */
template <typename T>
typename std::enable_if<std::is_arithmetic<T>::value, string>::type implode( const string& strSeparator, const vector< T >& iVector )
{
	vector< string > strVector;
	for( typename vector< T >::const_iterator it( iVector.begin( ) ); it != iVector.end( ); ++it )
		strVector.push_back( to_string( *it ) );
	
	return implode( strSeparator, strVector );
}

int string_to_int( const string &strNumber );
double string_to_double( const string &strNumber );

#endif