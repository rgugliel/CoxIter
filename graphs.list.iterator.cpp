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

#include "graphs.list.iterator.h"

GraphsListIterator::GraphsListIterator( GraphsList *gl )
: bLimitVerticesMax( false ), 
  graphsList( gl ),
  iVerticesCountMax( 0 ) // The value 0 has no effect because of bLimitVerticesMax( false )
{
	ptr = graphsList->begin( );
	if( ptr )
	{
		iVCount = ptr->iVertices.size( );
		iGraphIndex = 0;
	}
}

GraphsListIterator::GraphsListIterator( GraphsList *gl, const unsigned int& iVerticesCountMin, const unsigned int& iVerticesCountMax )
: iVerticesCountMax( iVerticesCountMax ), bLimitVerticesMax( false )
{
	if( iVerticesCountMax && iVerticesCountMin <= iVerticesCountMax )
		bLimitVerticesMax = true;
	
	graphsList = gl;
	
	if( graphsList->graphs.size( ) < iVerticesCountMin )
		throw( "Graphs of this size don't exist" );
	
	// ---------------------------------------------------
	// We look fot the first graph in the list
	unsigned int i( iVerticesCountMin );
	while( i < graphsList->graphs.size( ) && !graphsList->graphs[i].size( ) )
		i++;
	
	if( i < graphsList->graphs.size( ) )
	{
		ptr = graphsList->graphs[i].begin( );
		iVCount = ptr->iVertices.size( );
		iGraphIndex = 0;
	}
	else
		ptr = 0;
}

GraphsListIterator::GraphsListIterator( )
{
}

Graph* GraphsListIterator::next( )
{
	ptr = graphsList->next( iVCount, iGraphIndex );

	if( bLimitVerticesMax && iVCount > iVerticesCountMax )
		return 0;
	
	return ptr;
}

GraphsListIterator& GraphsListIterator::operator++( )
{
	if( !next( ) )
		ptr = 0;
	
	return *this;
}
