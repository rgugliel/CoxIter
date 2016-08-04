#ifndef MULTITOOLS_H
#define MULTITOOLS_H

#include <string>
#include <iostream>
#include <fstream>
#include <iterator>

#include "arithmeticity.number.h"
#include "../../../tools/regexp.h"
#include "../../../tools/string.h"
#include "../../coxiter.h"
#include "../../arithmeticity.h"

using namespace std;

class CycleVector
{
	public:
		ArithmeticityNumber an;
		unsigned int iVector;
		
		CycleVector( ArithmeticityNumber an, unsigned int iVector );
		
		bool operator==( const CycleVector& cv );
		friend bool operator<( const CycleVector& cv1, const CycleVector& cv2 );
};

class MultiTools
{
	private:
		string strError;
		
		vector< string > strGraphsFilename;
		
		unsigned int iMax;
		vector< unsigned int > iWeights;
		vector< vector< ArithmeticityNumber > > anGramMatix;
	public:
		bool bReadGraphsList( const string& strList );
		vector< string > get_strGraphsFilename( ) const;
		
		vector< CycleVector > cyclesVectors;
		void createCycles_main( ); // TODO
		void createCycles_enumerate( const unsigned int& iStartingIndex );
		void createCycles_tovector( );
};

#endif // MULTITOOLS_H
