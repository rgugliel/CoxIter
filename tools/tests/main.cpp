#include <iostream>

#ifdef WIN32
	#include <Windows.h>
#else
	#include <sys/time.h>
	#include <ctime>
#endif

using namespace std;

#include "tests.h"

typedef unsigned long long timestamp_t;

static timestamp_t get_timestamp( )
{
	#ifdef WIN32
		return GetTickCount( ) * 1000;
	#else
		struct timeval now;
		gettimeofday (&now, NULL);
		return  now.tv_usec + (timestamp_t)now.tv_sec * 1000000;
	#endif
}

int main( ) 
{
	Tests t;
	
	timestamp_t t0 = get_timestamp( );
	
	t.readGraphsFile( "../tests.txt" );
	t.runTests();
	
	timestamp_t t1 = get_timestamp( );
	cout << "Temps de calculs: " << ( t1 - t0 ) / 1000000.0L << endl;
	
	return 0;
}
