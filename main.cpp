#include <iostream>

#include "app.h"

using namespace std;

int main( int argc, char **argv )
{
	// -------------------------------------------------------
	App app;
	
	// lecture des paramètres donnés au programme
	if( !app.bReadMainParameters( argc, argv ) )
		return 0;
	
	if( app.strOutFilenameBasis == "" )
		app.bCoutFile = false;
	
	// si la sortie standard est redirigée dans un fichier (-cf)
	if( app.bCoutFile )
		cout << "Output is redirected to " << app.strOutFilenameBasis << ".output" << endl;
	
	if( app.strInFilename == "" )
	{
		cout << "No input file given" << endl;
		return 0;
	}
	
	app.run( );
	
	return 0;
}