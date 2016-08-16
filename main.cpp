#include <iostream>

#include "app.h"
#include "coxiter.h"

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
	
	app.run( );
	
	return 0;
}
