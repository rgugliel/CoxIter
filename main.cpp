/*
Copyright (C) 2013-2017
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

#include <iostream>

#include "app.h"
#include "coxiter.h"

using namespace std;

int main(int argc, char **argv) {
  // -------------------------------------------------------
  App app;

  // lecture des paramètres donnés au programme
  if (!app.readMainParameters(argc, argv))
    return 0;

  if (app.outFilenameBasis == "")
    app.bCoutFile = false;

  // si la sortie standard est redirigée dans un fichier (-cf)
  if (app.bCoutFile)
    cout << "Output is redirected to " << app.outFilenameBasis << ".output"
         << endl;

  app.run();

  return 0;
}
