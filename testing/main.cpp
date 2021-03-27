/*
Copyright (C) 2013, 2014, 2015, 2016, 2017
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

#include <chrono>
#include <iostream>

using namespace std;

#include "tests.h"

int main() {
  chrono::time_point<chrono::system_clock> timeStart, timeEnd;

  cout << "Test suite for CoxIter\n" << endl;

  Tests t;
  if (!t.readGraphsFile("../tests.txt")) {
    cout << "Error\n\tCannot read ../tests.txt" << endl;
    cout << "Aborting" << endl;
    return 0;
  }

  cout << "Number of graphs that will be read: " << t.get_testsCount() << endl;

  timeStart = chrono::system_clock::now();

  if (!t.runTests()) {
    cout << "Error\n\t" << t.get_error() << endl;
    cout << "Aborting" << endl;
    return 0;
  }

  timeEnd = chrono::system_clock::now();
  cout << "Computation time: "
       << chrono::duration<double, milli>(timeEnd - timeStart).count() / 1000
       << "s\n"
       << endl;

  return 0;
}
