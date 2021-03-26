/*
Copyright (C) 2013-2017
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

#include "string.h"

void str_replace(string &str, const string &from, const string &to) {
  size_t start_pos = 0;
  while ((start_pos = str.find(from, start_pos)) != std::string::npos) {
    str.replace(start_pos, from.length(), to);
    start_pos += to.length();
  }
}

vector<string> explode(const string &separator, string source) {
  size_t found;
  vector<string> results;

  found = source.find_first_of(separator);
  while (found != string::npos) {
    if (found > 0)
      results.push_back(source.substr(0, found));
    else
      results.push_back("");

    source = source.substr(found + 1);
    found = source.find_first_of(separator);
  }

  if (source.length() > 0)
    results.push_back(source);

  return results;
}

void explode(const string &separator, string source, vector<string> &results) {
  size_t found;
  results.clear();

  found = source.find_first_of(separator);
  while (found != string::npos) {
    if (found > 0)
      results.push_back(source.substr(0, found));
    else
      results.push_back("");

    source = source.substr(found + 1);
    found = source.find_first_of(separator);
  }

  if (source.length() > 0)
    results.push_back(source);
}

void explode(const string &separator, string source, vector<int> &results) {
  size_t found;
  results.clear();

  found = source.find_first_of(separator);
  while (found != string::npos) {
    if (found > 0)
      results.push_back(stoi(source.substr(0, found)));
    else
      results.push_back(0);

    source = source.substr(found + 1);
    found = source.find_first_of(separator);
  }

  if (source.length() > 0)
    results.push_back(stoi(source));
}

void explode(const string &separator, string source,
             vector<unsigned int> &results) {
  size_t found;
  results.clear();

  found = source.find_first_of(separator);
  while (found != string::npos) {
    if (found > 0)
      results.push_back(abs(stoi(source.substr(0, found))));
    else
      results.push_back(0);

    source = source.substr(found + 1);
    found = source.find_first_of(separator);
  }

  if (source.length() > 0)
    results.push_back(abs(stoi(source)));
}

string implode(const string &separator, const vector<string> &vector) {
  ostringstream oStr;
  string res;

  copy(vector.begin(), vector.end(),
       ostream_iterator<string>(oStr, separator.c_str()));

  res = oStr.str();
  return res.substr(0, res.size() - separator.size());
}
