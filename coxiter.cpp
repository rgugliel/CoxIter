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

#include "coxiter.h"

CoxIter::CoxIter()
    : checkCocompactness(false), checkCofiniteness(false), bCoutFile(false),
      debug(false), isGramMatrixFieldKnown(false),
      isGrowthSeriesComputed(false), hasBoldLine(false), hasDottedLine(false),
      hasDottedLineWithoutWeight(0), bWriteInfo(false), isGraphExplored(false),
      isGraphsProductsComputed(false), bUseOpenMP(true),
      brEulerCaracteristic(0), graphsList_spherical(nullptr),
      graphsList_euclidean(nullptr), dimension(0), euclideanMaxRankFound(0),
      sphericalMaxRankFound(0), isDimensionGuessed(false),
      fVectorAlternateSum(0), isArithmetic(-1), isCocompact(-2),
      isFiniteCovolume(-2), verticesAtInfinityCount(0), verticesCount(0),
      outCout(0), sBufOld(0), strError(""),
      strOuputMathematicalFormat("generic") {
#ifndef _OPENMP
  this->bUseOpenMP = false;
#endif
}

CoxIter::CoxIter(const vector<vector<unsigned int>> &matrix,
                 const unsigned int &dimension)
    : checkCocompactness(false), checkCofiniteness(false), bCoutFile(false),
      isGramMatrixFieldKnown(false), isGraphExplored(false),
      isGraphsProductsComputed(false), isGrowthSeriesComputed(false),
      hasBoldLine(false), hasDottedLine(false), hasDottedLineWithoutWeight(0),
      bWriteInfo(false), debug(false), bUseOpenMP(true),
      brEulerCaracteristic(0), graphsList_spherical(nullptr),
      graphsList_euclidean(nullptr), dimension(dimension),
      euclideanMaxRankFound(0), sphericalMaxRankFound(0),
      isDimensionGuessed(false), fVectorAlternateSum(0), isCocompact(-1),
      isFiniteCovolume(-1), verticesAtInfinityCount(0), verticesCount(0),
      outCout(0), sBufOld(0), strError(""), strOuputMathematicalFormat("") {
  verticesCount = matrix.size();

  initializations();

  coxeterMatrix = matrix;

  maximalSubgraphRank = dimension ? dimension : verticesCount;

#ifndef _OPENMP
  this->bUseOpenMP = false;
#endif
}

CoxIter::~CoxIter() {
  if (graphsList_spherical)
    delete graphsList_spherical;

  if (graphsList_euclidean)
    delete graphsList_euclidean;

  // if cout is redirected to a file
  if (bCoutFile) {
    outCout->close();
    cout.rdbuf(sBufOld); // we restore the cout
  }
}

bool CoxIter::bRunAllComputations() {
  if (!coxeterMatrix.size())
    return false;

  if (!isGraphExplored)
    exploreGraph();

  if (!isGraphsProductsComputed)
    computeGraphsProducts();

  if (!bEulerCharacteristicFVector())
    return false;

  if (checkCofiniteness)
    checkCovolumeFiniteness();

  if (checkCocompactness)
    isGraphCocompact();

  return true;
}

#ifndef _COMPILE_WITHOUT_REGEXP_
bool CoxIter::parseGraph(istream &streamIn) {
  string strLine;
  PCRERegexp regexp;
  PCREResult regexpRes;

  // loops variable, first vertice, second vertice, weight, number of vertices,
  // index of the current row
  unsigned int i, i1, i2, i3, verticesFileCount, rowIndex(1);

  vector<unsigned int>::const_iterator it;

  vector<unsigned int> orders; // orders found

  // ---------------------------------------------------------------------------
  // Reading the number of vertices and, eventually, dimension
  if (getline(streamIn, strLine)) {
    if (regexp.preg_match_all("([[:digit:]]+)[[:space:]]?([[:digit:]]*)",
                              strLine, regexpRes) == 1) {
      verticesFileCount = verticesCount = stoi(regexpRes[1][0]);
      dimension = regexpRes[2][0] != "" ? stoi(regexpRes[2][0]) : 0;
    } else {
      strError = "First line with number of vertices missing";
      return false;
    }
  } else {
    strError = "EMPTY_FILE";
    return false;
  }

  // ---------------------------------------------------------------------------
  // first line
  if (!getline(streamIn, strLine)) {
    strError = "EMPTY_FILE";
    return false;
  }

  // names of the vertices
  if (regexp.preg_match_all("^vertices labels:[[:space:]]?([[:alnum:]-_ ]+)$",
                            strLine, regexpRes)) {
    vector<string> strVL(explode(" ", regexpRes[1][0]));
    if (strVL.size() != verticesFileCount) {
      strError = "VERTICES_LABEL_COUNT";
      return false;
    }

    for (i = 0; i < verticesFileCount; i++) {
      map_vertices_labelToIndex[strVL[i]] = i;
      map_vertices_indexToLabel.push_back(strVL[i]);
    }

    if (map_vertices_labelToIndex.size() != verticesFileCount) {
      strError = "VERTICES_LABEL_COUNT";
      return false;
    }

    if (!getline(streamIn, strLine)) {
      strError = "EMPTY_FILE";
      return false;
    }
  } else {
    for (i = 0; i < verticesFileCount; i++) {
      map_vertices_labelToIndex[to_string(i + 1)] = i;
      map_vertices_indexToLabel.push_back(to_string(i + 1));
    }
  }

  bool bRemoveDottedEdges(false); // If we want to remove dotted edges

  // ---------------------------------------------------------------------------
  // removed vertices
  vector<unsigned int> verticesShift(verticesFileCount,
                                     0); // Shifts for the removed vertices
  unsigned int truncCount(0);

  if (strVertices.size()) // If we want to specify a subset of the vertices
  {
    auto strAllVertices(map_vertices_indexToLabel);
    sort(strAllVertices.begin(), strAllVertices.end());

    set_difference(strAllVertices.begin(), strAllVertices.end(),
                   strVertices.begin(), strVertices.end(),
                   std::back_inserter(strVerticesRemove));
    sort(strVerticesRemove.begin(), strVerticesRemove.end());
  }

  for (vector<string>::const_iterator itStr(strVerticesRemove.begin());
       itStr != strVerticesRemove.end(); ++itStr) {
    if (*itStr == "dotted" &&
        map_vertices_labelToIndex.find("dotted") ==
            map_vertices_labelToIndex.end()) // remove dotted edges?
    {
      bRemoveDottedEdges = true;
      continue;
    }

    if (map_vertices_labelToIndex.find(*itStr) ==
        map_vertices_labelToIndex.end()) {
      strError = "This vertex does not exist: " + *itStr;
      return false;
    }

    truncCount++;

    for (i = map_vertices_labelToIndex[*itStr]; i < verticesFileCount; i++)
      verticesShift[i]++;
  }
  verticesCount -= truncCount;

  // ---------------------------------------------------------------------------
  // initializations
  initializations(); // now that we know the real number of vertices

  // ---------------------------------------------------------------------------
  // reading the graph
  do {
    regexpRes.clear();

    // Usual row: "first vertice" "second vertice" "weight"
    if (regexp.preg_match_all(
            "([[:alnum:]_-]+)[[:space:]]([[:alnum:]_-]+)[[:space:]]([[:digit:]]"
            "+)([[:space:]]+#[[:space:]]*([^\n]+))?",
            strLine, regexpRes)) {
      if (map_vertices_labelToIndex.find(regexpRes[1][0]) ==
          map_vertices_labelToIndex.end()) {
        strError = "The following vertex is unknown: " + regexpRes[1][0];
        return false;
      }

      if (map_vertices_labelToIndex.find(regexpRes[2][0]) ==
          map_vertices_labelToIndex.end()) {
        strError = "The following vertex is unknown: " + regexpRes[2][0];
        return false;
      }

      i1 = map_vertices_labelToIndex[regexpRes[1][0]];
      i2 = map_vertices_labelToIndex[regexpRes[2][0]];
      i3 = stoi(regexpRes[3][0]);

      if (i3 == 1 && bRemoveDottedEdges)
        i3 = 2;

      rowIndex++;

      // Removed vertex?
      if (binary_search(strVerticesRemove.begin(), strVerticesRemove.end(),
                        regexpRes[1][0]) ||
          binary_search(strVerticesRemove.begin(), strVerticesRemove.end(),
                        regexpRes[2][0]))
        continue;

      // on tient compte du décalage lié à la troncation
      i1 -= verticesShift[i1];
      i2 -= verticesShift[i2];

      // on garde le poids (pour corps engendré par les coefficients de la
      // matrice de Gram)
      if (find(orders.begin(), orders.end(), i3) == orders.end())
        orders.push_back(i3);

      if (i3 == 1) // Weight of the dotted line given?
      {
        if (regexpRes.size() > 5) {
          unsigned int index(linearizationMatrix_index(min(i1, i2), max(i1, i2),
                                                       verticesCount));
          strWeights[index] = regexpRes[5][0];
        } else
          hasDottedLineWithoutWeight = 1;
      }

      // si on avait déjà cette arête avec un ordre différent
      if (coxeterMatrix[i1][i2] != 2 && coxeterMatrix[i1][i2] != i3) {
        strError = "Edge has multiple orders (" + regexpRes[1][0] + "," +
                   regexpRes[2][0] + ")";
        return false;
      }

      coxeterMatrix[i1][i2] = i3;
      coxeterMatrix[i2][i1] = i3;

      if (i3 == 1) // dotted
        hasDottedLine = true;
      else if (i3 == 0)
        hasBoldLine = true;
    } else if (strLine != "") {
      if (bWriteInfo)
        cout << "Unread line (incorrect format): "
             << "#" << strLine << "#" << rowIndex << endl;

      rowIndex++;
      continue;
    }
  } while (getline(streamIn, strLine));

  // ---------------------------------------------------------------------------
  // Labels and co
  auto v_ItL(map_vertices_indexToLabel);
  map_vertices_indexToLabel.clear();
  map_vertices_labelToIndex.clear();

  unsigned int j(0);
  for (unsigned int i(0); i < verticesFileCount; i++) {
    if ((!i && !verticesShift[i]) ||
        (i && verticesShift[i] == verticesShift[i - 1])) {
      map_vertices_indexToLabel.push_back(v_ItL[i]);
      map_vertices_labelToIndex[v_ItL[i]] = j++;
    }
  }

  maximalSubgraphRank = dimension ? dimension : verticesCount;

  // ---------------------------------------------------------------------------
  // some information
  if (bWriteInfo) {
    cout << "Reading graph: " << endl;
    cout << "\tNumber of vertices: " << verticesCount << endl;
    cout << "\tDimension: " << (dimension ? to_string(dimension) : "?") << endl;

    cout << "\tVertices: ";
    for (vector<string>::const_iterator itStr(
             map_vertices_indexToLabel.begin());
         itStr != map_vertices_indexToLabel.end(); ++itStr)
      cout << (itStr != map_vertices_indexToLabel.begin() ? ", " : "")
           << *itStr;
    cout << endl;
  }

  // ---------------------------------------------------------------------------
  // Field generated by the entries of the Gram matrix
  for (it = orders.begin(); it != orders.end(); ++it) {
    if (*it == 1) // dotted line
      break;
    else if (*it == 4)
      strGramMatrixField +=
          (strGramMatrixField == "" ? "sqrt(2)" : ", sqrt(2)");
    else if (*it == 5)
      strGramMatrixField +=
          (strGramMatrixField == "" ? "sqrt(5)" : ", sqrt(5)");
    else if (*it == 6)
      strGramMatrixField +=
          (strGramMatrixField == "" ? "sqrt(3)" : ", sqrt(3)");
    else if (*it >= 7) // le static_cast est là pour VC++
      strGramMatrixField +=
          (strGramMatrixField == ""
               ? "cos(pi/" + to_string(static_cast<long long>(*it)) + ")"
               : ", cos(pi/" + to_string(static_cast<long long>(*it)) + ")");
  }

  if (it == orders.end()) {
    strGramMatrixField =
        strGramMatrixField == "" ? "Q" : ("Q[" + strGramMatrixField + "]");
    isGramMatrixFieldKnown = true;

    if (bWriteInfo)
      cout << "\tField generated by the entries of the Gram matrix: "
           << strGramMatrixField << endl;
  } else {
    if (bWriteInfo)
      cout << "\tField generated by the entries of the Gram matrix: ?" << endl;
  }

  if (bWriteInfo)
    cout << "File read\n" << endl;

  return true;
}

bool CoxIter::bReadGraphFromFile(const string &strInputFilename) {
  // ---------------------------------------------------------------------------
  // try to open the file
  ifstream fileIn(strInputFilename.c_str());
  if (fileIn.fail()) {
    strError = "Cannot open file";
    return false;
  }

  if (!parseGraph(fileIn))
    return false;

  fileIn.close();

  return true;
}
#endif

void CoxIter::initializations() {
  // ------------------------------------------------------
  // de initializations
  if (graphsList_spherical)
    delete graphsList_spherical;

  if (graphsList_euclidean)
    delete graphsList_euclidean;

  graphsProductsCount_spherical.clear();
  graphsProductsCount_euclidean.clear();

  factorials.clear();
  powersOf2.clear();

  isGraphExplored = false;
  isGraphsProductsComputed = false;

  // ------------------------------------------------------
  // initializations
  coxeterMatrix = vector<vector<unsigned int>>(
      verticesCount, vector<unsigned int>(verticesCount, 2));
  bVerticesVisited = vector<bool>(verticesCount, false);
  bEdgesVisited =
      vector<vector<bool>>(verticesCount, vector<bool>(verticesCount, false));

  graphsList_spherical =
      new GraphsList(verticesCount, &map_vertices_indexToLabel);
  graphsList_euclidean =
      new GraphsList(verticesCount, &map_vertices_indexToLabel);

  graphsProductsCount_euclidean =
      vector<map<vector<vector<short unsigned int>>, unsigned int>>(
          verticesCount + 1,
          map<vector<vector<short unsigned int>>, unsigned int>());
  graphsProductsCount_spherical =
      vector<map<vector<vector<short unsigned int>>, unsigned int>>(
          verticesCount + 1,
          map<vector<vector<short unsigned int>>, unsigned int>());

  // ------------------------------------------------------------
  // sauvegarde de quelques calculs
  factorials = vector<mpz_class>(verticesCount + 2, 1);
  powersOf2 = vector<mpz_class>(verticesCount + 2, 1);
  for (unsigned int i(1); i <= verticesCount + 1; i++) {
    factorials[i] = factorials[i - 1] * (long int)i;
    powersOf2[i] = mpz_class(2) * powersOf2[i - 1];
  }
}

bool CoxIter::bWriteGraph(const string &strOutFilenameBasis) {
  if (strOutFilenameBasis == "") {
    strError = "No file specified for writing the graph";
    return false;
  }

  map_vertices_labels_create();

  string strFilename(strOutFilenameBasis + ".coxiter");
  ofstream out(strFilename.c_str());
  if (!out.is_open()) {
    strError = "Cannot open the file for writing the graph";
    return false;
  }

  out << verticesCount << (dimension ? " " + to_string(dimension) : "") << endl;
  out << "vertices labels: ";
  for (vector<string>::const_iterator it(map_vertices_indexToLabel.begin());
       it != map_vertices_indexToLabel.end(); ++it)
    out << (it == map_vertices_indexToLabel.begin() ? "" : " ") << *it;
  out << endl;

  for (unsigned int i(0); i < verticesCount; i++) {
    for (unsigned int j(0); j < i; j++) {
      if (coxeterMatrix[i][j] != 2)
        out << map_vertices_indexToLabel[j] << " "
            << map_vertices_indexToLabel[i] << " " << coxeterMatrix[i][j]
            << endl;
    }
  }

  out.close();

  return true;
}

void CoxIter::map_vertices_labels_create() {
  if (map_vertices_indexToLabel.size())
    return; // nothing to do

  for (unsigned int i(0); i < verticesCount; i++) {
    map_vertices_labelToIndex[to_string(i + 1)] = i;
    map_vertices_indexToLabel.push_back(to_string(i + 1));
  }
}

void CoxIter::map_vertices_labels_reinitialize() {
  map_vertices_labelToIndex.clear();
  map_vertices_indexToLabel.clear();

  for (unsigned int i(0); i < verticesCount; i++) {
    map_vertices_labelToIndex[to_string(i + 1)] = i;
    map_vertices_indexToLabel.push_back(to_string(i + 1));
  }
}

bool CoxIter::bWriteGraphToDraw(const string &strOutFilenameBasis) {
  unsigned int i, j;

  map_vertices_labels_create();

  // ----------------------------------------------------------------------
  // ouverture du fichier
  if (strOutFilenameBasis == "") {
    strError = "No file specified for writing the graph";
    return false;
  }

  string strFilename(strOutFilenameBasis + ".graphviz");
  ofstream out(strFilename.c_str());
  if (!out.is_open()) {
    strError = "Cannot open the file for writing the graph";
    return false;
  }

  // ----------------------------------------------------------------------
  // écriture à proprement parler
  out << "graph G { " << endl;
  for (i = 0; i < verticesCount; i++)
    out << "\t\"" << map_vertices_indexToLabel[i] << "\";" << endl;

  for (i = 0; i < verticesCount; i++) {
    for (j = i + 1; j < verticesCount; j++) {
      if (coxeterMatrix[i][j] > 3)
        out << "\t\"" << map_vertices_indexToLabel[i] << "\" -- \""
            << map_vertices_indexToLabel[j] << "\" [label=\""
            << coxeterMatrix[i][j] << "\"];" << endl;
      else if (coxeterMatrix[i][j] > 2)
        out << "\t\"" << map_vertices_indexToLabel[i] << "\" -- \""
            << map_vertices_indexToLabel[j] << "\";" << endl;
      else if (coxeterMatrix[i][j] == 1)
        out << "\t\"" << map_vertices_indexToLabel[i] << "\" -- \""
            << map_vertices_indexToLabel[j] << "\" [style=dotted];" << endl;
      else if (coxeterMatrix[i][j] == 0)
        out << "\t\"" << map_vertices_indexToLabel[i] << "\" -- \""
            << map_vertices_indexToLabel[j] << "\" [label=\"inf\"];" << endl;
    }
  }

  out << "}";
  out.close();

  return true;
}

void CoxIter::exploreGraph() {
  vector<short unsigned int> vertices;
  short unsigned int i, j, k, l;

  if (!verticesCount)
    throw(string("CoxIter::exploreGraph: No graph given"));

  if (isGraphExplored)
    return;

  // -------------------------------------------------------------------
  // pour chaque sommet, on cherche toutes les chaînes qui partent, ce qui donne
  // les An, Bn, Dn, En, Hn
  path.clear();

  for (i = 0; i < verticesCount; i++) {
    path.clear();
    bVerticesVisited = vector<bool>(verticesCount, false);
    bEdgesVisited =
        vector<vector<bool>>(verticesCount, vector<bool>(verticesCount, false));

    DFS(i, i);
  }

  // -------------------------------------------------------------------
  // recherche des A_1, G_2^k avec k >= 4, F_4
  vector<bool> bVerticesLinkable, bVerticesLinkableTemp;
  for (i = 0; i < verticesCount; i++) {
    bVerticesLinkable = vector<bool>(verticesCount, true);
    for (j = 0; j < verticesCount; j++) {
      if (coxeterMatrix[i][j] != 2)
        bVerticesLinkable[j] = false;
    }

    // ajout du sommet (A_1)
    bVerticesLinkable[i] = false;
    graphsList_spherical->addGraph(vector<short unsigned int>(1, i),
                                   bVerticesLinkable, 0, true);

    // on regarde si on trouve avec ce sommet: Gn, TA1 ,TC2
    for (j = 0; j < verticesCount; j++) {
      if (coxeterMatrix[i][j] >= 4 || !coxeterMatrix[i][j]) {
        // ------------------------------------------------------------------
        // G2 et TA1
        if (i < j) {
          vertices.clear();
          vertices.push_back(i);
          vertices.push_back(j);

          bVerticesLinkableTemp = bVerticesLinkable;
          for (k = 0; k < verticesCount; k++) {
            if (coxeterMatrix[j][k] != 2)
              bVerticesLinkableTemp[k] = false;
          }

          if (coxeterMatrix[i][j]) // ici, c'est un graphe sphérique
            graphsList_spherical->addGraph(vertices, bVerticesLinkableTemp, 6,
                                           true, 0, 0, coxeterMatrix[i][j]);
          else // ici, graphe euclidien (TA1)
            graphsList_euclidean->addGraph(vertices, bVerticesLinkableTemp, 0,
                                           false, 0, 0, 0);
        }

        // ------------------------------------------------------------------
        // TC2 = [ 4, 4 ]
        if (coxeterMatrix[i][j] == 4) {
          for (k = 0; k < verticesCount; k++) {
            if (coxeterMatrix[k][j] == 4 && i != k &&
                coxeterMatrix[i][k] == 2) {
              bVerticesLinkableTemp = bVerticesLinkable;
              for (l = 0; l < verticesCount; l++) {
                if (coxeterMatrix[k][l] != 2)
                  bVerticesLinkableTemp[l] = false;
                if (coxeterMatrix[j][l] != 2)
                  bVerticesLinkableTemp[l] = false;
              }

              graphsList_euclidean->addGraph(vector<short unsigned int>(1, j),
                                             bVerticesLinkableTemp, 2, false, i,
                                             k);
            }
          }
        }
      }
    }
  }

  isGraphExplored = true;
}

void CoxIter::DFS(unsigned int root, unsigned int from) {
  // -------------------------------------------------------------------
  // initializations
  bool subcall(false); // to know if we call DFS

  unsigned int i;

  /*
   * 	We don't want cycles
   * 		We mark neighbours of iFrom as visited (to avoir cycles)
   * 		We stock this in verticesVisited to restore it at the end
   */
  vector<unsigned int> visitedVertices;

  if (root != from) {
    for (i = 0; i < verticesCount; i++) {
      if (coxeterMatrix[from][i] != 2) {
        if (!bVerticesVisited[i])
          visitedVertices.push_back(i);

        bVerticesVisited[i] = true;
      }
    }
  } else
    bVerticesVisited[root] = true; // obviously...

  // -------------------------------------------------------------------
  // DFS
  path.push_back(root); // we add root to the path

  for (i = 0; i < verticesCount; i++) {
    // if we have an edge AND if i was not traversed AND if the edge (root,i)
    // was not traversed
    if (coxeterMatrix[root][i] == 3 && !bVerticesVisited[i] &&
        !bEdgesVisited[root][i]) {
      bEdgesVisited[root][i] = bEdgesVisited[i][root] = true;
      subcall = true;
      DFS(i, root);
    }
  }

  // -------------------------------------------------------------------
  // un-initializations
  for (const auto &visitedVertex : visitedVertices)
    bVerticesVisited[visitedVertex] = false;

  bVerticesVisited[root] = false;

  if (from != root) // If it is a recursive call
    bEdgesVisited[root][from] = bEdgesVisited[from][root] = false;

  // If DFS was not called, then the path is maximal
  if (!subcall)
    addGraphsFromPath();

  path.pop_back();
}

void CoxIter::addGraphsFromPath() {
  // sommets que l'on ne peut pas lier au graphe (n sommets);
  vector<bool> bVerticesLinkable(verticesCount, true);

  // sommets que l'on ne peut pas lier au graphe (n-1 sommets), sommets que l'on
  // ne peut pas lier au graphe (n-2 sommets)
  vector<bool> bVerticesLinkable_0_nMin1(verticesCount, true),
      bVerticesLinkable_0_nMin2(verticesCount, true);

  // sommets que l'on ne peut pas lier au graphe (1 --> n), sommets que l'on ne
  // peut pas lier au graphe (1 --> n-1), , sommets que l'on ne peut pas lier au
  // graphe (2 --> n)
  vector<bool> bVerticesLinkable_1_n(verticesCount, true),
      bVerticesLinkable_1_nMin1(verticesCount, true),
      bVerticesLinkable_2_n(verticesCount, true);

  // vecteur temporaire
  vector<bool> bVerticesLinkableTemp, bVerticesLinkableTempTemp;

  // chemin en cours de construction
  vector<short unsigned int> pathTemp;

  // i, j, k, l: variables de boucles
  short unsigned int i, j, k, l, max(path.size()), iOrder;

  for (i = 0; i < max; i++) {
    pathTemp.push_back(path[i]);

    // --------------------------------------------------------------------
    // mise à jour des voisinages occupés
    for (j = 0; j < verticesCount; j++) {
      if (coxeterMatrix[path[i]][j] != 2)
        bVerticesLinkable[j] = false;

      if (i >= 1 && coxeterMatrix[path[i]][j] != 2)
        bVerticesLinkable_1_n[j] = false;

      if (i >= 2 && coxeterMatrix[path[i]][j] != 2)
        bVerticesLinkable_2_n[j] = false;

      if (i >= 1 && coxeterMatrix[path[i - 1]][j] != 2)
        bVerticesLinkable_0_nMin1[j] = false;

      if (i >= 2 && coxeterMatrix[path[i - 1]][j] != 2)
        bVerticesLinkable_1_nMin1[j] = false;

      if (i >= 2 && coxeterMatrix[path[i - 2]][j] != 2)
        bVerticesLinkable_0_nMin2[j] = false;
    }

    // --------------------------------------------------------------------
    // An
    if (i != 0) // on ajoute pas les sommets
      graphsList_spherical->addGraph(pathTemp, bVerticesLinkable, 0, true);

    // --------------------------------------------------------------------
    // TAn, n >= 2
    if (i >= 1) {
      bVerticesLinkable_1_nMin1[pathTemp[i - 1]] = false;

      for (j = 0; j < verticesCount; j++) {
        // si la partie centrale ne pose pas de problème ET qu'on est lié à
        // chaque extrémité
        if (bVerticesLinkable_1_nMin1[j] && !bVerticesLinkable[j] &&
            coxeterMatrix[j][pathTemp[0]] == 3 &&
            coxeterMatrix[j][pathTemp[i]] == 3) {
          // mise à jour des linkables avec le somme trouvé
          bVerticesLinkableTemp = bVerticesLinkable;
          for (k = 0; k < verticesCount; k++) {
            if (coxeterMatrix[k][j] != 2)
              bVerticesLinkableTemp[k] = false;
          }

          graphsList_euclidean->addGraph(
              pathTemp, bVerticesLinkableTemp, 0, false, j, 0,
              1); // TODO OPTIMIZATION modifier ce 1 (relatif à une meilleure
                  // valeur que "0" par défaut pour les dernières variables)
        }
      }
    }

    // --------------------------------------------------------------------
    // Dn, TBn, TDn
    if (i >= 2) {
      bVerticesLinkable_0_nMin2[pathTemp[i - 2]] = false;

      // --------------------------------------------------------------------
      // Dn, TBn (n >= 4), TDn (n >= 4)
      // on regarde les voisins de l'avant dernier sommet
      for (j = 0; j < verticesCount; j++) {
        // si y'a une arête ET si on est pas déjà dans le chemin ET si pas
        // interdit ET pas lien entre deux extrémités
        if (coxeterMatrix[j][pathTemp[i - 1]] == 3 &&
            (pathTemp[i - 1] != j && pathTemp[i] != j) &&
            bVerticesLinkable_0_nMin2[j] &&
            coxeterMatrix[j][pathTemp[i]] == 2) {
          bVerticesLinkableTemp = bVerticesLinkable;
          for (k = 0; k < verticesCount; k++) {
            if (coxeterMatrix[k][j] != 2)
              bVerticesLinkableTemp[k] = false;
          }

          graphsList_spherical->addGraph(pathTemp, bVerticesLinkableTemp, 3,
                                         true, j); // Dn

          // --------------------------------------------------------------------
          // ici, on va tenter de trouver un TD_n (n >= 4) (i.e. prolonger par
          // une arrête au 2ème sommet)
          for (k = 0; k < verticesCount; k++) {
            if (k != pathTemp[0] && k != j && k != pathTemp[i] &&
                coxeterMatrix[pathTemp[1]][k] == 3 &&
                bVerticesLinkable_2_n[k] &&
                coxeterMatrix[pathTemp[0]][k] == 2 &&
                coxeterMatrix[k][j] == 2) {
              bVerticesLinkableTempTemp = bVerticesLinkableTemp;
              for (l = 0; l < verticesCount; l++) {
                if (coxeterMatrix[k][l] != 2)
                  bVerticesLinkableTempTemp[l] = false;
              }

              graphsList_euclidean->addGraph(
                  pathTemp, bVerticesLinkableTempTemp, 3, false, k, j); // TDn
            }
          }

          // --------------------------------------------------------------------
          // ici, on va tenter de trouver un TB_n (n >= 4) (i.e. prolonger par
          // une arête de poids 4 à gauche)
          for (k = 0; k < verticesCount; k++) {
            if (coxeterMatrix[pathTemp[0]][k] == 4 &&
                bVerticesLinkable_1_n[k] && coxeterMatrix[j][k] == 2) {
              bVerticesLinkableTempTemp = bVerticesLinkableTemp;
              for (l = 0; l < verticesCount; l++) {
                if (coxeterMatrix[k][l] != 2)
                  bVerticesLinkableTempTemp[l] = false;
              }

              graphsList_euclidean->addGraph(pathTemp,
                                             bVerticesLinkableTempTemp, 1,
                                             false, j, k, 1); // TBn
            }
          }
        }

        // --------------------------------------------------------------------
        // TB3
        if (i == 2) {
          if (coxeterMatrix[pathTemp[1]][j] == 4 &&
              coxeterMatrix[pathTemp[0]][j] == 2 &&
              coxeterMatrix[pathTemp[2]][j] == 2) {
            bVerticesLinkableTemp = bVerticesLinkable;
            for (k = 0; k < verticesCount; k++) {
              if (coxeterMatrix[k][j] != 2)
                bVerticesLinkableTemp[k] = false;
            }

            graphsList_euclidean->addGraph(pathTemp, bVerticesLinkableTemp, 1,
                                           false, j); // TB3
          }
        }
      }
    }

    // --------------------------------------------------------------------
    // E6, E7, E8, TE6, TE7, TE8
    if (i >= 4 && i <= 7) {
      AnToEn_AnToTEn(pathTemp, bVerticesLinkable);
    }

    // --------------------------------------------------------------------
    // Bn, F4 et Hn, TG2, TCn et TF4
    if (i >= 1) {
      for (j = 0; j < verticesCount;
           j++) // on regarde si on peut prolonger la chaîne de 1 avec une arête
                // de poids 4
      {
        /*
         * 	Le premier paquet de conditions donne:
         * 		Prolonger par une arrête de poids 4 --> Cn ; sphérique
         * 		Prolonger par une arrête de poids 5 si (2 ou 3 sommets)
         * --> (H3 ou H4) ; sphérique Prolonger par une arrête de poids 6 --> [
         * 3, 6 ] ; euclidien
         */
        if (((coxeterMatrix[path[i]][j] == 4) ||
             (coxeterMatrix[path[i]][j] == 5 && (i == 1 || i == 2)) ||
             (coxeterMatrix[path[i]][j] == 6 && i == 1)) &&
            bVerticesLinkable_0_nMin1[j]) {
          iOrder = coxeterMatrix[path[i]][j];
          bVerticesLinkableTemp = bVerticesLinkable;
          for (k = 0; k < verticesCount; k++) {
            if (coxeterMatrix[k][j] != 2)
              bVerticesLinkableTemp[k] = false;
          }

          if (coxeterMatrix[path[i]][j] < 6) // sphérique
            graphsList_spherical->addGraph(pathTemp, bVerticesLinkableTemp,
                                           (iOrder == 4 ? 1 : 7), true, j);
          else
            graphsList_euclidean->addGraph(pathTemp, bVerticesLinkableTemp, 6,
                                           false, j, 0, 1);

          auto bVerticesLinkableTemp_bck(
              bVerticesLinkableTemp); // Contains info for: pathTemp + j

          // ------------------------------------------
          // on va tenter de prolonger cela en un TCn, n \geq 3
          if (coxeterMatrix[path[i]][j] == 4) {
            for (k = 0; k < verticesCount; k++) {
              if (coxeterMatrix[k][pathTemp[0]] == 4 && k != j &&
                  bVerticesLinkable_1_n[k] && coxeterMatrix[k][j] == 2) {
                // Additional info for vertex k
                for (l = 0; l < verticesCount; l++) {
                  if (coxeterMatrix[k][l] != 2)
                    bVerticesLinkableTemp[l] = false;
                }

                graphsList_euclidean->addGraph(pathTemp, bVerticesLinkableTemp,
                                               2, false, k, j);

                bVerticesLinkableTemp =
                    bVerticesLinkableTemp_bck; // Restoring to info of pathTemp
                                               // and j
              }
            }
          }

          // ------------------------------------------
          // ici, on a un B3, que l'on va tenter de prolonger en F4 ou un B4 que
          // l'on va tenter de prolonger en un TF4
          if ((i == 1 || i == 2) && coxeterMatrix[path[i]][j] == 4)
            B3ToF4_B4ToTF4(bVerticesLinkable_0_nMin1, pathTemp, j);

        } // if (((coxeterMatrix[ path[i] ][j] == 4) || (coxeterMatrix[
          // path[i] ][j] == 5 && (i == 1 || i == 2)) || (coxeterMatrix[
          // path[i] ][j] == 6 && i == 1)) && bVerticesLinkable_0_nMin1[j])
      }
    }
  }
}

void CoxIter::AnToEn_AnToTEn(const vector<short unsigned int> &pathTemp,
                             const vector<bool> &bVerticesLinkable) {
  unsigned int pathSize(pathTemp.size());

  /*
   * 	2 pour le sommets 3 (i.e. cas sphérique: E6, E7, E8 ou car euclidien
   * \tilde E8) 3 pour le sommet 4 (i.e. cas euclidien: \tilde E7)
   */

  bool isSpherical(pathSize <= 7 ? true : false);
  unsigned int iStart(isSpherical || pathSize == 8 ? 2 : 3);

  // E6, E7, E8, \tilde E8
  AnToEn_AnToTEn(pathTemp, bVerticesLinkable, isSpherical, iStart);

  if (pathSize == 7) // \tile E7
    AnToEn_AnToTEn(pathTemp, bVerticesLinkable, false, 3);
}

void CoxIter::AnToEn_AnToTEn(const vector<short unsigned int> &pathTemp,
                             const vector<bool> &bVerticesLinkable,
                             const bool &isSpherical,
                             const short unsigned int &start) {
  unsigned int pathSize(pathTemp.size()), j, k, l;
  vector<bool> bVerticesLinkableTemp, bVerticesLinkableTempTemp;

  /*
   * 	Ici, on a donc un An (n=5, 6 ou 7) avec 1 -- 2 -- 3 -- 4 -- 5 ...
   * 	On va cherche si iStart a un voisin admissible
   */
  for (unsigned int i(0); i < verticesCount; i++) {
    // si le sommet est pas utilisbale (s'il l'est c'est qu'il n'est pas voisin
    // de la base) ET si y'a un lien
    if (false == bVerticesLinkable[i] &&
        coxeterMatrix[i][pathTemp[start]] == 3) {
      // on va chercher si c'est uniquement à cause d'un des sommets différents
      // de iStart que le sommet n'est pas admissible
      for (j = 0; j < pathSize; j++) {
        if (j != start && coxeterMatrix[pathTemp[j]][i] != 2)
          break;
      }

      if (j == pathSize) // admissible
      {
        bVerticesLinkableTemp = bVerticesLinkable;
        for (k = 0; k < verticesCount; k++) {
          if (coxeterMatrix[i][k] != 2)
            bVerticesLinkableTemp[k] = false;
        }

        if (isSpherical) {
          graphsList_spherical->addGraph(pathTemp, bVerticesLinkableTemp, 4,
                                         true, i); // En

          // on a un E6 qu'on va tenter de prolonger en un TE6
          if (pathSize == 5) {
            for (j = 0; j < verticesCount; j++) {
              if (coxeterMatrix[i][j] == 3 && bVerticesLinkable[j]) {
                bVerticesLinkableTempTemp = bVerticesLinkableTemp;
                for (l = 0; l < verticesCount; l++) {
                  if (coxeterMatrix[j][l] != 2)
                    bVerticesLinkableTempTemp[l] = false;
                }

                graphsList_euclidean->addGraph(
                    pathTemp, bVerticesLinkableTempTemp, 4, false, i, j); // En
              }
            }
          }
        } else
          graphsList_euclidean->addGraph(pathTemp, bVerticesLinkableTemp, 4,
                                         false, i); // TEn
      }
    }
  }
}

/*
 * 	B3ToF4_B4ToTF4
 * 		On tente de prolonger un B3 en un F4 ou un B4 en \tilde F4
 * (euclidien)
 * 										                  4
 * 			Le B3 donné est: pathTemp[0] --------- pathTemp[1]
 * ------------ endVertex
 *
 * 			Le B4 donné est: pathTemp[0] --------- pathTemp[1]
 * ------------ pathTemp[2] ------------ endVertex
 * 													            4
 *
 * 		Paramètres:
 * 			bVerticesBeginLinkable: ce qui est linkable à cause du
 * premier (deux premiers) sommet(s) pathTemp: deux(trois) premiers sommets
 * endVertex: 3ème(4ème) sommet
 * */
void CoxIter::B3ToF4_B4ToTF4(const vector<bool> &bVerticesBeginLinkable,
                             vector<short unsigned int> pathTemp,
                             const short unsigned int &endVertex) {
  bool isSpherical(pathTemp.size() == 2); // true si sphérique (on cherche F4),
                                          // false si euclidien (on cherche TF4)
  unsigned int i, j, iV2(pathTemp[1]);
  vector<bool> bVerticesLinkable;

  pathTemp.push_back(endVertex);

  // on va parcourir les sommets et regarder les voisins de poids 3 de endVertex
  for (i = 0; i < verticesCount; i++) {
    if (coxeterMatrix[endVertex][i] == 3 && bVerticesBeginLinkable[i] &&
        coxeterMatrix[i][iV2] == 2) {
      bVerticesLinkable = bVerticesBeginLinkable;

      if (!isSpherical &&
          coxeterMatrix[pathTemp[2]][i] !=
              2) // If i is connected to pathTemp[2], this won't work
        continue;

      for (j = 0; j < verticesCount; j++) {
        if (coxeterMatrix[j][iV2] != 2 || coxeterMatrix[j][endVertex] != 2 ||
            coxeterMatrix[j][i] != 2)
          bVerticesLinkable[j] = false;
      }

      if (isSpherical)
        graphsList_spherical->addGraph(pathTemp, bVerticesLinkable, 5, true,
                                       i); // Fn
      else
        graphsList_euclidean->addGraph(pathTemp, bVerticesLinkable, 5, false,
                                       i); // TFn
    }
  }
}

void CoxIter::printPath() {
  if (path.size() == 1)
    return;

  for (const auto &p : path)
    cout << p << " ; ";

  cout << endl;
}

/*! 	\fn vector2str
 * 	\brief Vector --> string
 * 	\param polynomial(const vector< mpz_class >& polynomial) Polynomial
 *
 * 	\remark: When mpz_class will have a to_string(mpz) we'll move that to
 * tools/polynomials.h TODO
 */
string vector2str(const vector<mpz_class> &polynomial) {
  bool bFirst(true);
  unsigned int iSize(polynomial.size());
  string strRes;
  mpz_class mpzTemp;

  for (unsigned int i(0); i < iSize; i++) {
    if (polynomial[i] != 0) {
      if (bFirst) {
        strRes += polynomial[i].get_str() +
                  (i ? " * x" + string(i > 1 ? "^" + to_string(i) : "") : "");
        bFirst = false;
      } else {
        if ((polynomial[i] != 1 && polynomial[i] != -1) || !i) {
          mpzTemp = abs(polynomial[i]);
          strRes += (polynomial[i] > 0 ? " + " : " - ") + mpzTemp.get_str() +
                    (i ? " * x" + string(i > 1 ? "^" + to_string(i) : "") : "");
        } else
          strRes += (polynomial[i] > 0 ? " + " : " - ") +
                    string(i > 1 ? "x^" + to_string(i) : "x");
      }
    }
  }

  return strRes;
}

const vector<vector<GraphsProductSet>> *
CoxIter::get_ptr_graphsProducts() const {
  return &graphsProducts;
}

string CoxIter::get_strGrowthSeries_raw() {
  string strGrowth;

  if (strOuputMathematicalFormat == "generic") {
    strGrowth = "g(x) = (" + growthSeries_raw + ")^-1;";
  } else if (strOuputMathematicalFormat == "mathematica") {
    strGrowth = "Symb[s_, x_] := Product[Sum[x^i, {i, 0, s[[i]] - 1}], {i, 1, "
                "Length[s]}];\n";
    strGrowth += "g[x_] := (" + growthSeries_raw + ")^-1;";
  } else if (strOuputMathematicalFormat == "pari") {
    strGrowth = "Symb = (S, y) -> prod(i=1, length(S), sum(i=0,S[i]-1,y^i));\n";
    strGrowth += "g(x) = (" + growthSeries_raw + ")^-1;";
  }

  return strGrowth;
}

string CoxIter::get_strGrowthSeries() {
  string strGrowth;

  if (strOuputMathematicalFormat == "generic") {
    strGrowth = "f(x) = ";
    if (growthSeries_cyclotomicNumerator.size()) {
      strGrowth += "C(";
      unsigned int max(growthSeries_cyclotomicNumerator.size());
      for (unsigned int i(0); i < max; i++)
        strGrowth +=
            (i ? "," : "") + to_string(growthSeries_cyclotomicNumerator[i]);

      strGrowth += ")";
    }

    strGrowth += "/(" + vector2str(growthSeries_polynomialDenominator) + ")";
  } else if (strOuputMathematicalFormat == "gap") {
    unsigned int cycloSize(growthSeries_cyclotomicNumerator.size());
    unsigned int denominatorSize(growthSeries_polynomialDenominator.size());

    strGrowth += "f := Product([";
    for (unsigned int i(0); i < cycloSize; i++)
      strGrowth +=
          (i ? "," : "") + to_string(growthSeries_cyclotomicNumerator[i]);
    strGrowth += "], i -> CyclotomicPolynomial(Rationals,i))/ValuePol([";
    for (unsigned int i(0); i < denominatorSize; i++)
      cout << (i ? "," : "") << growthSeries_polynomialDenominator[i];

    cout << "], X(Rationals));";
  } else if (strOuputMathematicalFormat == "mathematica") {
    unsigned int cycloSize(growthSeries_cyclotomicNumerator.size());

    strGrowth =
        "Cyclo[s_, x_] := Product[Cyclotomic[s[[i]], x], {i, 1, Length[s]}];";

    strGrowth += "f[x_] := Cyclo[{";
    for (unsigned int i(0); i < cycloSize; i++)
      strGrowth +=
          (i ? "," : "") + to_string(growthSeries_cyclotomicNumerator[i]);

    strGrowth +=
        "},x]/(" + vector2str(growthSeries_polynomialDenominator) + ");";
  } else if (strOuputMathematicalFormat == "pari") {
    unsigned int cycloSize(growthSeries_cyclotomicNumerator.size());
    strGrowth = "Cyclo = (S, y) -> prod(i=1, length(S), polcyclo(S[i],y));\n";
    strGrowth += "f(x) = Cyclo([";
    for (unsigned int i(0); i < cycloSize; i++)
      strGrowth +=
          (i ? "," : "") + to_string(growthSeries_cyclotomicNumerator[i]);

    strGrowth +=
        "],x)/(" + vector2str(growthSeries_polynomialDenominator) + ");";
  }

  return strGrowth;
}

void CoxIter::printGrowthSeries() {
  if (!isGrowthSeriesComputed)
    growthSeries();

  if (strOuputMathematicalFormat == "generic") {
    cout << "f(x) = ";
    if (growthSeries_cyclotomicNumerator.size()) {
      cout << "C(";
      unsigned int max(growthSeries_cyclotomicNumerator.size());
      for (unsigned int i(0); i < max; i++)
        cout << (i ? "," : "") << growthSeries_cyclotomicNumerator[i];

      cout << ")";
    }

    cout << "/(";
    Polynomials::polynomialDisplay(growthSeries_polynomialDenominator);
    cout << ")";

    if (debug)
      cout << "\ng(x) = (" << growthSeries_raw << ")^-1;";
  } else if (strOuputMathematicalFormat == "gap") {
    unsigned int cycloSize(growthSeries_cyclotomicNumerator.size());
    unsigned int denominatorSize(growthSeries_polynomialDenominator.size());

    cout << "f := Product([";
    for (unsigned int i(0); i < cycloSize; i++)
      cout << (i ? "," : "") << growthSeries_cyclotomicNumerator[i];

    cout << "], i -> CyclotomicPolynomial(Rationals,i))/ValuePol([";
    for (unsigned int i(0); i < denominatorSize; i++)
      cout << (i ? "," : "") << growthSeries_polynomialDenominator[i];

    cout << "], X(Rationals));";

    if (debug)
      cout << "\ng(x) = (" << growthSeries_raw << ")^-1;";
  } else if (strOuputMathematicalFormat == "mathematica") {
    unsigned int cycloSize(growthSeries_cyclotomicNumerator.size());

    cout
        << "Cyclo[s_, x_] := Product[Cyclotomic[s[[i]], x], {i, 1, Length[s]}];"
        << endl;
    if (debug)
      cout << "Symb[s_, x_] := Product[Sum[x^i, {i, 0, s[[i]] - 1}], {i, 1, "
              "Length[s]}];"
           << endl;

    cout << "f[x_] := Cyclo[{";
    for (unsigned int i(0); i < cycloSize; i++)
      cout << (i ? "," : "") << growthSeries_cyclotomicNumerator[i];
    cout << "},x]";

    cout << "/(";
    Polynomials::polynomialDisplay(growthSeries_polynomialDenominator);
    cout << ");";

    if (debug)
      cout << "\ng[x_] := (" << growthSeries_raw << ")^-1;";
  } else if (strOuputMathematicalFormat == "pari") {
    unsigned int cycloSize(growthSeries_cyclotomicNumerator.size());

    cout << "Cyclo = (S, y) -> prod(i=1, length(S), polcyclo(S[i],y));" << endl;
    if (debug)
      cout << "Symb = (S, y) -> prod(i=1, length(S), sum(i=0,S[i]-1,y^i));"
           << endl;

    cout << "f(x) = Cyclo([";
    for (unsigned int i(0); i < cycloSize; i++)
      cout << (i ? "," : "") << growthSeries_cyclotomicNumerator[i];
    cout << "],x)/(";
    Polynomials::polynomialDisplay(growthSeries_polynomialDenominator);
    cout << ");";

    if (debug)
      cout << "\ng(x) = (" << growthSeries_raw << ")^-1;";
  }
}

ostream &operator<<(ostream &o, const CoxIter &g) {
  o << "Graphes sphériques: " << endl;
  o << *g.graphsList_spherical;

  o << "Graphes euclidien: " << endl;
  o << *g.graphsList_euclidean;

  return o;
}

int CoxIter::isGraphCocompact() {
  if (isCocompact >= 0)
    return isCocompact;

  if (!isGraphsProductsComputed)
    computeGraphsProducts();

  if (!checkCocompactness || !graphsProducts.size()) {
    isCocompact = -1;
    return -1;
  }

  if (!graphsProducts[1].size() ||
      !graphsProducts[0].size()) // No vertices, no edges
  {
    isCocompact = 0;
    return 0;
  }

  if (hasBoldLine) {
    isCocompact = 0;
    return 0;
  }

  // ----------------------------------------------------
  // the test
  if (bUseOpenMP && verticesCount >= 15)
    isCocompact = isGraph_cocompact_finiteVolume_parallel(1) ? 1 : 0;
  else
    isCocompact = isGraph_cocompact_finiteVolume_sequential(1) ? 1 : 0;

  return isCocompact;
}

int CoxIter::checkCovolumeFiniteness() {
  if (isFiniteCovolume >= 0)
    return isFiniteCovolume;

  if (!isGraphsProductsComputed)
    computeGraphsProducts();

  // ----------------------------------------------------
  // some stupid tests
  if (!graphsProducts.size()) {
    isFiniteCovolume = -1;
    return -1;
  }

  if (!graphsProducts[2].size() && !graphsProducts[1].size()) // No vertices
  {
    isFiniteCovolume = 0;
    return 0;
  }

  if (!graphsProducts[0].size()) // No edges
  {
    isFiniteCovolume = 0;
    return 0;
  }

  /*
   * The duplication of the data is not beautiful but it allows
   * us to simplify the code of the function isGraph_cocompact_finiteVolume.
   * Moreover, we want the program to be quick but memory is not really an
   * issue.
   *
   * [2] Contains a combined list of finite and infinite vertices
   */
  graphsProducts[2].insert(graphsProducts[2].end(), graphsProducts[1].begin(),
                           graphsProducts[1].end());

  // ----------------------------------------------------
  // the test
  if (bUseOpenMP && verticesCount >= 15)
    isFiniteCovolume = isGraph_cocompact_finiteVolume_parallel(2) ? 1 : 0;
  else
    isFiniteCovolume = isGraph_cocompact_finiteVolume_sequential(2) ? 1 : 0;

  return isFiniteCovolume;
}

bool CoxIter::isGraph_cocompact_finiteVolume_sequential(unsigned int index) {
  unsigned int extendedCount;

  vector<Graph *> diffSubNotBig, diffBigNotSub;
  vector<Graph *>::const_iterator itGBig;

  bool isExtendable;

  for (const auto &sphericalProductsCodim1 : graphsProducts[0]) {
    extendedCount = 0;

    for (const auto &gpBig : graphsProducts[index]) {
      diffSubNotBig.clear();
      diffBigNotSub.clear();

      set_difference(sphericalProductsCodim1.graphs.begin(),
                     sphericalProductsCodim1.graphs.end(), gpBig.graphs.begin(),
                     gpBig.graphs.end(), back_inserter(diffSubNotBig));

      set_difference(gpBig.graphs.begin(), gpBig.graphs.end(),
                     sphericalProductsCodim1.graphs.begin(),
                     sphericalProductsCodim1.graphs.end(),
                     back_inserter(diffBigNotSub));

      isExtendable = true;
      for (const auto &graphSub : diffSubNotBig) {
        for (itGBig = diffBigNotSub.begin(); itGBig != diffBigNotSub.end();
             ++itGBig) {
          if (graphSub->bIsSubgraphOf(*itGBig))
            break;
        }

        if (itGBig == diffBigNotSub.end()) {
          isExtendable = false;
          break;
        }
      }

      if (isExtendable)
        extendedCount++;
    }

    if (extendedCount != 2) {
      if (debug) {
        cout << "----------------------------------------------------------"
             << endl;
        cout << (index == 1 ? "Compactness" : "Finite covolume") << " test"
             << endl;
        cout << "Trying to extend the product: " << endl;
        cout << sphericalProductsCodim1 << endl;
        cout << "Succeeded in " << extendedCount << " ways instead of 2"
             << endl;

        for (vector<GraphsProductSet>::const_iterator gpBig(
                 graphsProducts[index].begin());
             gpBig != graphsProducts[index].end(); ++gpBig) {
          if (sphericalProductsCodim1.areVerticesSubsetOf(*gpBig))
            cout << "Candidate: \n" << *gpBig << endl;
        }
        cout << "----------------------------------------------------------"
             << endl;
      }

      return 0;
    }
  }

  return 1;
}

bool CoxIter::isGraph_cocompact_finiteVolume_parallel(unsigned int index) {
  unsigned int extendedCount, max(graphsProducts[0].size()), i;

  vector<Graph *> diffSubNotBig, diffBigNotSub;
  vector<Graph *>::const_iterator itGBig;

  bool isExtendable, exit(false);

#pragma omp parallel if (bUseOpenMP && verticesCount >= 15)
  {
#pragma omp single nowait
    {
      for (i = 0; i < max && !exit; i++) {
#pragma omp task private(itGBig, isExtendable, diffSubNotBig, diffBigNotSub,   \
                         extendedCount) shared(max, index, exit)               \
    firstprivate(i)
        {
          extendedCount = 0;

          for (const auto &gpBig : graphsProducts[index]) {
            diffSubNotBig.clear();
            diffBigNotSub.clear();
            set_difference(graphsProducts[0][i].graphs.begin(),
                           graphsProducts[0][i].graphs.end(),
                           gpBig.graphs.begin(), gpBig.graphs.end(),
                           back_inserter(diffSubNotBig));

            set_difference(gpBig.graphs.begin(), gpBig.graphs.end(),
                           graphsProducts[0][i].graphs.begin(),
                           graphsProducts[0][i].graphs.end(),
                           back_inserter(diffBigNotSub));

            isExtendable = true;
            for (const auto &graphSub : diffSubNotBig) {
              for (itGBig = diffBigNotSub.begin();
                   itGBig != diffBigNotSub.end(); ++itGBig) {
                if (graphSub->bIsSubgraphOf(*itGBig))
                  break;
              }

              if (itGBig == diffBigNotSub.end()) {
                isExtendable = false;
                break;
              }
            }

            if (isExtendable)
              extendedCount++;
          }

          if (extendedCount != 2) {
            if (debug) {
#pragma omp critical
              {
                cout << "------------------------------------------------------"
                        "----"
                     << endl;
                cout << (index == 1 ? "Compactness" : "Finite covolume")
                     << " test" << endl;
                cout << "Trying to extend the product: " << endl;
                cout << graphsProducts[0][i] << endl;
                cout << "Succeeded in " << extendedCount << " ways instead of 2"
                     << endl;

                for (const auto &gpBig : graphsProducts[index]) {
                  if (graphsProducts[0][i].areVerticesSubsetOf(gpBig))
                    cout << "Candidate: \n" << gpBig << endl;
                }
                cout << "------------------------------------------------------"
                        "----"
                     << endl;
              }
            }

#pragma omp atomic write
            exit = true;
          }
        }

        if (exit)
          break;
      }
    }
  }

  return !exit;
}

void CoxIter::computeGraphsProducts() {
  if (isGraphsProductsComputed)
    return;

  if (!isGraphExplored)
    exploreGraph();

  if (debug) {
    cout << "Connected spherical graphs" << endl;
    cout << *this->graphsList_spherical << endl;
  }

  graphsProducts = vector<vector<GraphsProductSet>>(3);
  vector<bool> bGPVerticesNonLinkable(vector<bool>(verticesCount, false));
  GraphsProduct gp; ///< Current graphs product

  // --------------------------------------------------------------
  // produits de graphes sphériques
  GraphsListIterator grIt_spherical(this->graphsList_spherical);

#pragma omp parallel if (bUseOpenMP && verticesCount >= 15)
  {
#pragma omp single nowait
    while (grIt_spherical.ptr) {
#pragma omp task firstprivate(grIt_spherical, bGPVerticesNonLinkable, gp)
      {
        computeGraphsProducts(grIt_spherical, &graphsProductsCount_spherical,
                              true, gp, bGPVerticesNonLinkable);
      }

      ++grIt_spherical;
    }
  }

  // --------------------------------------------------------------
  // produits de graphes euclidiens
  if (debug) {
    cout << "Connected euclidean graphs" << endl;
    cout << *this->graphsList_euclidean;
  }

  GraphsListIterator grIt_euclidean(this->graphsList_euclidean);
#pragma omp parallel if (bUseOpenMP && verticesCount >= 15)
  {
#pragma omp single nowait
    while (grIt_euclidean.ptr) {
#pragma omp task firstprivate(grIt_euclidean, bGPVerticesNonLinkable, gp)
      {
        computeGraphsProducts(grIt_euclidean, &graphsProductsCount_euclidean,
                              false, gp, bGPVerticesNonLinkable);
      }

      ++grIt_euclidean;
    }
  }

  if (debug) {
    cout << "\nProduct of euclidean graphs" << endl;
    printEuclideanGraphsProducts(&graphsProductsCount_euclidean);
  }

  isGraphsProductsComputed = true;

  // ---------------------------------------------------------
  // We guess the dimension
  if (!dimension) {
    dimension = max(euclideanMaxRankFound + 1, sphericalMaxRankFound);
    isDimensionGuessed = true;

    if (euclideanMaxRankFound == sphericalMaxRankFound) {
      graphsProducts[0] = graphsProducts[1];
      graphsProducts[1].clear();
    } else if (euclideanMaxRankFound > sphericalMaxRankFound) {
      graphsProducts[0].clear();
      graphsProducts[1].clear();
    } else if (sphericalMaxRankFound > euclideanMaxRankFound + 1) {
      graphsProducts[2].clear();
    }
  }
}

void CoxIter::computeGraphsProducts(
    GraphsListIterator grIt,
    vector<map<vector<vector<short unsigned int>>, unsigned int>>
        *graphsProductsCount,
    const bool &isSpherical, GraphsProduct &gp,
    vector<bool> &bGPVerticesNonLinkable) {
  vector<short unsigned int>::iterator iIt;
  vector<short unsigned int> flaggedVertices;
  unsigned int graphRank(0);
  vector<vector<short unsigned int>> vFootPrintTest;

  while (grIt.ptr && (gp.rank + graphRank <= maximalSubgraphRank)) {
    // ---------------------------------------------------
    // est ce que le graphe est admissible?
    for (iIt = grIt.ptr->vertices.begin(); iIt != grIt.ptr->vertices.end();
         ++iIt) {
      if (bGPVerticesNonLinkable[*iIt]) // si pas linkable
        break;
    }

    // si le graphe est admissible
    if (iIt == grIt.ptr->vertices.end()) {
      // le graphe est ajouté au produit
      gp.graphs.push_back(grIt.ptr);

      // taille du graphe courant
      graphRank = isSpherical ? grIt.ptr->vertices.size()
                              : (grIt.ptr->vertices.size() - 1);
      gp.rank += graphRank;

      // Create the footprint of the product. The goal is to decide if we
      // already have this product
      vFootPrintTest = gp.createFootPrint();

#pragma omp critical
      {
        if (checkCocompactness || checkCofiniteness) {
          if (dimension) // If we know the dimension, everything is easier
          {
            if (isSpherical) {
              // Keeping track of spherical subgraphs
              if ((gp.rank == (dimension - 1) || gp.rank == dimension))
                graphsProducts[gp.rank + 1 - dimension].push_back(
                    GraphsProductSet(gp));
            }

            // Euclidean subgraphs
            if (!isSpherical && gp.rank == (dimension - 1) && checkCofiniteness)
              graphsProducts[2].push_back(GraphsProductSet(gp));
          } else {
            if (isSpherical) {
              if (gp.rank == sphericalMaxRankFound + 1) {
                graphsProducts[0] = graphsProducts[1];
                graphsProducts[1].clear();
                graphsProducts[1].push_back(GraphsProductSet(gp));
              } else if (gp.rank > sphericalMaxRankFound + 1) {
                graphsProducts[0].clear();
                graphsProducts[1].clear();
                graphsProducts[1].push_back(GraphsProductSet(gp));
              } else if (gp.rank + 1 >= sphericalMaxRankFound)
                graphsProducts[gp.rank + 1 - sphericalMaxRankFound].push_back(
                    GraphsProductSet(gp));
            } else {
              if (checkCofiniteness) {
                if (gp.rank > euclideanMaxRankFound)
                  graphsProducts[2].clear();

                if (gp.rank >= euclideanMaxRankFound)
                  graphsProducts[2].push_back(GraphsProductSet(gp));
              }
            }
          }
        }

        if (isSpherical && gp.rank >= sphericalMaxRankFound)
          sphericalMaxRankFound = gp.rank;

        if (!isSpherical && gp.rank >= euclideanMaxRankFound)
          euclideanMaxRankFound = gp.rank;

        if ((*graphsProductsCount)[gp.rank].find(vFootPrintTest) ==
            (*graphsProductsCount)[gp.rank].end())
          (*graphsProductsCount)[gp.rank][vFootPrintTest] = 1;
        else
          (*graphsProductsCount)[gp.rank][vFootPrintTest]++;
      }

      // mise à jour des sommets que l'on ne peut plus prendre
      for (unsigned int i = 0; i < verticesCount; i++) {
        if (!grIt.ptr->bVerticesLinkable[i] && !bGPVerticesNonLinkable[i]) {
          flaggedVertices.push_back(i);
          bGPVerticesNonLinkable[i] = true;
        }
      }

      // récursion
      computeGraphsProducts(++grIt, graphsProductsCount, isSpherical, gp,
                            bGPVerticesNonLinkable);

      // -----------------------------------------------
      // dé-initialisations

      // on remet la liste à son état d'avant la récursion
      for (iIt = flaggedVertices.begin(); iIt != flaggedVertices.end(); ++iIt)
        bGPVerticesNonLinkable[*iIt] = false;

      gp.rank -= graphRank;

      // le graphe est enlevé
      gp.graphs.pop_back();

      if (!gp.graphs.size())
        break;
    } else
      ++grIt;
  }
}

void CoxIter::IS_computations(const string &t0, const string &s0) {
  infSeq_t0 = get_vertexIndex(t0);
  infSeq_s0 = get_vertexIndex(s0);

  infSeqFVectorsUnits = vector<unsigned int>(dimension, 0);
  infSeqFVectorsPowers = vector<unsigned int>(dimension, 0);

  if (!coxeterMatrix.size())
    return;

  if (!isGraphExplored)
    exploreGraph();

  if (!isGraphsProductsComputed)
    computeGraphsProducts();

  vector<bool> bGPVerticesNonLinkable(vector<bool>(verticesCount, false));
  GraphsProduct gp; ///< Current graphs product

  // --------------------------------------------------------------
  // produits de graphes sphériques
  GraphsListIterator grIt_spherical(this->graphsList_spherical);

#pragma omp parallel if (bUseOpenMP && verticesCount >= 15)
  {
#pragma omp single nowait
    while (grIt_spherical.ptr) {
#pragma omp task firstprivate(grIt_spherical, bGPVerticesNonLinkable, gp)
      {
        computeGraphsProducts_IS(grIt_spherical, true, gp,
                                 bGPVerticesNonLinkable);
      }

      ++grIt_spherical;
    }
  }

  // --------------------------------------------------------------
  // produits de graphes euclidiens
  GraphsListIterator grIt_euclidean(this->graphsList_euclidean);
#pragma omp parallel if (bUseOpenMP && verticesCount >= 15)
  {
#pragma omp single nowait
    while (grIt_euclidean.ptr) {
#pragma omp task firstprivate(grIt_euclidean, bGPVerticesNonLinkable, gp)
      {
        computeGraphsProducts_IS(grIt_euclidean, false, gp,
                                 bGPVerticesNonLinkable);
      }

      ++grIt_euclidean;
    }
  }
}

void CoxIter::computeGraphsProducts_IS(GraphsListIterator grIt,
                                       const bool &isSpherical,
                                       GraphsProduct &gp,
                                       vector<bool> &bGPVerticesNonLinkable) {
  vector<short unsigned int>::iterator iIt;
  vector<short unsigned int> flaggedVertices;
  unsigned int graphRank(0);

  while (grIt.ptr && (gp.rank + graphRank <= maximalSubgraphRank)) {
    // ---------------------------------------------------
    // est ce que le graphe est admissible?
    for (iIt = grIt.ptr->vertices.begin(); iIt != grIt.ptr->vertices.end();
         ++iIt) {
      if (bGPVerticesNonLinkable[*iIt]) // si pas linkable
        break;
    }

    // si le graphe est admissible
    if (iIt == grIt.ptr->vertices.end()) {
      // le graphe est ajouté au produit
      gp.graphs.push_back(grIt.ptr);

      // taille du graphe courant
      graphRank = isSpherical ? grIt.ptr->vertices.size()
                              : (grIt.ptr->vertices.size() - 1);
      gp.rank += graphRank;

#pragma omp critical
      {
        if (isSpherical || gp.rank == (dimension - 1)) {
          bool bSpecialIn_t0(false), bSpecialIn_s0(false);
          unsigned int nonCommute_t0(0), nonCommute_s0(0);

          for (const auto &g : gp.graphs) {
            for (const auto &v : g->vertices) {
              if (v == infSeq_t0)
                bSpecialIn_t0 = true;
              else if (coxeterMatrix[v][infSeq_t0] != 2)
                nonCommute_t0++;

              if (v == infSeq_s0)
                bSpecialIn_s0 = true;
              else if (coxeterMatrix[v][infSeq_s0] != 2)
                nonCommute_s0++;
            }
          }

          if (!bSpecialIn_t0 && !bSpecialIn_s0) {
            if (!nonCommute_t0 && !nonCommute_s0)
              infSeqFVectorsUnits[isSpherical ? dimension - gp.rank : 0]++;
            else if (nonCommute_t0 && nonCommute_s0)
              infSeqFVectorsPowers[isSpherical ? dimension - gp.rank : 0] += 2;
            else if (nonCommute_t0 && !nonCommute_s0) {
              infSeqFVectorsUnits[isSpherical ? dimension - gp.rank : 0]++;
              infSeqFVectorsPowers[isSpherical ? dimension - gp.rank : 0]++;
            } else
              infSeqFVectorsPowers[isSpherical ? dimension - gp.rank : 0]++;
          } else if (!bSpecialIn_t0 && bSpecialIn_s0) {
            if (!nonCommute_s0)
              infSeqFVectorsUnits[isSpherical ? dimension - gp.rank : 0] += 2;
            else {
              infSeqFVectorsUnits[isSpherical ? dimension - gp.rank : 0]++;
              infSeqFVectorsPowers[isSpherical ? dimension - gp.rank : 0]++;
            }
          } else if (bSpecialIn_t0 && !bSpecialIn_s0) {
            if (nonCommute_t0)
              infSeqFVectorsPowers[isSpherical ? dimension - gp.rank : 0]++;
          } else
            infSeqFVectorsUnits[isSpherical ? dimension - gp.rank : 0]++;
        }
      }

      // mise à jour des sommets que l'on ne peut plus prendre
      for (unsigned int i = 0; i < verticesCount; i++) {
        if (!grIt.ptr->bVerticesLinkable[i] && !bGPVerticesNonLinkable[i]) {
          flaggedVertices.push_back(i);
          bGPVerticesNonLinkable[i] = true;
        }
      }

      // récursion
      computeGraphsProducts_IS(++grIt, isSpherical, gp, bGPVerticesNonLinkable);

      // -----------------------------------------------
      // dé-initialisations

      // on remet la liste à son état d'avant la récursion
      for (iIt = flaggedVertices.begin(); iIt != flaggedVertices.end(); ++iIt)
        bGPVerticesNonLinkable[*iIt] = false;

      gp.rank -= graphRank;

      // le graphe est enlevé
      gp.graphs.pop_back();

      if (!gp.graphs.size())
        break;
    } else
      ++grIt;
  }
}

bool CoxIter::bCanBeFiniteCovolume() {
  // -----------------------------------------------------------
  // Some verifications
  if (!dimension)
    throw(string("CoxIter::bCanBeFiniteCovolume: Dimension not specified"));

  if (!isGraphExplored)
    exploreGraph();

  // -----------------------------------------------------------
  // Initializations
  graphsProducts_bCanBeFiniteCovolume = vector<vector<GraphsProductSet>>(1);

  GraphsListIterator grIt_euclidean(this->graphsList_euclidean);
  vector<bool> bGPVerticesNonLinkable(vector<bool>(verticesCount, false));
  GraphsProduct gp; ///< Current graphs product

// -----------------------------------------------------------
// We find the products of euclidean graphs
#pragma omp parallel if (bUseOpenMP && verticesCount >= 15)
  {
#pragma omp single nowait
    while (grIt_euclidean.ptr) {
#pragma omp task firstprivate(grIt_euclidean, bGPVerticesNonLinkable, gp)
      {
        bCanBeFiniteCovolume_computeGraphsProducts(grIt_euclidean, gp,
                                                   bGPVerticesNonLinkable);
      }

      ++grIt_euclidean;
    }
  }

  // -----------------------------------------------------------
  // Check the condition: every connected affine graph of rank at least 2 is a
  // subgraph of an affine graph of rank n-1
  grIt_euclidean = GraphsListIterator(
      GraphsListIterator(this->graphsList_euclidean, 3)); // TODO: partir à 2?
  bool bCanBeExtended;

  while (grIt_euclidean.ptr) {
    bCanBeExtended = false;
    for (auto graphProd :
         graphsProducts_bCanBeFiniteCovolume[0]) // TODO OPTIMIZATION
                                                 // paralleliser?
    {
      for (auto gr : graphProd.graphs) {
        // For two affine graphs G1 and G2, G1 is a subgraph of G2 iff G1=G2
        if (*grIt_euclidean.ptr == *gr) {
          bCanBeExtended = true;
          break;
        }
      }

      if (bCanBeExtended)
        break;
    }

    if (!bCanBeExtended) {
      if (debug) {
        cout << "Can be of finite covolume: no" << endl;
        cout << "\tCannot extend the affine graph: " << endl;
        cout << "\t" << *grIt_euclidean.ptr << "\n" << endl;
      }

      return false;
    }

    ++grIt_euclidean;
  }

  return true;
}

void CoxIter::bCanBeFiniteCovolume_computeGraphsProducts(
    GraphsListIterator grIt, GraphsProduct &gp,
    vector<bool> &bGPVerticesNonLinkable) {
  vector<short unsigned int>::iterator iIt;
  vector<short unsigned int> flaggedVertices;
  unsigned int graphGrank(0);

  while (grIt.ptr && (gp.rank + graphGrank <= verticesCount)) {
    // ---------------------------------------------------
    // est ce que le graphe est admissible?
    for (iIt = grIt.ptr->vertices.begin(); iIt != grIt.ptr->vertices.end();
         ++iIt) {
      if (bGPVerticesNonLinkable[*iIt]) // si pas linkable
        break;
    }

    // si le graphe est admissible
    if (iIt == grIt.ptr->vertices.end()) {
      // le graphe est ajouté au produit
      gp.graphs.push_back(grIt.ptr);

      // taille du graphe courant
      graphGrank = grIt.ptr->vertices.size() - 1;
      gp.rank += graphGrank;

#pragma omp critical
      {
        if (gp.rank == (dimension - 1))
          graphsProducts_bCanBeFiniteCovolume[0].push_back(
              GraphsProductSet(gp));
      }

      // mise à jour des sommets que l'on ne peut plus prendre
      for (unsigned int i = 0; i < verticesCount; i++) {
        if (!grIt.ptr->bVerticesLinkable[i] && !bGPVerticesNonLinkable[i]) {
          flaggedVertices.push_back(i);
          bGPVerticesNonLinkable[i] = true;
        }
      }

      // récursion
      bCanBeFiniteCovolume_computeGraphsProducts(++grIt, gp,
                                                 bGPVerticesNonLinkable);

      // -----------------------------------------------
      // dé-initialisations

      // on remet la liste à son état d'avant la récursion
      for (iIt = flaggedVertices.begin(); iIt != flaggedVertices.end(); ++iIt)
        bGPVerticesNonLinkable[*iIt] = false;

      gp.rank -= graphGrank;

      // le graphe est enlevé
      gp.graphs.pop_back();

      if (!gp.graphs.size())
        break;
    } else
      ++grIt;
  }
}

vector<vector<short unsigned int>> CoxIter::bCanBeFiniteCovolume_complete() {
  // -----------------------------------------------------------
  // Some verifications
  if (!dimension)
    throw(string(
        "CoxIter::bCanBeFiniteCovolume_complete: Dimension not specified"));

  if (!isGraphExplored)
    exploreGraph();

  // -----------------------------------------------------------
  // Initializations
  graphsProducts_bCanBeFiniteCovolume =
      vector<vector<GraphsProductSet>>(dimension + 1);

  GraphsListIterator grIt_euclidean(this->graphsList_euclidean);
  vector<bool> bGPVerticesNonLinkable(vector<bool>(verticesCount, false));
  GraphsProduct gp; ///< Current graphs product

// -----------------------------------------------------------
// We find the products of euclidean graphs
#pragma omp parallel if (bUseOpenMP && verticesCount >= 15)
  {
#pragma omp single nowait
    while (grIt_euclidean.ptr) {

#pragma omp task firstprivate(grIt_euclidean, bGPVerticesNonLinkable, gp)
      {
        bCanBeFiniteCovolume_complete_computeGraphsProducts(
            grIt_euclidean, gp, bGPVerticesNonLinkable);
      }

      ++grIt_euclidean;
    }
  }

  vector<vector<short unsigned int>> graphsNotExtendable(0);

  // -----------------------------------------------------------
  // Check the condition: every connected affine graph of rank at least 2 is a
  // subgraph of an affine graph of rank n-1
  bool bCanBeExtended, bSubproduct, bFound;

  for (unsigned int i(2); i < dimension - 1; i++) {
    for (auto gpSmall : graphsProducts_bCanBeFiniteCovolume[i]) {
      // If gpSmall is not a subgraph of gpBig for every gpBig (affine of rank
      // dimension-1), then the graph is not of finite covolume
      bCanBeExtended = false;

      for (auto gpBig : graphsProducts_bCanBeFiniteCovolume[dimension - 1]) {
        // First test
        if (!gpSmall.areVerticesSubsetOf(gpBig))
          continue;

        bSubproduct = true;
        for (auto gSmall : gpSmall.graphs) {
          bFound = false;
          for (auto gBig : gpBig.graphs) {
            if (*gSmall == *gBig) {
              bFound = true;
              break;
            }
          }

          if (!bFound) {
            bSubproduct = false;
            break;
          }
        }

        if (bSubproduct) {
          bCanBeExtended = true;
          break;
        }
      }

      if (!bCanBeExtended) {
        graphsNotExtendable.push_back(gpSmall.get_vertices());
      }
    }
  }

  return graphsNotExtendable;
}

void CoxIter::bCanBeFiniteCovolume_complete_computeGraphsProducts(
    GraphsListIterator grIt, GraphsProduct &gp,
    vector<bool> &bGPVerticesNonLinkable) {
  vector<short unsigned int>::iterator iIt;
  vector<short unsigned int> flaggedVertices;
  unsigned int graphRank(0);

  while (grIt.ptr && (gp.rank + graphRank <= verticesCount)) {
    // ---------------------------------------------------
    // est ce que le graphe est admissible?
    for (iIt = grIt.ptr->vertices.begin(); iIt != grIt.ptr->vertices.end();
         ++iIt) {
      if (bGPVerticesNonLinkable[*iIt]) // si pas linkable
        break;
    }

    // si le graphe est admissible
    if (iIt == grIt.ptr->vertices.end()) {
      // le graphe est ajouté au produit
      gp.graphs.push_back(grIt.ptr);

      // taille du graphe courant
      graphRank = grIt.ptr->vertices.size() - 1;
      gp.rank += graphRank;

#pragma omp critical
      {
        if (2 <= gp.rank && gp.rank <= (dimension - 1))
          graphsProducts_bCanBeFiniteCovolume[gp.rank].push_back(
              GraphsProductSet(gp));
      }

      // mise à jour des sommets que l'on ne peut plus prendre
      for (unsigned int i = 0; i < verticesCount; i++) {
        if (!grIt.ptr->bVerticesLinkable[i] && !bGPVerticesNonLinkable[i]) {
          flaggedVertices.push_back(i);
          bGPVerticesNonLinkable[i] = true;
        }
      }

      // récursion
      bCanBeFiniteCovolume_complete_computeGraphsProducts(
          ++grIt, gp, bGPVerticesNonLinkable);

      // -----------------------------------------------
      // dé-initialisations

      // on remet la liste à son état d'avant la récursion
      for (iIt = flaggedVertices.begin(); iIt != flaggedVertices.end(); ++iIt)
        bGPVerticesNonLinkable[*iIt] = false;

      gp.rank -= graphRank;

      // le graphe est enlevé
      gp.graphs.pop_back();

      if (!gp.graphs.size())
        break;
    } else
      ++grIt;
  }
}

mpz_class CoxIter::i_orderFiniteSubgraph(const unsigned int &type,
                                         const unsigned int &dataSupp) {
  if (type == 0) // A_n
    return factorials[dataSupp + 1];
  else if (type == 1) // Bn
    return (factorials[dataSupp] * powersOf2[dataSupp]);
  else if (type == 3) // Dn
    return (factorials[dataSupp] * powersOf2[dataSupp - 1]);
  else if (type == 4) {
    if (dataSupp == 6)
      return 51840;
    else if (dataSupp == 7)
      return 2903040;
    else if (dataSupp == 8)
      return 696729600;
  } else if (type == 5) // F4
    return 1152;
  else if (type == 6) // G_2^n
    return (2 * dataSupp);
  else if (type == 7) {
    if (dataSupp == 3)
      return 120;
    else if (dataSupp)
      return 14400;
  } else
    throw(0);

  return 0;
}

void CoxIter::growthSeries_mergeTerms(vector<mpz_class> &polynomial,
                                      vector<unsigned int> &symbol,
                                      vector<mpz_class> tempPolynomial,
                                      const vector<unsigned int> &tempSymbol,
                                      mpz_class biTemp) {
  unsigned int symbolMax(symbol.size() ? symbol.size() - 1 : 0);
  unsigned int tempSymbolMax(tempSymbol.size() - 1);

  vector<unsigned int> tempSymbolDenominator(tempSymbol);

  if (tempSymbolMax < symbolMax)
    tempSymbolDenominator.insert(tempSymbolDenominator.end(),
                                 symbolMax - tempSymbolMax, 0);

  // First step to compute the lcm of the two symbols
  for (unsigned int i(1); i <= tempSymbolMax; i++) {
    const auto jMin = i > symbolMax ? 1 : symbol[i] + 1;

    for (unsigned int j(jMin); j <= tempSymbol[i]; j++)
      Polynomials::polynomialDotSymbol(polynomial, i);
  }

  // Second step to compute the lcm of the two symbols
  for (unsigned int i(1); i <= symbolMax; i++) {
    tempSymbolDenominator[i] =
        i > tempSymbolMax ? symbol[i] : max(tempSymbol[i], symbol[i]);

    const auto jMin = i > tempSymbolMax ? 1 : tempSymbol[i] + 1;

    for (unsigned int j(jMin); j <= symbol[i]; j++)
      Polynomials::polynomialDotSymbol(tempPolynomial, i);
  }

  // we eventually add some zeroes
  if (polynomial.size() < tempPolynomial.size())
    polynomial.insert(polynomial.end(),
                      tempPolynomial.size() - polynomial.size(), 0);

  unsigned int tempPolynomialDegree(tempPolynomial.size() - 1);

  // Addition of the two numerators
  for (unsigned int i(0); i <= tempPolynomialDegree; i++)
    polynomial[i] += tempPolynomial[i] * biTemp;

  // ----------------------------------------------------
  // Final stuff
  symbol = tempSymbolDenominator;

  // We remove final 0
  while (polynomial[tempPolynomialDegree] == 0)
    tempPolynomialDegree--;
  polynomial.erase(polynomial.begin() + tempPolynomialDegree + 1,
                   polynomial.end());
}

void CoxIter::growthSeries() {
  if (!bUseOpenMP || verticesCount < 10)
    growthSeries_sequential();
  else
    growthSeries_parallel();

  if (debug)
    growthSeries_details();
}

void CoxIter::growthSeries_details() {
  if (!isGraphExplored)
    exploreGraph();

  if (!isGraphsProductsComputed)
    computeGraphsProducts();

  unsigned int sizeMax(graphsProductsCount_spherical.size());

  unsigned int exponent; // Temporary exponent

  growthSeries_raw = "1";

  string strGrowth_temp, strSymbol;

  for (unsigned int size(1); size < sizeMax; size++) {
    strGrowth_temp = "";

    for (const auto &product :
         graphsProductsCount_spherical[size]) // For each product of that size
    {
      // ----------------------------------------------------
      // Preliminary stuff

      // We compute the symbol and the exponent of this product
      growthSeries_symbolExponentFromProduct(product.first, strSymbol,
                                             exponent);

      mpz_class biTemp((int)product.second * ((size % 2) ? -1 : 1));

      strGrowth_temp +=
          (strGrowth_temp == "" ? "" : " + ") +
          (product.second == 1 ? "" : to_string(product.second) + " * ") +
          "x^" + to_string(exponent) + "/";

      if (strOuputMathematicalFormat == "mathematica")
        strGrowth_temp += "Symb[{" + strSymbol + "},x]";
      else if (strOuputMathematicalFormat == "pari")
        strGrowth_temp += "Symb([" + strSymbol + "],x)";
      else
        strGrowth_temp += "[" + strSymbol + "]";
    }

    if (strGrowth_temp != "")
      growthSeries_raw += ((size % 2) ? " - (" : " + (") + strGrowth_temp + ")";
  }
}

void CoxIter::growthSeries_sequential() {
  if (!isGraphExplored)
    exploreGraph();

  if (!isGraphsProductsComputed)
    computeGraphsProducts();

  unsigned int sizeMax(graphsProductsCount_spherical.size());

  vector<unsigned int> symbol; // Temporary symbol
  unsigned int symbolMax;      // Size of the temporary symbol
  unsigned int exponent;       // Temporary exponent

  vector<unsigned int> growthSeries_symbolNumerator;
  growthSeries_polynomialDenominator = vector<mpz_class>(1, 1);
  growthSeries_cyclotomicNumerator.clear();
  growthSeries_isFractionReduced = true;

  vector<unsigned int> symbolDenominatorTemp;
  unsigned int symbolDenominatorMax(0);

  for (unsigned int size(1); size < sizeMax; size++) {
    for (const auto &product : graphsProductsCount_spherical[size]) {
      // ----------------------------------------------------
      // Preliminary stuff

      // We compute the symbol and the exponent of this product
      growthSeries_symbolExponentFromProduct(product.first, symbol, exponent);

      symbolMax = symbol.size() - 1;
      symbolDenominatorTemp = symbol;

      mpz_class biTemp((int)product.second * ((size % 2) ? -1 : 1));

      vector<mpz_class> tempPolynomial(vector<mpz_class>(exponent, 0));
      tempPolynomial.push_back(1); // x^exponent

      for (unsigned int i(symbolMax + 1); i <= symbolDenominatorMax;
           i++) // we add some zeroes
        symbolDenominatorTemp.push_back(0);

      // ----------------------------------------------------
      // Update
      /*
       * We have here the current rational function:
       * growthSeries_polynomialDenominator / symbolDenominator We want to add
       * the rational function: (-1)^iSize * product.second * x^exponent /
       * symbol
       */

      // First step to compute the lcm of the two symbols symbolDenominator and
      // symbol
      for (unsigned int i(1); i <= symbolMax; i++) {
        const auto jMin =
            i > symbolDenominatorMax ? 1 : growthSeries_symbolNumerator[i] + 1;
        for (unsigned int j(jMin); j <= symbol[i]; j++)
          Polynomials::polynomialDotSymbol(growthSeries_polynomialDenominator,
                                           i);
      }

      // Second step to compute the lcm of the two symbols symbolDenominator
      // and symbol
      for (unsigned int i(1); i <= symbolDenominatorMax; i++) {
        symbolDenominatorTemp[i] =
            max(symbol[i], growthSeries_symbolNumerator[i]);

        for (unsigned int j(i > symbolMax ? 1 : symbol[i] + 1);
             j <= growthSeries_symbolNumerator[i]; j++)
          Polynomials::polynomialDotSymbol(tempPolynomial, i);
      }

      // we eventually add some zeroes
      if (growthSeries_polynomialDenominator.size() < tempPolynomial.size())
        growthSeries_polynomialDenominator.insert(
            growthSeries_polynomialDenominator.end(),
            tempPolynomial.size() - growthSeries_polynomialDenominator.size(),
            0);

      unsigned int tempPolynomialDegree(tempPolynomial.size() - 1);

      // Addition of the two numerators
      for (unsigned int i(0); i <= tempPolynomialDegree; i++)
        growthSeries_polynomialDenominator[i] += tempPolynomial[i] * biTemp;

      // ----------------------------------------------------
      // Final stuff
      growthSeries_symbolNumerator = symbolDenominatorTemp;
      symbolDenominatorMax = growthSeries_symbolNumerator.size() - 1;

      // We remove final 0
      while (growthSeries_polynomialDenominator[tempPolynomialDegree] == 0)
        tempPolynomialDegree--;
      growthSeries_polynomialDenominator.erase(
          growthSeries_polynomialDenominator.begin() + tempPolynomialDegree + 1,
          growthSeries_polynomialDenominator.end());
    }
  }

  // --------------------------------------------------------------
  // Symbols --> Cyclotomic polynomials
  vector<unsigned int> cyclotomicTemp;

  for (unsigned int i(symbolDenominatorMax); i >= 2; i--) {
    if (growthSeries_symbolNumerator[i]) {
      auto divisors(listDivisors(i, true));
      divisors.push_back(i);

      for (unsigned int j(1); j <= growthSeries_symbolNumerator[i]; j++)
        cyclotomicTemp.insert(cyclotomicTemp.end(), divisors.begin(),
                              divisors.end());
    }
  }

  // --------------------------------------------------------------
  // Simplifications
  unsigned int cyclotomicTempSize(cyclotomicTemp.size()),
      cyclotomicMax(Polynomials::cyclotomicPolynomials.size() - 1);

  for (unsigned int i(0); i < cyclotomicTempSize; i++) {
    if (cyclotomicMax < cyclotomicTemp[i] ||
        !Polynomials::dividePolynomialByPolynomial(
            growthSeries_polynomialDenominator,
            Polynomials::cyclotomicPolynomials[cyclotomicTemp[i]]))
      growthSeries_cyclotomicNumerator.push_back(cyclotomicTemp[i]);

    if (cyclotomicMax < cyclotomicTemp[i])
      growthSeries_isFractionReduced = false;
  }

  // --------------------------------------------------------------
  // Final stuff

  // We remove final 0
  while (!growthSeries_symbolNumerator[symbolDenominatorMax])
    symbolDenominatorMax--;
  growthSeries_symbolNumerator.erase(growthSeries_symbolNumerator.begin() +
                                         symbolDenominatorMax + 1,
                                     growthSeries_symbolNumerator.end());

  sort(growthSeries_cyclotomicNumerator.begin(),
       growthSeries_cyclotomicNumerator.end());

  isGrowthSeriesComputed = true;
}

void CoxIter::growthSeries_parallel() {
  if (!isGraphExplored)
    exploreGraph();

  if (!isGraphsProductsComputed)
    computeGraphsProducts();

  unsigned int sizeMax(graphsProductsCount_spherical.size());

  growthSeries_polynomialDenominator.clear();
  growthSeries_cyclotomicNumerator.clear();
  growthSeries_isFractionReduced = true;

  // -----------------------------------------------------------------
  // Private and local variables
  vector<unsigned int> temp_symbolDenominatorTemp;
  vector<unsigned int> symbol; // Temporary symbol
  unsigned int exponent;       // Temporary exponent
  unsigned int threadId;

  // -----------------------------------------------------------------
  // Shared and local variables
  int iOMPMaxThreads(omp_get_max_threads());

  vector<vector<mpz_class>> gs_polynomialDenominator(iOMPMaxThreads,
                                                     vector<mpz_class>(1, 0));
  vector<vector<unsigned int>> gs_symbolNumerator(iOMPMaxThreads,
                                                  vector<unsigned int>(0));

  gs_polynomialDenominator[0][0] =
      1; // Master thread, empty set -> trivial subgroup

#pragma omp parallel for default(none)                                         \
    shared(sizeMax, gs_symbolNumerator, gs_polynomialDenominator) private(     \
        exponent, symbol, threadId, temp_symbolDenominatorTemp)                \
        schedule(static, 1)
  for (unsigned int size = sizeMax - 1; size >= 1; size--) // For each size
  {
    threadId = omp_get_thread_num();

    for (auto product :
         graphsProductsCount_spherical[size]) // For each product of that size
    {
      // ----------------------------------------------------
      // Preliminary stuff

      // We compute the symbol and the exponent of this product
      growthSeries_symbolExponentFromProduct(product.first, symbol, exponent);
      mpz_class biTemp((int)product.second * ((size % 2) ? -1 : 1));

      vector<mpz_class> temp_polynomial(vector<mpz_class>(exponent, 0));
      temp_polynomial.push_back(1); // x^exponent

      growthSeries_mergeTerms(gs_polynomialDenominator[threadId],
                              gs_symbolNumerator[threadId], temp_polynomial,
                              symbol, biTemp);
    }
  }

  // --------------------------------------------------------------
  // Reduction
  growthSeries_polynomialDenominator = gs_polynomialDenominator[0];
  for (int thread(1); thread < iOMPMaxThreads; thread++) {
    if (gs_symbolNumerator[thread].size())
      growthSeries_mergeTerms(
          growthSeries_polynomialDenominator, gs_symbolNumerator[0],
          gs_polynomialDenominator[thread], gs_symbolNumerator[thread]);
  }

  // --------------------------------------------------------------
  // Symbols --> Cyclotomic polynomials
  vector<unsigned int> cyclotomicTemp;
  unsigned int symbolDenominatorMax(gs_symbolNumerator[0].size() - 1);

  for (unsigned int i(symbolDenominatorMax); i >= 2; i--) {
    if (gs_symbolNumerator[0][i]) {
      auto divisors(listDivisors(i, true));
      divisors.push_back(i);

      for (unsigned int j(1); j <= gs_symbolNumerator[0][i]; j++)
        cyclotomicTemp.insert(cyclotomicTemp.end(), divisors.begin(),
                              divisors.end());
    }
  }

  // --------------------------------------------------------------
  // Simplifications
  unsigned int cyclotomicTempSize(cyclotomicTemp.size()),
      cyclotomicMax(Polynomials::cyclotomicPolynomials.size() - 1);
  for (unsigned int i(0); i < cyclotomicTempSize; i++) {
    if (cyclotomicMax < cyclotomicTemp[i] ||
        !Polynomials::dividePolynomialByPolynomial(
            growthSeries_polynomialDenominator,
            Polynomials::cyclotomicPolynomials[cyclotomicTemp[i]]))
      growthSeries_cyclotomicNumerator.push_back(cyclotomicTemp[i]);

    if (cyclotomicMax < cyclotomicTemp[i])
      growthSeries_isFractionReduced = false;
  }

  // --------------------------------------------------------------
  // Final stuff
  sort(growthSeries_cyclotomicNumerator.begin(),
       growthSeries_cyclotomicNumerator.end());
  isGrowthSeriesComputed = true;
}

void CoxIter::get_growthSeries(vector<unsigned int> &cyclotomicNumerator,
                               vector<mpz_class> &polynomialDenominator,
                               bool &isReduced) {
  if (!isGrowthSeriesComputed)
    growthSeries();

  cyclotomicNumerator = growthSeries_cyclotomicNumerator;
  polynomialDenominator = growthSeries_polynomialDenominator;
  isReduced = growthSeries_isFractionReduced;
}

bool CoxIter::get_isGrowthSeriesReduced() {
  if (!isGrowthSeriesComputed)
    growthSeries();

  return growthSeries_isFractionReduced;
}

vector<mpz_class> CoxIter::get_growthSeries_denominator() {
  if (!isGrowthSeriesComputed)
    growthSeries();

  return growthSeries_polynomialDenominator;
}

void CoxIter::growthSeries_symbolExponentFromProduct(
    const vector<vector<short unsigned int>> &product, string &strSymbol,
    unsigned int &exponent) const {
  short unsigned int j, jMax, k;

  vector<unsigned int> symbol;
  exponent = 0;

  for (unsigned int i(0); i < 8;
       i++) // For each type of irreducible spherical graph
  {
    jMax = product[i].size();

    // pour chaque taille
    for (j = 0; j < jMax; j++) {
      if (!product[i][j])
        continue;

      vector<unsigned int> symbolTemp;

      switch (i) {
      case 0: // An
        exponent += product[i][j] * (j + 1) * (j + 2) / 2;
        for (k = 2; k <= j + 2; k++)
          symbolTemp.push_back(k);
        break;
      case 1: // Bn
        exponent += product[i][j] * (j + 1) * (j + 1);
        for (k = 1; k <= j + 1; k++)
          symbolTemp.push_back(2 * k);
        break;
      case 3: // Dn
        exponent += product[i][j] * j * (j + 1);
        for (k = 1; k <= j; k++)
          symbolTemp.push_back(2 * k);
        symbolTemp.push_back(j + 1);
        break;
      case 4:       // En:
        if (j == 5) // E6
        {
          symbolTemp = vector<unsigned int>{2, 5, 6, 8, 9, 12};
          exponent += product[i][j] * 36;
        } else if (j == 6) {
          symbolTemp = vector<unsigned int>{2, 6, 8, 10, 12, 14, 18};
          exponent += product[i][j] * 63;
        } else if (j == 7) {
          symbolTemp = vector<unsigned int>{2, 8, 12, 14, 18, 20, 24, 30};
          exponent += product[i][j] * 120;
        }
        break;
      case 5: // F4
        symbolTemp = vector<unsigned int>{2, 6, 8, 12};
        exponent += product[i][j] * 24;
        break;
      case 6: // G_2^m
        exponent += product[i][j] * (j + 1);
        symbolTemp.push_back(2);
        symbolTemp.push_back(j + 1);
        break;
      case 7:
        if (j == 2) // H3
        {
          symbolTemp = vector<unsigned int>{2, 6, 10};
          exponent += product[i][j] * 15;
        } else if (j == 3) {
          symbolTemp = vector<unsigned int>{2, 12, 20, 30};
          exponent += product[i][j] * 60;
        }
        break;
      }

      for (k = 0; k < product[i][j]; k++)
        symbol.insert(symbol.end(), symbolTemp.begin(), symbolTemp.end());
    }
  }

  sort(symbol.begin(), symbol.end());

  strSymbol = implode(",", symbol);
}

void CoxIter::growthSeries_symbolExponentFromProduct(
    const vector<vector<short unsigned int>> &product,
    vector<unsigned int> &symbol, unsigned int &exponent) const {
  unsigned int j, jMax, k;

  symbol.clear();
  exponent = 0;

  for (unsigned int i(0); i < 8;
       i++) // For each type of irreducible spherical graph
  {
    jMax = product[i].size();

    // pour chaque taille
    for (j = 0; j < jMax; j++) {
      if (!product[i][j])
        continue;

      vector<unsigned int> symbolTemp;

      switch (i) {
      case 0: // An
        exponent += product[i][j] * (j + 1) * (j + 2) / 2;
        for (k = 2; k <= j + 2; k++)
          symbolTemp.push_back(k);
        break;
      case 1: // Bn
        exponent += product[i][j] * (j + 1) * (j + 1);
        for (k = 1; k <= j + 1; k++)
          symbolTemp.push_back(2 * k);
        break;
      case 3: // Dn
        exponent += product[i][j] * j * (j + 1);
        for (k = 1; k <= j; k++)
          symbolTemp.push_back(2 * k);
        symbolTemp.push_back(j + 1);
        break;
      case 4:       // En:
        if (j == 5) // E6
        {
          symbolTemp = vector<unsigned int>{2, 5, 6, 8, 9, 12};
          exponent += product[i][j] * 36;
        } else if (j == 6) {
          symbolTemp = vector<unsigned int>{2, 6, 8, 10, 12, 14, 18};
          exponent += product[i][j] * 63;
        } else if (j == 7) {
          symbolTemp = vector<unsigned int>{2, 8, 12, 14, 18, 20, 24, 30};
          exponent += product[i][j] * 120;
        }
        break;
      case 5: // F4
        symbolTemp = vector<unsigned int>{2, 6, 8, 12};
        exponent += product[i][j] * 24;
        break;
      case 6: // G_2^m
        exponent += product[i][j] * (j + 1);
        symbolTemp.push_back(2);
        symbolTemp.push_back(j + 1);
        break;
      case 7:
        if (j == 2) // H3
        {
          symbolTemp = vector<unsigned int>{2, 6, 10};
          exponent += product[i][j] * 15;
        } else if (j == 3) {
          symbolTemp = vector<unsigned int>{2, 12, 20, 30};
          exponent += product[i][j] * 60;
        }
        break;
      }

      for (auto symb : symbolTemp) {
        for (k = symbol.size(); k <= symb; k++)
          symbol.push_back(0);

        symbol[symb] += product[i][j];
      }
    }
  }
}

bool CoxIter::bEulerCharacteristicFVector() {
  // variables de boucles
  size_t i, j, k, max;
  map<vector<vector<short unsigned int>>, unsigned int>::iterator itMap;

  mpz_class biTemp, biOrderTemp;
  MPZ_rational brAlternateTemp;

  bool bPositive(true);

  fVector = vector<unsigned int>(dimension + 1, 0);
  int fVectorIndex(dimension);

  unsigned int currentVerticesCount(0);

  fVectorAlternateSum = 0;
  brEulerCaracteristic = 1;
  strEulerCharacteristic_computations = "1";

  if (debug)
    cout << "\nProducts of spherical graphs" << endl;

  fVector[dimension] = 1;

  // par taille de nombre de sommets
  for (vector<map<vector<vector<short unsigned int>>, unsigned int>>::iterator
           itMaps(graphsProductsCount_spherical.begin());
       itMaps != graphsProductsCount_spherical.end(); ++itMaps) {
    brAlternateTemp = 0;

    // on parcourt les produits pour la taille donnée
    for (itMap = itMaps->begin(); itMap != itMaps->end(); ++itMap) {
      biTemp = 1;

      if (debug)
        cout << "\t" << currentVerticesCount << ": ";

      // pour chaque type de graphe
      for (i = 0; i < 8; i++) {
        max = (itMap->first[i]).size();

        // pour chaque taille
        for (j = 0; j < max; j++) {
          if (itMap->first[i][j]) {
            if (debug)
              cout << (char)(i + 65) << "_" << (j + 1) << "^"
                   << itMap->first[i][j] << " | ";

            biOrderTemp = i_orderFiniteSubgraph(i, j + 1);
            for (k = 1; k <= itMap->first[i][j]; k++)
              biTemp *= biOrderTemp;
          }
        }
      }

      brAlternateTemp += MPZ_rational(itMap->second, biTemp);

      if (dimension) {
        if (fVectorIndex < 0)
          return false;

        fVector[fVectorIndex] += itMap->second;
      }
      if (debug)
        cout << "N: " << itMap->second << " / Order: " << biTemp.get_str()
             << endl;
    }

    currentVerticesCount++;

    if (bPositive)
      brEulerCaracteristic += brAlternateTemp;
    else
      brEulerCaracteristic -= brAlternateTemp;

    bPositive = !bPositive;
    fVectorIndex--;
  }

  if (dimension && graphsProductsCount_euclidean.size() + 1 > dimension) {
    // si la dimension est spécifiée, on va mettre à jour le f-vecteur et la
    // somme alternée avec le nombre de sommets à l'infini
    verticesAtInfinityCount = 0;
    for (itMap = graphsProductsCount_euclidean[dimension - 1].begin();
         itMap != graphsProductsCount_euclidean[dimension - 1].end(); ++itMap)
      verticesAtInfinityCount += itMap->second;

    fVector[0] += verticesAtInfinityCount;

    for (i = 0; i < dimension; i++)
      fVectorAlternateSum += ((i % 2) ? -1 : 1) * fVector[i];
  }

  return true;
}

// ##################################################################################################################################3
// Affichages

void CoxIter::printEuclideanGraphsProducts(
    vector<map<vector<vector<short unsigned int>>, unsigned int>>
        *graphsProductsCount) {
  // variables de boucles
  size_t i, j, max;
  map<vector<vector<short unsigned int>>, unsigned int>::iterator itMap;

  // par taille de nombre de sommets
  for (vector<map<vector<vector<short unsigned int>>, unsigned int>>::iterator
           itMaps(graphsProductsCount->begin());
       itMaps != graphsProductsCount->end(); ++itMaps) {
    // on parcourt les produits pour la taille donnée
    for (itMap = itMaps->begin(); itMap != itMaps->end(); ++itMap) {
      cout << "\t";
      // pour chaque type de graphe
      for (i = 0; i < 8; i++) {
        max = (itMap->first[i]).size();

        // pour chaque taille
        for (j = 0; j < max; j++) {
          if (itMap->first[i][j])
            cout << "T" << (char)(i + 65) << "_" << j << "^"
                 << itMap->first[i][j] << " | ";
        }
      }

      cout << "N: " << itMap->second << endl;
    }
  }
}

void CoxIter::printCoxeterMatrix() {
  cout << "Coxeter matrix" << endl;

  cout << "\tVertices: ";
  for (vector<string>::const_iterator it(map_vertices_indexToLabel.begin());
       it != map_vertices_indexToLabel.end(); ++it)
    cout << (it == map_vertices_indexToLabel.begin() ? "" : ", ") << *it;
  cout << endl;

  unsigned int i, j;
  if (strOuputMathematicalFormat == "mathematica" ||
      strOuputMathematicalFormat == "gap") {
    cout << "\t[";
    for (i = 0; i < verticesCount; i++) {
      cout << (i ? "," : "") << "[";
      for (j = 0; j < verticesCount; j++) {
        cout << (j ? "," : "")
             << (i == j ? 1
                        : (coxeterMatrix[i][j] < 2 ? 0 : coxeterMatrix[i][j]));
      }
      cout << "]";
    }
    cout << "]" << endl;
  } else {
    for (i = 0; i < verticesCount; i++) {
      cout << "\t";
      for (j = 0; j < verticesCount; j++) {
        cout << (j ? "," : "")
             << (i == j ? 1
                        : (coxeterMatrix[i][j] < 2 ? 0 : coxeterMatrix[i][j]));
      }
      cout << endl;
    }
  }
}

void CoxIter::printCoxeterGraph() {
  cout << "Coxeter graph:\n\t[" << get_strCoxeterGraph() << "]\n" << endl;
}

void CoxIter::printGramMatrix() {
  if (strOuputMathematicalFormat == "gap")
    printGramMatrix_GAP();
  else if (strOuputMathematicalFormat == "latex")
    printGramMatrix_LaTeX();
  else if (strOuputMathematicalFormat == "mathematica")
    printGramMatrix_Mathematica();
  else if (strOuputMathematicalFormat == "pari")
    printGramMatrix_PARI();
  else
    cout << "Gram matrix  \n\t" << get_strGramMatrix() << "\n" << endl;

  unsigned int i, j;
  for (i = 0; i < verticesCount; i++) {
    for (j = 0; j < i; j++) {
      if (coxeterMatrix[i][j] == 1 &&
          strWeights.find(linearizationMatrix_index(j, i, verticesCount)) ==
              strWeights.end()) {
        cout << "l" << j << "m" << i
             << ": weight of the dotted line between hyperplanes "
             << map_vertices_indexToLabel[j] << " and "
             << map_vertices_indexToLabel[i] << endl;
      }
    }
  }

  cout << endl;
}

void CoxIter::printGramMatrix_GAP() {
  cout << "Gram matrix (GAP): \n\t" << get_strGramMatrix_GAP() << "\n" << endl;
}

void CoxIter::printGramMatrix_Mathematica() {
  cout << "Gram matrix (Mathematica): \n\t" << get_strGramMatrix_Mathematica()
       << "\n"
       << endl;
}

void CoxIter::printGramMatrix_PARI() {
  cout << "Gram matrix (PARI): \n\t" << get_strGramMatrix_PARI() << "\n"
       << endl;
}

void CoxIter::printGramMatrix_LaTeX() {
  cout << "Gram matrix (LaTeX): \n\t" << get_strGramMatrix_LaTeX() << "\n"
       << endl;
}

void CoxIter::printEdgesVisitedMatrix() {
  unsigned int i, j;
  cout << "Matrix of visited edges" << endl;

  for (i = 0; i < verticesCount; i++) {
    for (j = 0; j < verticesCount; j++)
      cout << (bEdgesVisited[i][j] ? 1 : 0) << " ";
    cout << endl;
  }
}

bool CoxIter::bIsVertexValid(const string &strVertexLabel) const {
  return (find(map_vertices_indexToLabel.begin(),
               map_vertices_indexToLabel.end(),
               strVertexLabel) != map_vertices_indexToLabel.end());
}

void CoxIter::map_vertices_labels_removeReference(const unsigned int &index) {
  if (index > map_vertices_indexToLabel.size())
    return;

  for (unsigned int i(index + 1); i < verticesCount; i++)
    map_vertices_labelToIndex[map_vertices_indexToLabel[i]]--;

  map_vertices_labelToIndex.erase(map_vertices_indexToLabel[index]);
  map_vertices_indexToLabel.erase(map_vertices_indexToLabel.begin() + index);
}

void CoxIter::map_vertices_labels_addReference(const string &strLabel) {
  map_vertices_labelToIndex[strLabel] = map_vertices_indexToLabel.size();
  map_vertices_indexToLabel.push_back(strLabel);
}

// ##################################################################################################################################3
// Accesseurs
unsigned int CoxIter::get_vertexIndex(const string &strVertexLabel) const {
  map<string, unsigned int>::const_iterator it(
      map_vertices_labelToIndex.find(strVertexLabel));

  if (it == map_vertices_labelToIndex.end())
    throw(string("CoxIter::get_bertexIndex: Invalid vertex name: " +
                 strVertexLabel));

  return it->second;
}

string CoxIter::get_vertexLabel(const unsigned int &vertex) const {
  if (vertex >= verticesCount)
    throw(string("CoxIter::get_strVertexLabel: Invalid vertex index"));

  return map_vertices_indexToLabel[vertex];
}

vector<string> CoxIter::get_str_map_vertices_indexToLabel() const {
  return map_vertices_indexToLabel;
}

vector<vector<unsigned int>> CoxIter::get_coxeterMatrix() const {
  return coxeterMatrix;
}

unsigned int CoxIter::get_coxeterMatrixEntry(const unsigned int &i,
                                             const unsigned int &j) const {
  if (i >= verticesCount || j >= verticesCount)
    throw(string("CoxIter::get_coxeterMatrixEntry: This entry does not exist"));

  return coxeterMatrix[i][j];
}

std::map<unsigned int, string> CoxIter::get_strWeights() const {
  return strWeights;
}

string CoxIter::get_strCoxeterMatrix() const {
  string strCox;

  for (vector<vector<unsigned int>>::const_iterator itRow(
           coxeterMatrix.begin());
       itRow != coxeterMatrix.end(); ++itRow) {
    if (itRow != coxeterMatrix.begin())
      strCox += "\n";

    for (vector<unsigned int>::const_iterator itCol(itRow->begin());
         itCol != itRow->end(); ++itCol)
      strCox += (itCol == itRow->begin() ? "" : ",") + to_string(*itCol);
  }

  return strCox;
}

vector<vector<string>> CoxIter::get_array_str_2_GramMatrix() const {
  size_t i, j;
  vector<vector<string>> strGramMatrix(
      vector<vector<string>>(verticesCount, vector<string>(verticesCount, "")));

  for (i = 0; i < verticesCount; i++) {
    for (j = 0; j <= i; j++) {
      if (i == j)
        strGramMatrix[i][i] = "2";
      else if (coxeterMatrix[i][j] == 0)
        strGramMatrix[i][j] = "-2";
      else if (coxeterMatrix[i][j] == 1)
        strGramMatrix[i][j] =
            "2*l" + to_string(static_cast<long long>(min(i, j))) + "m" +
            to_string(static_cast<long long>(max(i, j)));
      else {
        if (coxeterMatrix[i][j] == 2)
          strGramMatrix[i][j] = "0";
        else if (coxeterMatrix[i][j] == 3)
          strGramMatrix[i][j] = "-1";
        else if (coxeterMatrix[i][j] == 4)
          strGramMatrix[i][j] = "-sqrt(2)";
        else if (coxeterMatrix[i][j] == 5)
          strGramMatrix[i][j] = "-(1+sqrt(5))/2";
        else if (coxeterMatrix[i][j] == 6)
          strGramMatrix[i][j] = "-sqrt(3)";
        else
          strGramMatrix[i][j] =
              "-2*cos(%pi/" +
              to_string(static_cast<long long>(coxeterMatrix[i][j])) + ")";
      }

      strGramMatrix[j][i] = strGramMatrix[i][j];
    }
  }

  return strGramMatrix;
}

string CoxIter::get_strCoxeterGraph() const {
  unsigned int i, j;
  vector<unsigned int> usedVertices;

  string strCoxeterGraph(""), strTemp;

  for (i = 0; i < verticesCount; i++) {
    strTemp = "";
    for (j = i + 1; j < verticesCount; j++) {
      if (coxeterMatrix[i][j] != 2) {
        strTemp += (strTemp == "" ? "[" : ",[") + to_string(j + 1) + "," +
                   to_string(coxeterMatrix[i][j]) + "]";

        auto it(lower_bound(usedVertices.begin(), usedVertices.end(), j));
        if (it == usedVertices.end() || !(*it == j))
          usedVertices.insert(it, j);
      }
    }

    if (strTemp != "") {
      auto it(lower_bound(usedVertices.begin(), usedVertices.end(), i));
      if (it == usedVertices.end() || !(*it == i))
        usedVertices.insert(it, i);

      strCoxeterGraph += (strCoxeterGraph == "" ? "[" : ",[") +
                         to_string(i + 1) + "," + strTemp + "]";
    }
  }

  for (i = 0; i < verticesCount;
       i++) // We display the non-used (i.e. disconnected) vertices
  {
    auto it(lower_bound(usedVertices.begin(), usedVertices.end(), i));
    if (it == usedVertices.end() || !(*it == i))
      strCoxeterGraph += ",[" + to_string(i + 1) + "]";
  }

  return strCoxeterGraph;
}

string CoxIter::get_strGramMatrix() const {
  size_t i, j;
  string strGramMatrix("");

  for (i = 0; i < verticesCount; i++) {
    strGramMatrix += (i ? ", [" : "[ ");
    for (j = 0; j < verticesCount; j++) {
      if (j > 0)
        strGramMatrix += ", ";

      if (i == j)
        strGramMatrix += "1";
      else if (coxeterMatrix[i][j] == 0)
        strGramMatrix += "-1";
      else if (coxeterMatrix[i][j] == 1) {
        map<unsigned int, string>::const_iterator itF(strWeights.find(
            linearizationMatrix_index(min(i, j), max(i, j), verticesCount)));

        if (itF != strWeights.end())
          strGramMatrix += itF->second;
        else
          strGramMatrix += "l" + to_string(static_cast<long long>(min(i, j))) +
                           "m" + to_string(static_cast<long long>(max(i, j)));
      } else {
        if (coxeterMatrix[i][j] == 2)
          strGramMatrix += "0";
        else if (coxeterMatrix[i][j] == 3)
          strGramMatrix += "-1/2";
        else if (coxeterMatrix[i][j] == 4)
          strGramMatrix += "-sqrt(2)/2";
        else if (coxeterMatrix[i][j] == 5)
          strGramMatrix += "-(1+sqrt(5))/4";
        else if (coxeterMatrix[i][j] == 6)
          strGramMatrix += "-sqrt(3)/2";
        else
          strGramMatrix +=
              "-cos(%pi/" +
              to_string(static_cast<long long>(coxeterMatrix[i][j])) + ")";
      }
    }
    strGramMatrix += "]";
  }

  return strGramMatrix;
}

string CoxIter::get_strGramMatrix_LaTeX() const {
  size_t i, j;

  string strGramMatrix("G = \\left(\\begin{array}{*{" +
                       to_string(verticesCount) + "}{c}}");
  for (i = 0; i < verticesCount; i++) {
    strGramMatrix += (i ? "\\\\" : "");
    for (j = 0; j < verticesCount; j++) {
      if (j > 0)
        strGramMatrix += " & ";

      if (i == j)
        strGramMatrix += "1";
      else if (coxeterMatrix[i][j] == 0)
        strGramMatrix += "-1";
      else if (coxeterMatrix[i][j] == 1) {
        map<unsigned int, string>::const_iterator itF(strWeights.find(
            linearizationMatrix_index(min(i, j), max(i, j), verticesCount)));

        if (itF != strWeights.end())
          strGramMatrix += itF->second;
        else
          strGramMatrix += "l" + to_string(static_cast<long long>(min(i, j))) +
                           "m" + to_string(static_cast<long long>(max(i, j)));
      } else {
        if (coxeterMatrix[i][j] == 2)
          strGramMatrix += "0";
        else if (coxeterMatrix[i][j] == 3)
          strGramMatrix += "-\\frac{1}{2}";
        else if (coxeterMatrix[i][j] == 4)
          strGramMatrix += "-\\frac{\\sqrt 2}{2}";
        else if (coxeterMatrix[i][j] == 6)
          strGramMatrix += "-\\frac{\\sqrt 3}{2}";
        else
          strGramMatrix +=
              "-\\cos\\big(\\frac{\\pi}{" +
              to_string(static_cast<long long>(coxeterMatrix[i][j])) +
              "}\\big)";
      }
    }
  }

  strGramMatrix += "\\end{array} \\right)";

  return strGramMatrix;
}

string CoxIter::get_strGramMatrix_Mathematica() const {
  size_t i, j;

  string strGramMatrix("G := {");
  for (i = 0; i < verticesCount; i++) {
    strGramMatrix += (i ? ", {" : "{ ");
    for (j = 0; j < verticesCount; j++) {
      if (j > 0)
        strGramMatrix += ", ";

      if (i == j)
        strGramMatrix += "1";
      else if (coxeterMatrix[i][j] == 0)
        strGramMatrix += "-1";
      else if (coxeterMatrix[i][j] == 1) {
        map<unsigned int, string>::const_iterator itF(strWeights.find(
            linearizationMatrix_index(min(i, j), max(i, j), verticesCount)));

        if (itF != strWeights.end())
          strGramMatrix += itF->second;
        else
          strGramMatrix += "l" + to_string(static_cast<long long>(min(i, j))) +
                           "m" + to_string(static_cast<long long>(max(i, j)));
      } else {
        if (coxeterMatrix[i][j] == 2)
          strGramMatrix += "0";
        else if (coxeterMatrix[i][j] == 3)
          strGramMatrix += "-1/2";
        else if (coxeterMatrix[i][j] == 4)
          strGramMatrix += "-Sqrt[2]/2";
        else if (coxeterMatrix[i][j] == 6)
          strGramMatrix += "-Sqrt[3]/2";
        else
          strGramMatrix +=
              "-Cos[Pi/" +
              to_string(static_cast<long long>(coxeterMatrix[i][j])) + "]";
      }
    }
    strGramMatrix += "}";
  }

  strGramMatrix += "};";

  return strGramMatrix;
}

string CoxIter::get_strGramMatrix_PARI() const {
  size_t i, j;
  string strGramMatrix("G = [");

  for (i = 0; i < verticesCount; i++) {
    strGramMatrix += (i ? "; " : " ");
    for (j = 0; j < verticesCount; j++) {
      if (j > 0)
        strGramMatrix += ", ";

      if (i == j)
        strGramMatrix += "1";
      else if (coxeterMatrix[i][j] == 0)
        strGramMatrix += "-1";
      else if (coxeterMatrix[i][j] == 1) {
        map<unsigned int, string>::const_iterator itF(strWeights.find(
            linearizationMatrix_index(min(i, j), max(i, j), verticesCount)));

        if (itF != strWeights.end())
          strGramMatrix += itF->second;
        else
          strGramMatrix += "l" + to_string(static_cast<long long>(min(i, j))) +
                           "m" + to_string(static_cast<long long>(max(i, j)));
      } else {
        if (coxeterMatrix[i][j] == 2)
          strGramMatrix += "0";
        else if (coxeterMatrix[i][j] == 3)
          strGramMatrix += "-1/2";
        else if (coxeterMatrix[i][j] == 4)
          strGramMatrix += "-sqrt(2)/2";
        else if (coxeterMatrix[i][j] == 5)
          strGramMatrix += "-(1+sqrt(5))/4";
        else if (coxeterMatrix[i][j] == 6)
          strGramMatrix += "-sqrt(3)/2";
        else
          strGramMatrix +=
              "-cos(Pi/" +
              to_string(static_cast<long long>(coxeterMatrix[i][j])) + ")";
      }
    }
  }

  return (strGramMatrix + "];");
}

string CoxIter::get_strGramMatrix_GAP() const {
  size_t i, j;
  string strGramMatrix("G := [ [");

  for (i = 0; i < verticesCount; i++) {
    strGramMatrix += (i ? "], [" : " ");
    for (j = 0; j < verticesCount; j++) {
      if (j > 0)
        strGramMatrix += ", ";

      if (i == j)
        strGramMatrix += "1";
      else if (coxeterMatrix[i][j] == 0)
        strGramMatrix += "-1";
      else if (coxeterMatrix[i][j] == 1) {
        map<unsigned int, string>::const_iterator itF(strWeights.find(
            linearizationMatrix_index(min(i, j), max(i, j), verticesCount)));

        if (itF != strWeights.end())
          strGramMatrix += itF->second;
        else
          strGramMatrix += "l" + to_string(static_cast<long long>(min(i, j))) +
                           "m" + to_string(static_cast<long long>(max(i, j)));
      } else {
        if (coxeterMatrix[i][j] == 2)
          strGramMatrix += "0";
        else if (coxeterMatrix[i][j] == 3)
          strGramMatrix += "-1/2";
        else if (coxeterMatrix[i][j] == 4)
          strGramMatrix += "-Sqrt(2)/2";
        else if (coxeterMatrix[i][j] == 5)
          strGramMatrix += "-(1+Sqrt(5))/4";
        else if (coxeterMatrix[i][j] == 6)
          strGramMatrix += "-Sqrt(3)/2";
        else
          strGramMatrix +=
              "-Cos(FLOAT.PI/" +
              to_string(static_cast<long long>(coxeterMatrix[i][j])) + ")";
      }
    }
  }

  return (strGramMatrix + "] ];");
}

string CoxIter::get_strGramMatrixField() const {
  return (isGramMatrixFieldKnown ? strGramMatrixField : "");
}

MPZ_rational CoxIter::get_brEulerCaracteristic() const {
  return brEulerCaracteristic;
}

string CoxIter::get_strEulerCaracteristic() const {
  return brEulerCaracteristic.to_string();
}

string CoxIter::get_strEulerCharacteristic_computations() const {
  return strEulerCharacteristic_computations;
}

int CoxIter::get_fVectorAlternateSum() const { return fVectorAlternateSum; }

bool CoxIter::get_bWriteInfo() const { return bWriteInfo; }

void CoxIter::set_bWriteInfo(const bool &bNewValue) { bWriteInfo = bNewValue; }

bool CoxIter::get_debug() const { return debug; }

vector<unsigned int> CoxIter::get_fVector() const { return fVector; }

vector<unsigned int> CoxIter::get_infSeqFVectorsUnits() const {
  return infSeqFVectorsUnits;
}

vector<unsigned int> CoxIter::get_infSeqFVectorsPowers() const {
  return infSeqFVectorsPowers;
}

unsigned int CoxIter::get_verticesAtInfinityCount() const {
  return verticesAtInfinityCount;
}

unsigned int CoxIter::get_iIrreducibleSphericalGraphsCount() const {
  return graphsList_spherical->totalGraphsCount;
}

int CoxIter::get_isCocompact() {
  if (isCocompact == -2 && checkCocompactness)
    return isGraphCocompact();

  return isCocompact;
}

int CoxIter::get_isFiniteCovolume() {
  if (isFiniteCovolume == -2 && checkCofiniteness)
    return checkCovolumeFiniteness();

  return isFiniteCovolume;
}

int CoxIter::get_isArithmetic() const { return isArithmetic; }

unsigned int CoxIter::get_dimension() const { return dimension; }

bool CoxIter::get_dimensionGuessed() const { return isDimensionGuessed; }

string CoxIter::get_strError() const { return strError; }

unsigned int CoxIter::get_verticesCount() const { return verticesCount; }

bool CoxIter::get_bHasDottedLine() const { return hasDottedLine; }

int CoxIter::get_iHasDottedLineWithoutWeight() const {
  return hasDottedLineWithoutWeight;
}

GraphsList *CoxIter::get_gl_graphsList_spherical() const {
  return graphsList_spherical;
}

GraphsList *CoxIter::get_gl_graphsList_euclidean() const {
  return graphsList_euclidean;
}

bool CoxIter::get_hasSphericalGraphsOfRank(const unsigned int &rank) const {
  if (rank > verticesCount)
    return false;

  return (graphsProductsCount_spherical[rank].size() != 0);
}

bool CoxIter::get_hasEuclideanGraphsOfRank(const unsigned int &rank) const {
  if (rank > verticesCount)
    return false;

  return (graphsProductsCount_euclidean[rank - 1].size() != 0);
}

void CoxIter::set_isArithmetic(const unsigned int &arithmetic) {
  isArithmetic = arithmetic;
}

void CoxIter::set_checkCocompactness(const bool &value) {
  checkCocompactness = value;
}

void CoxIter::set_checkCofiniteness(const bool &value) {
  checkCofiniteness = value;
}

void CoxIter::set_debug(const bool &value) { debug = value; }

void CoxIter::set_useOpenMP(const bool &value) {
#ifdef _OPENMP
  bUseOpenMP = value;
#endif
}

void CoxIter::set_sdtoutToFile(const string &strFilename) {
  string strOutputCoutFilename(strFilename);
  outCout = new ofstream(strOutputCoutFilename.c_str());

  if (outCout->is_open()) {
    sBufOld = cout.rdbuf(outCout->rdbuf());
    bCoutFile = true;
  }
}

void CoxIter::set_dimension(const unsigned int &dimension_) {
  dimension = dimension_;
  maximalSubgraphRank = dimension ? dimension : verticesCount;
}

void CoxIter::set_coxeterMatrix(const vector<vector<unsigned int>> &iMat) {
  strWeights.clear();
  verticesCount = iMat.size();

  initializations();
  coxeterMatrix = iMat;

  isGraphExplored = false;
  isGraphsProductsComputed = false;
  isGrowthSeriesComputed = false;

  euclideanMaxRankFound = 0;
  sphericalMaxRankFound = 0;
  isDimensionGuessed = false;
  hasDottedLineWithoutWeight = -1;
  hasBoldLine = false;

  maximalSubgraphRank = dimension ? dimension : verticesCount;
}

void CoxIter::set_strOuputMathematicalFormat(const string &strO) {
  strOuputMathematicalFormat = strO;
}

void CoxIter::set_strVerticesToConsider(
    const vector<string> &strVerticesToConsider) {
  strVertices = strVerticesToConsider;
  sort(strVertices.begin(), strVertices.end());
  strVertices = vector<string>(strVertices.begin(),
                               unique(strVertices.begin(), strVertices.end()));
}

void CoxIter::set_strVerticesToRemove(
    const vector<string> &strVerticesRemove_) {
  strVerticesRemove = strVerticesRemove_;
  sort(strVerticesRemove.begin(), strVerticesRemove.end());
  strVerticesRemove = vector<string>(
      strVerticesRemove.begin(),
      unique(strVerticesRemove.begin(), strVerticesRemove.end()));
}
