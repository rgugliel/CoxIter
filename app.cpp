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

#include "app.h"

App::App()
    : bCoutFile(false), bOutputGraphToDraw(false), bOutputGraph(false),
      checkCanBeFiniteCovolume(false), checkCocompacity(false),
      checkFiniteCovolume(false), checkArithmeticity(false), computeEuler(true),
      computeGrowthRate(false), computeGrowthSeries(false),
      computeSignature(false), debug(false), bIndex2(false), useOpenMP(true),
      printCoxeterGraph(false), printCoxeterMatrix(false),
      printGramMatrix(false), bPrintHelp(false),
      ouputMathematicalFormat("generic") {}

bool App::readMainParameters(int argc, char **argv) {
  string temp, prevType, complete;

  if (argc == 1 &&
      isatty(fileno(stdin)) ==
          0) // If no parameter was given but a file was passed via stdin
  {
    checkCocompacity = true;
    checkFiniteCovolume = true;
    computeGrowthSeries = true;
#ifdef _COMPILE_WITH_PARI_
    computeGrowthRate = true;
    computeSignature = true;
#endif
    return true;
  }

  for (int i = 0; i < argc; ++i) {
    temp = std::string(argv[i]);

    complete += temp + " ";

    if (temp == "-i") // input
      prevType = "i";
    else if (temp == "-arithmeticity" || temp == "-a") {
      checkCocompacity = true;
      checkArithmeticity = true;
      prevType = "arithmeticity";
    } else if (temp == "-c" || temp == "-compacity" || temp == "-compactness" ||
               temp == "-compact" || temp == "-cocompact") {
      checkCocompacity = true;
      prevType = "compacity";
    } else if (temp == "-cf" || temp == "-of") // force le cout dans un fichier
    {
      bCoutFile = true;
      prevType = "cf";
    } else if (temp == "-debug") {
      debug = true;
      prevType = "debug";
    } else if (temp == "-drawgraph" || temp == "-dg") // output
    {
      bOutputGraphToDraw = true;
      prevType = "drawgraph";
    } else if (temp == "-drop" || temp == "-remove") // we remove a vertex
      prevType = "drop";
    else if (temp == "-ffv") {
      checkCanBeFiniteCovolume = true;
      prevType = "ffv";
    } else if (temp == "-fcv") {
      checkFiniteCovolume = true;
      prevType = "fv";
    } else if (temp == "-full") {
      checkCocompacity = true;
      checkFiniteCovolume = true;
      computeGrowthSeries = true;
#ifdef _COMPILE_WITH_PARI_
      computeGrowthRate = true;
      computeSignature = true;
#endif
      prevType = "full";
    } else if (temp == "-fv") {
      checkFiniteCovolume = true;
      prevType = "fv";
    } else if (temp == "-g" || temp == "-growth" || temp == "-poincarre") {
      computeGrowthSeries = true;
      prevType = "growth";
    } else if (temp == "-growthrate" || temp == "-gr") {
#ifdef _COMPILE_WITH_PARI_
      computeGrowthSeries = true;
      computeGrowthRate = true;
      prevType = "growthrate";
#endif
    } else if (temp == "-help") {
      bPrintHelp = true;
      prevType = "help";
    } else if (temp == "-ne" ||
               temp == "-neuler") // don't compute the Euler characteristic
    {
      computeEuler = false;
      prevType = "ne";
    } else if (temp == "-nopenmp" || temp == "-nparallel") {
      useOpenMP = false;
      prevType = "npenmp";
    } else if (temp == "-o") // Files for the output
    {
      prevType = "o";
    } else if (temp == "-pcg") // print Coxeter matrix
    {
      printCoxeterGraph = true;
      prevType = "pcg";
    } else if (temp == "-pcm") // print Coxeter matrix
    {
      printCoxeterMatrix = true;
      prevType = "pcm";
    } else if (temp == "-pgm") // print Gram matrix?
    {
      printGramMatrix = true;
      prevType = "pgm";
    } else if (temp == "-oformat" || temp == "-outputformat") {
      prevType = "oformat";
    } else if (temp == "-s" || temp == "-signature") {
#ifdef _COMPILE_WITH_PARI_
      computeSignature = true;
      prevType = "signature";
#endif
    } else if (temp == "-writegraph" || temp == "-wg") // write the graph
    {
      bOutputGraph = true;
      prevType = "wg";
    } else {
      if (prevType == "i")
        inFilename = temp;
      else if (prevType == "drop")
        verticesToRemove.push_back(temp);
      else if (prevType == "o")
        outFilenameBasis = temp;
      else if (prevType == "oformat") {
        transform(temp.begin(), temp.end(), temp.begin(), ::tolower);
        if (temp == "gap" || temp == "latex" || temp == "mathematica" ||
            temp == "pari")
          ouputMathematicalFormat = temp;
      }

      prevType = "";
    }
  }

  // --------------------------------------------------
  // Complete parameters
  string params;
  for (int i(0); i < argc; ++i)
    params += " " + std::string(argv[i]);

  // --------------------------------------------------
  // Choosing the vertices
  PCRERegexp regexp;
  PCREResult regexpRes;
  int regexpCount;

  if ((regexpCount =
           regexp.preg_match_all("-S=\\{([[:alnum:][:space:],_-]+)\\}", params,
                                 regexpRes, PCRE_CASELESS)) > 0) {
    str_replace(regexpRes[1][0], " ", "");
    vector<string> elements(explode(",", regexpRes[1][0]));

    for (const auto &vertex : elements)
      vertices.push_back(vertex);
  } else if ((regexpCount =
                  regexp.preg_match_all("-S=([[:alnum:],_-]+)", params,
                                        regexpRes, PCRE_CASELESS)) > 0) {
    for (int i(0); i < regexpCount; i++)
      vertices.push_back(regexpRes[1][i]);
  }

  // --------------------------------------------------
  // Index2
  if ((regexpCount = regexp.preg_match_all("-is=\\[?([[:alnum:],_-]+)\\]?",
                                           params, regexpRes, PCRE_CASELESS)) >
      0) {
    str_replace(regexpRes[1][0], " ", "");
    vector<string> elements(explode(",", regexpRes[1][0]));
    index2vertex_t0 = elements[0];
    if (elements.size() > 1)
      index2vertex_s0 = elements[1];

    str_replace(params, regexpRes[0][0], "");
    bIndex2 = true;
  } else if ((regexpCount = regexp.preg_match_all(
                  "(-gbd|-index2) ([[:alnum:]_-]+)", params, regexpRes,
                  PCRE_CASELESS)) > 0) {
    index2vertex_t0 = regexpRes[2][0];
    str_replace(params, regexpRes[0][0], "");
    bIndex2 = true;
  }

  return true;
}

// TODO: si exposant > 50, afficher erreur pour growth rate (car facteurs pas
// forcément premiers entre eux)

void App::extractIndex2Subgroup(CoxIter &ci) {
  unsigned int verticesCountStart(ci.get_verticesCount());

  if (index2vertex_t0 == "") {
    cout << "Error: A vertex must be given" << endl;
    return;
  }

  Index2 idx2(&ci);

  if (index2vertex_s0 != "") {
    cout << "This is an experimental feature (to be properly tested)" << endl;
    cout << "------------------------------------------------------\n" << endl;

    cout << "Infinite sequence:" << endl;

    if (!idx2.isVertexAdmissible(index2vertex_t0)) {
      cout << "\tError: " << idx2.get_error() << endl;
      return;
    }

    if (!idx2.isVertexAdmissible(index2vertex_s0)) {
      cout << "\tError: " << idx2.get_error() << endl;
      return;
    }

    if (!ci.get_dimension()) {
      cout << "\tError: The dimension must be specified" << endl;
      return;
    }

    if (ci.get_coxeterMatrixEntry(ci.get_vertexIndex(index2vertex_t0),
                                  ci.get_vertexIndex(index2vertex_s0)) >= 2) {
      cout << "\tError: The two hyperplanes must be (ultra)parallel " << endl;
      return;
    }

    ci.IS_computations(index2vertex_t0, index2vertex_s0);

    auto fVUnits(ci.get_infSeqFVectorsUnits());
    auto fVPowers(ci.get_infSeqFVectorsPowers());

    cout << "\tf-vector after n doubling:\n\t(" << implode(", ", fVUnits)
         << ", 1) + 2^(n-1)*(" << implode(", ", fVPowers) << ", 0)" << endl;
  } else
    cout << "Index two subgroup:" << endl;

  // -----------------------------------------
  // Doing the GBD
  if (!idx2.removeVertex(index2vertex_t0))
    cout << "\tError: " << idx2.get_error() << endl;
  else
    cout << "\tNumber of new hyperplanes after first doubling: "
         << (ci.get_verticesCount() - verticesCountStart) << endl;

  cout << endl;
}

void App::run() {
  if (bPrintHelp) {
    printHelp();
    return;
  }

  chrono::time_point<std::chrono::system_clock> timeStart, timeEnd;

  bool isEulerSuccess(true), canBeFiniteCovolume(false),
      isSignatureComputed(false);

  CoxIter ci;
  ci.set_checkCocompactness(checkCocompacity);
  ci.set_checkCofiniteness(checkFiniteCovolume);
  ci.set_debug(debug);
  ci.set_bWriteInfo(true);
  ci.set_ouputMathematicalFormat(ouputMathematicalFormat);
  ci.set_verticesToConsider(vertices);
  ci.set_verticesToRemove(verticesToRemove);

  if (bCoutFile)
    ci.set_sdtOutToFile(outFilenameBasis + ".output");

  Arithmeticity arithmeticity;

#ifdef _COMPILE_WITH_PARI_
  GrowthRate_Result grr;
  grr.isComputed = false;
  grr.perron = -1;
  grr.pisot = -1;
  grr.salem = -1;

  array<unsigned int, 3> signature;
#endif

  if (isatty(fileno(stdin)) == 0) {
    if (!ci.parseGraph(std::cin)) {
      cout << "Error while reading graph: " << ci.get_error() << endl;
      return;
    }
  } else if (inFilename != "") {
    // Reading of the graph
    if (!ci.readGraphFromFile(inFilename)) {
      cout << "Error while reading file: " << ci.get_error() << endl;
      return;
    }
  } else {
    cout << "No input file given\n" << endl;

    printHelp();
    return;
  }

  if (bIndex2)
    extractIndex2Subgroup(ci);

  if (printGramMatrix)
    ci.printGramMatrix();

  if (printCoxeterMatrix)
    ci.printCoxeterMatrix();

  if (printCoxeterGraph)
    ci.printCoxeterGraph();

  if (bOutputGraph && !ci.writeGraph(outFilenameBasis))
    cout << "Error while writing file: " << ci.get_error() << endl;

  if (bOutputGraphToDraw) {
    if (ci.writeGraphToDraw(outFilenameBasis)) {
      string command("dot -Tjpg -o\"" + outFilenameBasis + ".jpg\" \"" +
                     outFilenameBasis + ".graphviz\"");
#ifdef _DOT_PROGRAM_FOUND_
      FILE *fin;
      if ((fin = popen(command.c_str(), "r"))) {
        cout << "Image created: \n\t" << outFilenameBasis << ".jpg\n" << endl;
        pclose(fin);
      } else
        cout << "GraphViz command: \n\t" << command << "\n" << endl;
#else
      cout << "GraphViz command: \n\t" << command << "\n" << endl;
#endif
    } else
      cout << "Error while writing file: " << ci.get_error() << endl;
  }

  if (bIndex2)
    return;

  // -----------------------------------------------------------------
  // composantes connexes sphériques et euclidiennes
  cout << "Finding connected subgraphs......" << endl;
  timeStart = chrono::system_clock::now();
  ci.exploreGraph();

  try {
    if (checkCanBeFiniteCovolume)
      canBeFiniteCovolume = ci.canBeFiniteCovolume();
  } catch (const string &ex) {
    checkCanBeFiniteCovolume = false;
    cout << "\nError:\n\t" << ex << "\n" << endl;
  }

  // -----------------------------------------------------------------
  // calcul des produits (la majorité du temps de calcul concerne ce bloc)
  if (computeEuler || computeGrowthSeries || checkCocompacity ||
      checkFiniteCovolume) {
    cout << "Finding graphs products......" << endl;
    ci.computeGraphsProducts();
  }

  unsigned int dimension(ci.get_dimension());

  if (checkFiniteCovolume)
    ci.checkCovolumeFiniteness();

  // -----------------------------------------------------------------
  // calcul de la caractéristique d'Euler, f-vecteur et compacité
  cout << "Computations......" << endl;
  if (computeEuler && !ci.computeEulerCharacteristicFVector()) {
    isEulerSuccess = false;
    cout << "\n\n##############################################################"
            "########"
         << endl;
    cout << "\tAn error occurred." << endl;
    cout << "\tCheck the graph encoding." << endl;
    cout << "\tYou can run CoxIter with option '-debug' to see if the graph "
            "contains a spherical subgraph which has too big rank."
         << endl;
    cout << "##################################################################"
            "####\n"
         << endl;
  }

  if (checkCocompacity)
    ci.isGraphCocompact();

  if (checkArithmeticity)
    arithmeticity.test(ci, true);

  if (computeGrowthSeries) {
    ci.growthSeries();
  }

  if (ci.get_hasDottedLineWithoutWeight() != 0)
    computeSignature = false;

  if (computeSignature) {
#ifdef _COMPILE_WITH_PARI_
    try {
      Signature s;
      signature = s.computeSignature(ci.get_gramMatrix_PARI());
      isSignatureComputed = true;
    } catch (const string &ex) {
      cout << "\n---------------------------------------------------------"
           << endl;
      cout << "Error while computing the signature:\n\t" << ex << endl;
      cout << "---------------------------------------------------------\n"
           << endl;
    }
#endif
  }

  if (computeGrowthRate) {
#ifdef _COMPILE_WITH_PARI_
    try {
      GrowthRate gr;
      grr = gr.grrComputations(ci.get_growthSeries_denominator());
    } catch (const string &ex) {
      cout << "\n---------------------------------------------------------"
           << endl;
      cout << "Error while computing the growth rate:\n\t" << ex << endl;
      cout << "---------------------------------------------------------\n"
           << endl;
    }
#endif
  }

  timeEnd = chrono::system_clock::now();
  cout << "\tComputation time: "
       << chrono::duration<double, milli>(timeEnd - timeStart).count() / 1000
       << "s\n"
       << endl;

  // -----------------------------------------------------------------
  // Affichage des informations
  cout << "Information" << endl;

  if (ci.get_dimensionGuessed())
    cout << "\tGuessed dimension: " << ci.get_dimension() << endl;

  cout << "\tCocompact: "
       << (ci.get_isCocompact() >= 0
               ? (ci.get_isCocompact() == 0 ? "no" : "yes")
               : "?")
       << endl;
  if (checkCanBeFiniteCovolume)
    cout << "\tCan be of finite covolume: "
         << (canBeFiniteCovolume ? "yes" : "no") << endl;

  cout << "\tFinite covolume: "
       << (ci.get_isFiniteCovolume() >= 0
               ? (ci.get_isFiniteCovolume() == 0 ? "no" : "yes")
               : "?")
       << endl;
  if (checkArithmeticity) {
    cout << "\tArithmetic: ";
    if (ci.get_isArithmetic() == 1)
      cout << "yes" << endl;
    else if (ci.get_isArithmetic() == 0)
      cout << "no" << endl;
    else
      cout << "?"
           << (arithmeticity.get_error() != ""
                   ? "(" + arithmeticity.get_error() + ")"
                   : (ci.get_hasDottedLine() ? " (GRAPH HAS DOTTED EDGE)" : ""))
           << endl;
  }

  // f-vector, alternating sum of the components of the f-vector
  if (computeEuler && isEulerSuccess) {
    if (dimension) {
      const auto fVector(ci.get_fVector());

      cout << "\tf-vector: (";
      for (unsigned int i(0); i <= dimension; i++)
        cout << (i ? ", " : "") << fVector[i];
      cout << ")" << endl;

      cout << "\tNumber of vertices at infinity: "
           << ci.get_verticesAtInfinityCount() << endl;

      cout << "\tAlternating sum of the components of the f-vector: "
           << ci.get_fVectorAlternateSum() << endl;
    }

    cout << "\tEuler characteristic: " << ci.get_brEulerCaracteristic() << endl;
  }

  // volume
  if (computeEuler && dimension && isEulerSuccess && !(dimension % 2) &&
      ci.get_isFiniteCovolume() == 1) {
    cout << "\tCovolume: ";

    MPZ_rational cov((dimension / 2) % 2 ? -1 : 1);
    for (unsigned int i(1); i <= dimension; i++) {
      cov *= 2;
      cov /= i;
      if (i <= (dimension / 2))
        cov *= i;
    }

    cout << "pi^" << (dimension / 2) << " * "
         << cov * ci.get_brEulerCaracteristic() << endl;
  }

  if (isSignatureComputed) {
#ifdef _COMPILE_WITH_PARI_
    cout << "\tSignature (numerically): " << signature[0] << "," << signature[1]
         << "," << signature[2] << endl;
#endif
  }

  if (computeGrowthSeries) {
    cout << "\nGrowth series: " << endl;
    ci.printGrowthSeries();
    cout << endl;

#ifdef _COMPILE_WITH_PARI_
    if (computeGrowthRate && grr.isComputed && ci.get_isGrowthSeriesReduced()) {
      cout << "\nGrowth rate: " << grr.growthRate << endl;
      cout << "\tPerron number: "
           << (grr.perron < 0 ? "?" : (grr.perron > 0 ? "yes" : "no")) << endl;
      cout << "\tPisot number: "
           << (grr.pisot < 0 ? "?" : (grr.pisot > 0 ? "yes" : "no")) << endl;
      cout << "\tSalem number: "
           << (grr.salem < 0 ? "?" : (grr.salem > 0 ? "yes" : "no")) << endl;
    }
#endif
  }

  if (ci.get_isArithmetic() == -1) {
    vector<string> cycles(arithmeticity.get_allCycles());

    if (cycles.size()) {
      cout << "\nThe group is arithmetic if and only if all the following "
              "values lie in Z: \n";
      auto weights(ci.get_weights());

      cout << implode("\n", cycles) << endl;
      if (weights.size()) {
        cout << "with" << endl;
        for (auto it : weights)
          cout << "l"
               << linearizationMatrix_row(it.first, ci.get_verticesCount())
               << "m"
               << linearizationMatrix_col(it.first, ci.get_verticesCount())
               << " = " << it.second << endl;
      }
    }
  }

  cout << endl;
}

void App::printHelp() const {
  cout << "  _____          _____ _            \n"
          " / ____|        |_   _| |\n"
          "| |     _____  __ | | | |_ ___ _ __ \n"
          "| |    / _ \\ \\/ / | | | __/ _ \\ '__|\n"
          "| |___| (_) >  < _| |_| ||  __/ |\n"
          " \\_____\\___/_/\\_\\_____|\\__\\___|_|   \n"
       << endl;

  cout << "CoxIter is a program to compute invariants of hyperbolic Coxeter "
          "groups\n\n"
          "The basis usage is as follows:\n\t./coxiter < "
          "file-describing-the-graph\n\n"
          "For example, you could do\n\t./coxiter < "
          "../graphs/8-Bugaenko.coxiter\n\n"
          "You can also choose to compute/test only some of the invariants:\n"
          "\t-a  : test whether the group is arithmetic\n"
          "\t-c  : test whether the group is cocompact\n"
          "\t-fv : test whether the group has finite covolume\n"
          "\t-g  : growth series\n"
          "\t-gr : growth rate\n"
          "\t-s  : signature\n\n"
          "There are many more options regarding the format of the output and "
          "the possible computations.\nThe full documentation is available "
          "here:\n\thttps://rgugliel.github.io/CoxIter/"
       << endl;
}
