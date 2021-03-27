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

#include "signature.h"

Signature::Signature() {
  pari_init(50000000, 2);

  /*
   * Note: gEpsilon must be BIG compared to the precision up to which we compute
   * the eigenvalues
   * */
  gEpsilon = dbltor(1e-40);
}

Signature::~Signature() { pari_close(); }

std::array<unsigned int, int(3)> Signature::computeSignature(string matrix) {
  array<unsigned int, int(3)> signature({0, 0, 0});
  GEN gMatrix;

  long prec;
  setrealprecision(57, &prec); // increase precision for the gp_read_str

  pari_CATCH(CATCH_ALL) {
    throw(string(
        "Signature::computeSignature: Incorrect matrix; check the weights"));
  }
  pari_TRY { gMatrix = gp_read_str(matrix.c_str()); }
  pari_ENDCATCH

      GEN gEigenvalues(gel(jacobi(gMatrix, 20), 1)),
      gTemp;
  long int iEigenvaluesCount(lg(gEigenvalues));

  for (long int i(1); i < iEigenvaluesCount; i++) {
    gTemp = gel(gEigenvalues, i);

    if (mpcmp(absr(gTemp), gEpsilon) < 0)
      signature[2]++;
    else if (mpcmp(gadd(gTemp, gEpsilon), gen_0) < 0)
      signature[1]++;
    else if (mpcmp(gsub(gTemp, gEpsilon), gen_0) > 0)
      signature[0]++;
  }

  return signature;
}
