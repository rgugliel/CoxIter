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

#include "mpz_rational.h"

template <typename type> type ugcd(type u, type v) {
  while (v != 0) {
    type r = u % v;
    u = v;
    v = r;
  }

  return u;
}

MPZ_rational::MPZ_rational() : a(0), b(1) { update(); }

MPZ_rational::MPZ_rational(const int &i) : a(i), b(1) { update(); }

MPZ_rational::MPZ_rational(mpz_class a) : a(a), b(1) { update(); }

MPZ_rational::MPZ_rational(mpz_class a, mpz_class b) : a(a), b(b) { update(); }

#ifndef _COMPILE_WITHOUT_REGEXP_
MPZ_rational::MPZ_rational(string rational) {
  str_replace(rational, " ", "");

  PCRERegexp reg;
  PCREResult results;
  int iRegexpCount;

  // rationnel sous la forme: a
  iRegexpCount = reg.preg_match_all("^[\\+]{0,1}([\\-]{0,1})([[:digit:]]+)$",
                                    rational, results);
  if (iRegexpCount == -1)
    throw(0);
  else if (iRegexpCount > 0) {
    a = mpz_class(results[2][0].c_str(), 10);

    b = 1;

    if (b == 0)
      throw(0);

    if (results[1][0] == "-")
      a *= -1;

    update();
    return;
  }

  // rationnel sous la forme: a/b
  iRegexpCount = reg.preg_match_all(
      "^[\\+]{0,1}([\\-]{0,1})([[:digit:]]+)\\/([[:digit:]]+)$", rational,
      results);
  if (iRegexpCount == -1)
    throw(0);
  else if (iRegexpCount > 0) {
    a = mpz_class(results[2][0], 10);
    b = mpz_class(results[3][0], 10);

    if (b == 0)
      throw(0);

    if (results[1][0] == "-")
      a *= -1;

    update();
    return;
  }

  throw(0); // raté
}
#endif

void MPZ_rational::update() {
  mpz_class iGCD;

  mpz_gcd(iGCD.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());

  a /= iGCD;
  b /= iGCD;

  if (b < 0) {
    a *= -1;
    b *= -1;
  }

  isZero = a == 0;
  isInt = b == 1;
  isOne = isInt && a == 1;
  isMinusOne = isInt && a == -1;
}

bool MPZ_rational::isInteger() const { return isInt; }

bool MPZ_rational::isCOInteger() const { return (a == 1); }

MPZ_rational MPZ_rational::operator+(MPZ_rational const &n) const {
  return MPZ_rational(a * n.b + b * n.a, b * n.b);
}

MPZ_rational &MPZ_rational::operator+=(MPZ_rational const &n) {
  a = a * n.b + b * n.a;
  b = b * n.b;
  update();

  return *this;
}

MPZ_rational MPZ_rational::operator-(MPZ_rational const &n) const {
  return MPZ_rational(a * n.b - b * n.a, b * n.b);
}

void MPZ_rational::opp(MPZ_rational *&_c) const {
  _c = new MPZ_rational(-a, -b);
}

MPZ_rational &MPZ_rational::operator-=(MPZ_rational const &n) {
  a = a * n.b - b * n.a;
  b = b * n.b;
  update();

  return *this;
}

MPZ_rational MPZ_rational::operator*(MPZ_rational const &n) const {
  return MPZ_rational(a * n.a, b * n.b);
}

MPZ_rational &MPZ_rational::operator*=(MPZ_rational const &n) {
  a *= n.a;
  b *= n.b;
  update();

  return *this;
}

MPZ_rational MPZ_rational::operator/(MPZ_rational const &n) const {
  if (n.isZero)
    throw(0);

  return MPZ_rational(a * n.b, b * n.a);
}

MPZ_rational &MPZ_rational::operator/=(MPZ_rational const &n) {
  if (n.isZero)
    throw(0);

  a *= n.b;
  b *= n.a;
  update();

  return *this;
}

MPZ_rational &MPZ_rational::operator=(long int i) {
  a = i;
  b = 1;

  isZero = a == 0;
  isInt = b == 1;
  isOne = isInt && a == 1;
  isMinusOne = isInt && a == -1;

  return *this;
}

bool MPZ_rational::operator>=(const int &ni) const {
  MPZ_rational n(*this - ni);
  return (n.a >= 0);
}

bool MPZ_rational::operator==(const int &i) const { return (a == i && b == 1); }

bool MPZ_rational::operator==(MPZ_rational const &n) const {
  return (a == n.a && b == n.b);
}

bool MPZ_rational::operator==(const mpz_class &n) const {
  if (!isInt)
    return false;

  return (a == n);
}

bool MPZ_rational::operator!=(MPZ_rational const &n) const {
  return (a != n.a || b != n.b);
}

void MPZ_rational::print(ostream &o) const {
  if (isInt)
    o << a.get_str();
  else
    o << a.get_str() << "/" << b.get_str();
}

string MPZ_rational::to_string() const {
  if (isInt)
    return a.get_str();
  else
    return (a.get_str() + "/" + b.get_str());
}

ostream &operator<<(ostream &o, MPZ_rational const &n) {
  n.print(o);

  return o;
}

/*
BigInteger gcd(const BigInteger &n, const BigInteger &m)
{
        return (m == 0 ? n : gcd(m, n % m));
}*/
