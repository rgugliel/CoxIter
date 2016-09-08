CC = g++
CXX = g++
CXXFLAGS= -std=c++11 -Wall -Wno-reorder -Wno-unknown-pragmas -pedantic -O3 -fopenmp -Ic:/pcre/include -Ic:/gmp/include
LDFLAGS = -LC:/gmp/lib -LC:/pcre/lib -lgmpxx -lgmp -lpcre
SRCPREFIX=..

all: default

default: main.o graph.o coxiter.o mpz_rational.o number_template.o string.o graphs.list.o graphs.list.n.o graph.o graphs.product.set.o graphs.product.o graphs.list.iterator.o polynomials.o regexp.o app.o gbd.o arithmeticity.o
	$(CXX) $(CXXFLAGS) main.o coxiter.o mpz_rational.o number_template.o string.o graphs.list.o graphs.list.n.o graph.o graphs.product.set.o graphs.product.o graphs.list.iterator.o polynomials.o regexp.o app.o gbd.o arithmeticity.o -o coxiter $(LDFLAGS)
	
main.o : $(SRCPREFIX)/main.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@
	
coxiter.o : $(SRCPREFIX)/coxiter.cpp $(SRCPREFIX)/coxiter.h $(SRCPREFIX)/graphs.list.iterator.h $(SRCPREFIX)/graphs.list.h $(SRCPREFIX)/graphs.product.set.h $(SRCPREFIX)/graphs.product.h $(SRCPREFIX)/lib/string.h $(SRCPREFIX)/lib/numbers/mpz_rational.h $(SRCPREFIX)/lib/polynomials.h
	$(CXX) $(CXXFLAGS) -c $< -o $@
	
app.o : $(SRCPREFIX)/app.cpp $(SRCPREFIX)/app.h $(SRCPREFIX)/coxiter.h $(SRCPREFIX)/gbd.h $(SRCPREFIX)/arithmeticity.h
	$(CXX) $(CXXFLAGS) -c $< -o $@
	
gbd.o : $(SRCPREFIX)/gbd.cpp $(SRCPREFIX)/gbd.h $(SRCPREFIX)/coxiter.h
	$(CXX) $(CXXFLAGS) -c $< -o $@
	
arithmeticity.o : $(SRCPREFIX)/arithmeticity.cpp $(SRCPREFIX)/arithmeticity.h $(SRCPREFIX)/coxiter.h
	$(CXX) $(CXXFLAGS) -c $< -o $@
	
# ----------------------------------------------------
# Graphs
graph.o : $(SRCPREFIX)/graph.cpp $(SRCPREFIX)/graph.h
	$(CXX) $(CXXFLAGS) -c $< -o $@
	
graphs.list.o : $(SRCPREFIX)/graphs.list.cpp $(SRCPREFIX)/graphs.list.h $(SRCPREFIX)/graphs.list.n.h 
	$(CXX) $(CXXFLAGS) -c $< -o $@
	
graphs.list.iterator.o : $(SRCPREFIX)/graphs.list.iterator.cpp $(SRCPREFIX)/graphs.list.iterator.h $(SRCPREFIX)/graphs.list.h 
	$(CXX) $(CXXFLAGS) -c $< -o $@
	
graphs.list.n.o : $(SRCPREFIX)/graphs.list.n.cpp $(SRCPREFIX)/graphs.list.n.h $(SRCPREFIX)/graph.h
	$(CXX) $(CXXFLAGS) -c $< -o $@
	
graphs.product.o : $(SRCPREFIX)/graphs.product.cpp $(SRCPREFIX)/graphs.product.h $(SRCPREFIX)/graph.h
	$(CXX) $(CXXFLAGS) -c $< -o $@
	
graphs.product.set.o : $(SRCPREFIX)/graphs.product.set.cpp $(SRCPREFIX)/graphs.product.set.h $(SRCPREFIX)/graph.h $(SRCPREFIX)/graphs.product.h
	$(CXX) $(CXXFLAGS) -c $< -o $@
	
# ----------------------------------------------------
# Other libs	
string.o : $(SRCPREFIX)/lib/string.cpp $(SRCPREFIX)/lib/string.h
	$(CXX) $(CXXFLAGS) -c $< -o $@
	
regexp.o : $(SRCPREFIX)/lib/regexp.cpp $(SRCPREFIX)/lib/regexp.h
	$(CXX) $(CXXFLAGS) -c $< -o $@
	
# ----------------------------------------------------
# Maths
mpz_rational.o : $(SRCPREFIX)/lib/numbers/mpz_rational.cpp $(SRCPREFIX)/lib/numbers/mpz_rational.h $(SRCPREFIX)/lib/string.h $(SRCPREFIX)/lib/numbers/number_template.h
	$(CXX) $(CXXFLAGS) -c $< -o $@
	
number_template.o : $(SRCPREFIX)/lib/numbers/number_template.cpp $(SRCPREFIX)/lib/numbers/number_template.h
	$(CXX) $(CXXFLAGS) -c $< -o $@
	
polynomials.o : $(SRCPREFIX)/lib/polynomials.cpp $(SRCPREFIX)/lib/polynomials.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

# ----------------------------------------------------
# Varia	
clean:
	rm *.o

distclean: clean
	rm -rf Makefile

doc:
	gap -A -q -T < makedocrel.g

.PHONY: all default clean distclean doc
