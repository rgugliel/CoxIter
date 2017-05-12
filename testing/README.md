This folder contains the necessary files to test CoxIter.

# General information
The file `tests.txt` contains data about a large set of hyperbolic Coxeter groups. CoxIter will be run for each of these graph, and the computed values will be compared to expected values.

# Building
In `build/` do:
```
cmake ../
make
```

# Use
* First, extract the archive `graphs/testing_index2.tar.gz`
* In `build/` do:
```
./tests
```

The tests should complete in a few minutes. A summary is displayed and details are saved in the file `tests.txt.output`


# The file tests.txt
Each line corresponds to one graph, with the following format:
* The first part is the path to the .coxiter file
* Euler characteristic: if a rational number is found, it will be assumed to be the Euler characteristic
* Compactness: "compact" or "non-compact"
* Cofiniteness: "non-fv" if non-cofininite (otherwise it is assumed that the group is cofinite)
* Arithmeticity: "arithmetic" if arithmetic 
* f-vector: given between parentheses, eg "(65,131,67,1)"
* Growth series: an example for the syntax is given by "f(x) = C(2,2,2,2,2,3,3,4,4,6,6,8,12)/(1 - x - x^2 + x^3 - 2 * x^5 + x^6 + x^7 + x^10", where the first coefficients denote Cyclotomic polynomial
* Growth rate: an exemple is given by "tau=1.1762808182599175065440703384740350507;"
