This folder contains the necessary files to test CoxIter.

# General information
The file `tests.txt` contains data about a large set of hyperbolic Coxeter groups. CoxIter will be run for each of these graph, and the computed values will be compared to expected values.

# Building
In build/, do:
```
cmake ../
make
```

# Use

# The file tests.txt
Each line corresponds to one graph, with the following format:
* The first part is the path to the .coxiter file
* Euler characteristic: if a rational number is found, it will be assumed to be the Euler characteristic
* Compactness: "compact" or "non-compact"
* Cofiniteness: "non-fv" if non-cofininite (otherwise it is assumed that the group is cofinite)
* Arithmeticity: "arithmetic" if arithmetic 
* f-vector: given between parentheses, eg "(65,131,67,1)"
