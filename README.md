
# User Guide

There are two types of files used here. There are python files, from which polynomial coefficients are generated, and MATLAB scripts, which generate DC decompositions and run the convex-concave procedure.

## Python Files


The file polynomial.py contains an implementation for symbolic computation of multi variable polynomials which are addition, subtraction, multiplication, and exponentiation.
The file im.py generates the coefficients of the polynomial for the shape for shading problem.
Run the file as 

```
$ python im.py <path to image>.png
```
The image is expected to be greyscale. After running the above command, a folder with files {i}_{j}.txt will be generated.
Each file contains the coefficients for each f_r which is in binary order of exponents.

These coefficients are calculated according to Equation (12) in [1].

Requirements:

python==3.10.0
Pillow==10.0.0
numpy==1.25.2
tqdm==4.66.1

## MATLAB files

The MATLAB files should be compatible with any MATLAB version post- 2017b. There are a large number of files, and more details to follow soon.

### Requirements

Minimum Requirements are an LP solver and an SDP solver (for DD-Local and SDP-Local). Gurobi is the fastest, but requires an academic license. Linprog, SeDuMi, and CVX are alternatives, and we will provide versions supporting them shortly.


### Utilities

There are several utility functions, including:

* lex_index_nh.m
* lex_exps.m
* get_hesscoeffs.m

Full list to follow.

### DCD constructors

