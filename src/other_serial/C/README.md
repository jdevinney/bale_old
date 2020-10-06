> Copyright (c) 2020, Institute for Defense Analyses
> 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500\n
> All rights reserved.\n
> This file is part of Bale.   For licence information see the
> LICENSE file in the top level dirctory of the distribution.

# C_bale: the bale apps written in serial C serial
### In one sentence
A textbook, C, implementation of the apps in bale.

### The elevator pitch

The bale effort is, first and foremost, 
a vehicle for discussion for parallel programming productivity.  
We have included a simple C version of the apps in 
bale_classic as a concrete description of the apps written in a familiar language.
Some of the apps (like histo and ig) are trivial as serial apps.
It might be useful to view the serial version of the more complicated apps
(like toposort and sssp) before dealing with the parallelism and buffer communication 
of the bale_classic apps.  The C version of the way we implement a 
compressed row format data structure for a sparse matrix is also simpler than
the parallel version.  Finally, we also have examples of applications that
are efficient as serial codes, but have no known parallel implementations.

This is a self contained directory that does not depend 
on build process for the rest of bale.

## Nitty Gritty

### Where does it run?
This is written in generic C, hopefully it will run in any C environment.
It does depend on the ``argp`` library for command line argument parsing,
doxygen for documentation and ``pytest`` for unit testing.

### What is included here:

- Makefile - simple explicit makefile for all the apps.
- APPS:
  - [histo](histo.md) -- creates a large histogram (random stores) (see doxygen: histo.c)
  - [ig](ig.md) -- a large gather (random loads) (see doxygen: ig.c)
  - [randperm](randperm.md) -- creates a random permutation (see doxygen: randperm.c)
  - transpose_matrix.c -- computes the transpose of a sparse matrix [details](transpose_matrix.c)
  - [permute_matrix](prermute_matrix.md) -- applies row and column permutations to a sparse matrix (see doxygen: permute_matrix.c)
  - triangle.c -- counts the number triangles in a graph [details](triangle.md)
  - [toposort](toposort.md) -- performs a toposort (matrix) sort of a morally upper triangular matrix (see doxygen: toposort.c)
  - [sssp](sssp.md) -- solves the single source shortest path problem on a graph (see doxygen: sssp.c)
  - unionfind.c -- uses the union-find data structure to find connected components in a graph [details](unionfind.md)
- Other:
- spmat_utils.h, spmat_utils.c -- the sparse matrix library and some support functions [details](spmat_utils.md)
- std_options.h, std_options.c -- the command line parsing routines [details](std_options.md)

## Build Instructions
This is meant to be basic C, with a simple Makefile, so hopefully just typing 'make' will work.
```
    make  # make the apps
    make test # run the unit tests
    make doc # make the doxygen documentation
```
However, we do use the argp library from the GNU standard library. 
This is usually present by default on most linux systems. 
On Mac, you will probably have to install it by hand 
(one way to do this is 'brew install argp-standalone') and then mess with your
LD_LIBRARY_PATH and LDFLAGS variables.

## Testing
For bale 3.0 we have started using pytest for the unit testing.
One can run the test specified in file ``tests/test_all.py`` by hand with the command:

```
    pytest -s
```
For more details on how to modify the tests see [pytest](pytest.md)

## Documentation
serial_C is documented using Doxygen. 
If doxygen is installed on your system, you should be able to type ``make doc``
and then load ``html/index.html`` into a browser.

## Discussion



