# bale / serial_C
### In one sentence
A *textbook*, C, implementation of the apps in bale.

### The elevator pitch

The bale effort is, first and foremost, a vehicle for discussion for parallel programming productivity.  We have included a simple C version of the apps as a concrete description of the apps written in a familiar language.

This is a self contained directory that does not depend on build process for the rest of bale.

## Nitty Gritty

### Where does it run?
We have been using gcc, but this should run in C environment.

### What is included here:

- README.md  - this file
- Makefile   - simple explicit makefile for all the apps.
- spmat_utils.c/h - the only library or support code.
   It contains the sparse matrix library and few helpful routines for debugging and timing.
- demo_spmat.c - essentially a script used to test the development of spat_utils. It might be good way to familiarize onesself we our version of a compressed row format for a sparse matrix.
- apps:
  - histo.c            -- creates a large histogram (random stores)
  - ig.c               -- a large gather (random loads)
  - randperm.c         -- creates a random permutation
  - transpose_matrix.c -- part of spmat_utils, but interesting in it own right.
  - permute_matrix.c   -- part of spmat_utils, but interesting in it own right.
  - triangle.c         -- counts triangles in a graph
  - toposort.c         -- performs a toposort (matrix) sort on a morally upper triangular matrix
  - sssp.c             -- single source shortest path in a graph
  - unionfind.c        -- uses the union-find data structure to find connected components in a graph
- runall.sh - a test script that runs all applications on default parameters
- Doxyfile - the main file for building the documentation with doxygen
- mainpage.h - The main documentation page.

## Build Instructions
Should only need to type make.

## Testing
We are working on a new test script that uses pytest as a harness. To run this test, go to
the apps directory and run
    pytest -s -P=<path/to/bale/binaries> --node_range=0,10,2 -M 15

## Documentation
serial_C is documented using Doxygen. 
Should be able to type make doc

## Discussion

