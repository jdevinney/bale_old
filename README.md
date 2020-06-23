# bale
### In one sentence
A collection of buffered communication libraries and some mini-applications.

### The elevator pitch

The bale effort is, first and foremost, a vehicle for discussion for parallel programming productivity.  

One pillar of bale is a directory of "apps" that exhibit interesting communication patterns and programming demands. The apps can be written with aggregated communication as opposed to fine-grained point-to-point communication. We think aggregation is and will remain vital to getting top performance on parallel computers, but we don't like how difficult it is to write programs that use aggregation.

The other pillar of bale are three libraries that provide the programmer with an API to aggregate communications within application code. These libraries are: exstack, exstack2, and conveyors.

We hope that bale can lead to improved parallel programmer productivity (including existing and/or new programming models) and performance. 


## Nitty Gritty

### What language? Where does it run?
Bale is written in C and can be compiled with UPC or link against OpenSHMEM 1.4. Bale has been tested on a variety of architectures including Cray XC30, clusters with Infiniband, and SMP Linux.

The main components are :

- README.md  - this file
- INSTALL  - instructions for building bale
- libgetput  - parallel programming utility library. libgetput is a library that can be compiled on top of UPC or SHMEM and implements puts, gets, collectives, and atomics. Everything in bale except conveyors is built on top of libgetput.
- [exstack](exstack/README.md)   - The exstack and exstack2 libraries for aggregating communication
- convey - The conveyor library for aggregating communication
- spmat  -  a sparse matrix library
- apps  -  the applications directory. Includes histogram, 
   indexgather, toposort, triangle counting, randperm, etc.
- install.sh - the build and install script
- runall.sh - a test script that runs all applications
- Doxyfile - the main file for building the documentation with doxygen
- mainpage.h - The main documentation page.
- clang_upc_run.sh - a script to add "-fupc-threads-N" as appropriate when using clang-upc

All of Bale is supported and has been tested on:

- Cray UPC (cce 8.7.3)
- GNU UPC (5.2.0.1)
- Clang UPC (3.9.1-1)
- (not working yet) OpenMPI (OSHMEM) 3.0.6 
- (all but conveyors) Sandia OpenSHMEM (SOS) (1.4.5)
- Cray SHMEM (7.7.2)
- Cray openshmemX ?

Each of the applications in bale are implemented in multiple ways to showcase the pros and cons of each. In general the models used are:

1. agi      : standard PGAS model that uses puts, gets, and atomics
2. exstack  : a bulk synchronous buffering model
3. exstack2 : a asynchronous variant of exstack
4. conveyors: a more mature and sophisticated asynchronous model that is independent of exstack and exstack2.
5. In some applications we have included other variants in the "alternatives" directory.   

## Build Instructions
There is an install script (called install.sh) to make building easier for most people. Follow the directions in INSTALL. 

## Documentation
Bale is documented using Doxygen. See INSTALL for directions on how to generate the documentation.
