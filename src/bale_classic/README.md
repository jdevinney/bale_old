# bale (classic)
### In one sentence
A collection of buffered communication libraries and some mini-applications.

### The elevator pitch

The bale effort is, first and foremost, *a vehicle for discussion for parallel programming productivity*.  

One pillar of bale is a directory of [apps](apps/README.md) that exhibit interesting communication patterns and programming demands. The apps can be written with aggregated communication as opposed to fine-grained point-to-point communication. We think aggregation is and will remain vital to getting top performance on parallel computers, but we don't like how difficult it is to write programs that use aggregation. Each app is written in multiple ways to demonstrate the pros and cons of each.

The other pillar of bale are three libraries that provide the programmer with an API to aggregate communications within application code. These libraries are: exstack, exstack2, and conveyors. And while their API's are quite similar, they differ in their underlying implementations and behaviors.

We hope that bale can lead to improved parallel programmer productivity (including existing and/or new programming models) and performance. 

Note: We call this directory bale classic because this collection of code is the evolution of the original bale release. bale as a repository now contains lots of "flavors" of bale that are inspired by bale classic.


## Nitty Gritty of Bale Classic

### What language? Where does it run?
Bale is written in C and can be compiled with UPC or linked against OpenSHMEM 1.4. Bale has been tested on a variety of architectures including Cray XC30, clusters with Infiniband, and SMP Linux.

The main components are :

- README.md  - this file
- [INSTALL](INSTALL)  - instructions for building bale
- [libgetput](libgetput/README.md)  - parallel programming utility library. libgetput is a library that can be compiled on top of UPC or SHMEM and implements puts, gets, collectives, and atomics. Everything in bale except conveyors is built on top of libgetput.
- [exstack](exstack/README.md)   - The exstack and exstack2 libraries for aggregating communication
- [convey](convey/README.md) - The conveyor library for aggregating communication
- [spmat](spmat/README.md)  -  a sparse matrix library
- [apps](apps/README.md)  -  the applications directory. Includes [histogram](apps/histo_src/README.md), 
   [indexgather](apps/ig_src/README.md), [toposort](apps/topo_src/README.md), [transpose_matrix](apps/transpose_matrix_src/README.md), [triangle counting](apps/triangle_src/README.md), [randperm](apps/randperm_src/README.md), etc.
- install.sh - the build and install script
- runall.sh - a demo script that runs all applications
- Doxyfile - the main file for building the documentation with doxygen
- mainpage.h - The main documentation page.
- clang_upc_run.sh - a script to add "-fupc-threads-N" as appropriate when using clang-upc

All of Bale is supported and has been tested on:

- Cray UPC (cce 8.7.3)
- GNU UPC (5.2.0.1)
- Clang UPC (3.9.1-1)
- OpenMPI (OSHMEM) 3.0.6 with UCX on infiniband
- Cray SHMEM (7.7.2)
- Cray openshmemX ?

There are problems with the following...
- OpenMPI (OSHMEM) 3.0.6 with UCX on an SMP (progress issues)
- Sandia OpenSHMEM (SOS) (1.4.5) without MCA (all but conveyors)

### Build Instructions
There is an install script (called install.sh) to make building easier for most people. Follow the directions in INSTALL. 

### Testing
We are working on a new test script that uses pytest as a harness. To run this test, go to
the apps directory and run

    pytest -s -P=<path/to/bale/binaries> --node_range=0,10,2 -M 15

### Documentation
Bale is documented using Doxygen. See INSTALL for directions on how to generate the documentation.

### 
