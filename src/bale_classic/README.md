# bale (classic)

### Table of contents

* [What is bale?](#What-is-bale)
* [From the Book (formerly known as AGI)](#From-the-Book)
* [What is in bale?](#What-is-in-bale)
* [System Requirements](#System-Requirements)
* [Documentation](#Documentation)
* [Build Instructions](#Build-Instructions)
* [Run Instructions](#Running)
* [Testing](#Testing)
* [Version History](#Version-History)

### What is bale?

#### In one sentence
A collection of buffered communication libraries and some mini-applications.

#### The elevator pitch

The bale effort is, first and foremost, *a vehicle for discussion* for parallel programming productivity.  

The bale effort attempts to:

- demonstrate some challenges of implementing interesting (i.e. irregular) scalable distributed parallel applications.

- demonstrate an approach to achieving high performance for the internode communication in such applications

- explore concepts that make it easier to write, maintain, and get top performance from such applications


bale doesn’t claim to have the answers to better parallel programming. We use bale to evolve our thinking on parallel programming in the effort to make parallel programming easier and more productive and… fun! Yes, we think making it fun is a worthy goal. 

Note: We call this directory bale classic because this collection of code is the evolution of the original bale release. bale as a repository now contains other "cousins" of bale that are inspired by bale classic.

#### Pillars of bale: apps and aggregation

One pillar of bale is a directory of [apps](apps/README.md) that exhibit interesting communication patterns and programming demands. The apps can be written with aggregated communication as opposed to fine-grained point-to-point communication. We think aggregation is and will remain vital to getting top performance on parallel computers, but we don't like how difficult it is to write programs that use aggregation. Each app is written in multiple ways to demonstrate the pros and cons of each. 

The other pillar of bale are three libraries that provide the programmer with an API to aggregate communications within application code. These libraries are: exstack, exstack2, and conveyors. And while their API's are quite similar, they differ in their underlying implementations and behaviors.

### From the Book

One of the questions we ask ourselves in bale is what is the "best" version of an app? Obviously subjective, we consider ease of reading, ease of understanding what is happening when the code runs, ease of writing, and performance. We call this elusive version, "From the Book" or FTB in honor of [Paul Erdos](https://en.wikipedia.org/wiki/Proofs_from_THE_BOOK).

### What is in bale?

The main components are :

- [libgetput](libgetput/README.md)  - parallel programming allowing simple remote gets, puts, and atomics. libgetput is a library that can be compiled on top of UPC or SHMEM. Everything in bale except conveyors is built on top of libgetput.
- [exstack](exstack/README.md)   - The exstack and exstack2 libraries for aggregating communication
- [convey](convey/README.md) - The conveyor library for aggregating communication. More mature and sophisticated aggregation library than the exstacks.
- [spmat](spmat/README.md)  -  a sparse matrix library
- [std_options](std_options/README.md) - options parsing library
- [apps](apps/README.md)  -  the applications directory.

### System Requirements
Bale is written in C and can be compiled with UPC or linked against OpenSHMEM 1.4. Bale has been tested on a variety of architectures including Cray XC30, clusters with Infiniband, and SMP Linux. See [INSTALL.md](INSTALL.md) for detailed build instructions and [DEMO.md](DEMO.md) for some quick installation demos.

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

#### Other Dependencies

bale also depends on:

  - argp (part of glibc)

  - autoconf (2.69)

  - automake (1.13.4)

### Documentation

Bale is documented using Doxygen. To generate the html documentation you must have 
doxygen on your system (I recommend version 1.8.6). To generate the docs, simply change
to the $BALEDIR directory and type:
    `doxygen`
Then navigate your browser to $BALEDIR/html/index.html.

### Build Instructions
There is an install script (called install.sh) to make building easier for most people. Follow the directions in [INSTALL.md](INSTALL.md). 

### Running

All bale apps have a common set of standard options. In addition, bale apps that work on matrices or graphs have a common set graph input options. Run any app with '--help' for more information.

We have included a python script (run_apps.py) to make it easier to run a suite of tests in bale. Run this script with '--help' option for more information.

We have also included a Jupyter notebook to enable visualization and analysis of the results of bale runs. This notebook is in plot_results.ipynb.

### Testing
We have a unit test script that uses pytest as a harness. To run this test, go to
the apps directory and run

    pytest -s -P=<path/to/bale/binaries> --node_range=0,10,2 -M 15

### Version History

* May 2018: Initial Release v 1.0.0 
* Dec. 2018: bale 2.0.0 
  * New apps: transpose, randperm, permute_matrix, write_sparse_matrix

* Aug. 2018: bale 2.1.0 
  * update conveyors to version 0.5.0
* bale 3.0.0
  * Added cousins: Rust, Serial C, and Chapel to bale
  * New graph model: Geometric graphs
  * New app: SSSP
  * arg_parse replaced getopt
  * unit tests with pytest
  * docker files