# bale (classic)

### Table of contents

* [What is bale?](#What-is-bale)
* [What is in bale?](#What-is-in-bale)
* [System Requirements](#System-Requirements)
* [Documentation](#Documentation)
* [Build Instructions](#Build-Instructions)
* [Run Instructions](#Running)
* [Testing](#Testing)

### What is bale?

#### In one sentence

The bale effort is, first and foremost, *a vehicle for discussion* for parallel programming productivity.  

#### Some more detail

The bale effort attempts to:

- demonstrate some challenges of implementing interesting (i.e. irregular) scalable distributed parallel applications.

- demonstrate an approach to achieving high performance for the internode communication in such applications

- explore concepts that make it easier to write, maintain, and get top performance from such applications

bale doesn’t claim to have the answers to better parallel programming.

bale is not a collection of benchmarks.

We use bale to evolve our thinking on parallel programming in the effort to make parallel programming easier, more productive, and more fun. Yes, we think making it fun is a worthy goal!

Note: We call this directory bale classic because this collection of code is the evolution of the original bale release. bale as a repository now contains other "cousins" of bale that are inspired by bale classic.

#### The two pillars of bale: aggregation and apps

**Aggregation**

bale_classic contains three libraries that provide the programmer with an API to aggregate communications within application code. These libraries are: exstack, exstack2, and conveyors. And while their API's are quite similar, they differ in their underlying implementations and behaviors. We think aggregation is and will remain vital to getting top performance on parallel computers, but we don't like how difficult it is to write programs that use aggregation

**Apps**

bale_classic also contains a directory of [apps](apps/README.md) that exhibit interesting communication patterns and programming demands. The apps can be written with aggregated communication as opposed to fine-grained point-to-point communication. Each app is written in multiple ways to demonstrate the pros and cons of each. 

One of the questions we ask ourselves in bale is what is the "best" version of an app? Obviously subjective, we consider ease of reading, ease of understanding what is happening when the code runs, ease of writing, and performance. We call this elusive version, "**From the Book**" or FTB in honor of [Paul Erdos](https://en.wikipedia.org/wiki/Proofs_from_THE_BOOK).

### What is in bale?

The main components are :

- [libgetput](libgetput/README.md)  - parallel programming library allowing simple remote gets, puts, and atomics. libgetput can be compiled on top of UPC or SHMEM. Everything in bale except conveyors is built on top of libgetput.
- [exstack](exstack/README.md)   - The exstack and exstack2 libraries for aggregating communication
- [convey](convey/README.md) - The conveyor library for aggregating communication. More mature and sophisticated aggregation library than the exstacks.
- [spmat](spmat/README.md)  -  a sparse matrix library
- [std_options](std_options/README.md) - options parsing library
- [apps](apps/README.md)  -  the applications directory.

### System Requirements
bale_classic is written in C and can be compiled with UPC or linked against OpenSHMEM 1.4. bale_classic has been tested on a variety of architectures including Cray XC30, clusters with Infiniband, and SMP Linux. See [INSTALL.md](INSTALL.md) for detailed build instructions and [DEMO.md](DEMO.md) for some quick installation demos.

All of bale_classic is supported and has been tested on:

- Cray UPC (cce 8.7.3)
- GNU UPC (5.2.0.1)
- Clang UPC (3.9.1-1)
- OpenMPI (OSHMEM) 4.0.2 with UCX on infiniband
- Cray SHMEM (7.7.2)
- Cray openshmemX (8.0.1)

There are problems with the following...
- OpenMPI (OSHMEM) 4.0.3 with UCX on an SMP (progress issues with exstack2 and conveyors)

- Sandia OpenSHMEM (SOS) (1.4.5) with CMA on an SMP (conveyors)

#### Other Dependencies

bale_classic also depends on:

  - argp (part of glibc)

  - autoconf (2.69)

  - automake (1.13.4)

  - doxygen (optional)

  - python 3.5 (for unit tests and make_bale script, python 3.7 is required for run_apps.py script)

### Documentation

bale_classic is documented using Doxygen. To generate the html documentation you must have 
doxygen on your system (I recommend version 1.8.6). To generate the docs, simply change
to the $BALEDIR directory and type:
    `doxygen`
Then navigate your browser to $BALEDIR/html/index.html.

### Build Instructions
We include a build script (called make_bale) to make building easier for most people. Detailed instructions are found in [INSTALL.md](INSTALL.md). Also check out some quick start instructions in [DEMO.md](DEMO.md) 

### Running

All bale_classic apps have a common set of standard options. In addition, bale apps that work on matrices or graphs have a common set graph input options. Run any app with '--help' for more information.

We have included a python script (run_apps.py) to make it easier to run a suite of tests in bale. Run this script with '--help' option for more information.

We have also included a Jupyter notebook to enable visualization and analysis of the results of bale runs. This notebook is in plot_results.ipynb.

### Testing
We have a unit test script that uses pytest as a harness. To run this test, first build bale, then navigate to the apps directory and run

    pytest -s -P=<path/to/bale/binaries> --node_range=0,10,2 -M 15

If you get an error that claims that pytest cannot load the conftest.py file, I have found that deleting the
__pycache__ directories from the apps and apps/tests directories fixes this.

