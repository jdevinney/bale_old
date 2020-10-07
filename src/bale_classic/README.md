# bale (classic)

### Table of contents

* [What is bale?](#What-is-bale)
* [What is in bale?](#What-is-in-bale)
* [System Requirements](#System-Requirements)
* [Documentation](#Documentation)
* [Build Instructions](#Build-Instructions)
* [Run Instructions](#Running)
* [Testing](#Testing)
* [Version History](#Version-History)
* [Contact](#Contact)

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

bale contains three libraries that provide the programmer with an API to aggregate communications within application code. These libraries are: exstack, exstack2, and conveyors. And while their API's are quite similar, they differ in their underlying implementations and behaviors. We think aggregation is and will remain vital to getting top performance on parallel computers, but we don't like how difficult it is to write programs that use aggregation

**Apps**

bale also contains a directory of [apps](apps/README.md) that exhibit interesting communication patterns and programming demands. The apps can be written with aggregated communication as opposed to fine-grained point-to-point communication. Each app is written in multiple ways to demonstrate the pros and cons of each. 

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
Bale is written in C and can be compiled with UPC or linked against OpenSHMEM 1.4. Bale has been tested on a variety of architectures including Cray XC30, clusters with Infiniband, and SMP Linux. See [INSTALL.md](INSTALL.md) for detailed build instructions and [DEMO.md](DEMO.md) for some quick installation demos.

All of Bale is supported and has been tested on:

- Cray UPC (cce 8.7.3)
- GNU UPC (5.2.0.1)
- Clang UPC (3.9.1-1)
- OpenMPI (OSHMEM) 4.0.2 with UCX on infiniband
- Cray SHMEM (7.7.2)
- Cray openshmemX ?

There are problems with the following...
- OpenMPI (OSHMEM) 4.0.3 with UCX on an SMP (progress issues)

- Sandia OpenSHMEM (SOS) (1.4.5) without MCA (all but conveyors)

#### Other Dependencies

bale also depends on:

  - argp (part of glibc)

  - autoconf (2.69)

  - automake (1.13.4)

  - doxygen (optional)

  - python 3.5 (for unit tests and make_bale script, python 3.7 is required for run_apps.py script)

#### Docker files
bale now comes with some Docker files to assist you in getting a bale-friendly environment on desktop linux. These files are in the *docker* directory. There are 4 sub-directories here:
- cupc (Clang UPC)
- gupc (GNU UPC)
- oshmem (OpenMPI OpenSHMEM)
- sos (Sandia OpenSHMEM)

Each directory has a file "Dockerfile" that when run will build an evironment that is capable of building all of bale and running bale apps.

### Documentation

Bale is documented using Doxygen. To generate the html documentation you must have 
doxygen on your system (I recommend version 1.8.6). To generate the docs, simply change
to the $BALEDIR directory and type:
    `doxygen`
Then navigate your browser to $BALEDIR/html/index.html.

### Build Instructions
There is an build script (called make_bale) to make building easier for most people. Follow the directions in [INSTALL.md](INSTALL.md). 

### Running

All bale apps have a common set of standard options. In addition, bale apps that work on matrices or graphs have a common set graph input options. Run any app with '--help' for more information.

We have included a python script (run_apps.py) to make it easier to run a suite of tests in bale. Run this script with '--help' option for more information.

We have also included a Jupyter notebook to enable visualization and analysis of the results of bale runs. This notebook is in plot_results.ipynb.

### Testing
We have a unit test script that uses pytest as a harness. To run this test, first build bale, then navigate to the apps directory and run

    pytest -s -P=<path/to/bale/binaries> --node_range=0,10,2 -M 15

### Version History

* May 2018: Initial Release version 1.0.0 

* Dec. 2018: version 2.0.0 
  * New apps: transpose, randperm, permute_matrix, write_sparse_matrix

* Aug. 2018: version 2.1.0 
  * update conveyors to version 0.5.0

* Nov. 2020: version 3.0.0
  * Added cousins: Rust, Serial C, and Chapel to bale
  * New graph model: Geometric graphs
  * New app: SSSP
  * replaced write_sparse_matrix with sparse_matrix_io
  * arg_parse replaced getopt
  * unit tests with pytest
  * new make_bale script
  * new run_apps script
  * docker files
  * AGP (Atomics, Gets, and Puts) replaces AGI: simple PGAS style programming
  * FTB (From the Book) replaces As God Intended.

### Contact

bale@super.org